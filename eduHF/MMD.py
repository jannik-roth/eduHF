import numpy as np
from scipy.special import erf

from .basisfunction import ContractedGaussianFunction
from .geometry import Atom

def boys(T : float, n : int):
    if np.abs(T) < 1e-11:
        return 1.0 / (2.0 * n + 1.0)
    elif (n == 0):
        return 0.5 * np.sqrt(np.pi / T) * erf(np.sqrt(T))
    else:
        work = np.zeros(n+1)
        work[0] = 0.5 * np.sqrt(np.pi / T) * erf(np.sqrt(T))
        for i in range(1,n+1):
            work[i] = ((2*i-1) * work[i-1] - np.exp(-T)) / (2.0 * T)
        return work[n]

def E(n : int, m : int, T : int,
      AB : float, alpha : float, beta : float):
    
    p = alpha + beta
    q = alpha * beta / p

    if ((T < 0) or (T > (n+m))):
        return 0
    elif ((T == 0) and (n == 0) and (m == 0)):
        return np.exp(-q * AB ** 2)
    elif (n == 0):
        return (1.0 / (2.0 * p) * E(n, m-1, T-1, AB, alpha, beta)
                + alpha * AB / p * E(n, m-1, T, AB, alpha, beta)
                + (T+1.0) * E(n, m-1, T+1, AB, alpha, beta))
    else:
        return (1.0 / (2.0 * p) * E(n-1, m, T-1, AB, alpha, beta)
                - beta * AB / p * E(n-1, m, T, AB, alpha, beta)
                + (T+1.0) * E(n-1, m, T+1, AB, alpha, beta))
    
def E_der(n : int, m : int, T : int,
            AB : float, alpha : float, beta : float,
            center : int):
    res = 0.0
    
    if (center == 0):
        res = 2.0 * alpha * E(n+1, m, T, AB, alpha, beta)
        if (n > 0):
            res += -n * E(n-1, m, T, AB, alpha, beta)
    elif (center == 1):
        res = 2.0 * beta * E(n, m+1, T, AB, alpha, beta)
        if (m > 0):
            res += - m * E(n, m-1, T, AB, alpha, beta)
    return res

def R(t : int, u : int, v : int, n : int,
      p : float, xyz_p : np.array, xyz_c : np.array):
    
    diff_pc = xyz_p - xyz_c
    R_pc2 = np.sum(np.square(diff_pc))
    # print("R_pc2", R_pc2, xyz_p, xyz_c)
    tmp = 0.0
    
    if ((t == 0) and (u == 0) and (v == 0)):
        return (-2.0 * p) ** n * boys(p * R_pc2, n)
    elif ((t == 0) and (u == 0)):
        if (v > 1):
            tmp = (v - 1) * R(t, u, v-2, n+1, p, xyz_p, xyz_c)
        return tmp + diff_pc[2] * R(t, u, v-1, n+1, p, xyz_p, xyz_c)
    elif (t == 0):
        if (u > 1):
            tmp = (u - 1) * R(t, u-2, v, n+1, p, xyz_p, xyz_c)
        return tmp + diff_pc[1] * R(t, u-1, v, n+1, p, xyz_p, xyz_c)
    else:
        if (t > 1):
            tmp = (t - 1) * R(t-2, u, v, n+1, p, xyz_p, xyz_c)
        return tmp + diff_pc[0] * R(t-1, u, v, n+1, p, xyz_p, xyz_c)
    
def _overlap(alpha, L1, xyz1, beta, L2, xyz2):
    return (E(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta)
            * E(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta)
            * E(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta)
            * (np.pi / (alpha + beta)) ** 1.5)

def overlap_item(bf1 : ContractedGaussianFunction,
                 bf2 : ContractedGaussianFunction):
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            res += c1 * c2 * _overlap(a1, bf1.l_vec, bf1.xyz, 
                                      a2, bf2.l_vec, bf2.xyz)
    return res

def kinetic_item(bf1 : ContractedGaussianFunction,
                 bf2 : ContractedGaussianFunction):
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            res += c1 * c2 * (a2 * (2 * np.sum(bf2.l_vec) + 3) * _overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec, bf2.xyz)
                              - 2.0 * a2 ** 2 * (_overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec+np.array([2,0,0]), bf2.xyz)
                                                 + _overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec+np.array([0,2,0]), bf2.xyz)
                                                 + _overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec+np.array([0,0,2]), bf2.xyz))
                              - 0.5 * (bf2.l_vec[0] * (bf2.l_vec[0]-1) * _overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec+np.array([-2,0,0]), bf2.xyz)
                                       + bf2.l_vec[1] * (bf2.l_vec[1]-1) * _overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec+np.array([0,-2,0]), bf2.xyz)
                                       + bf2.l_vec[2] * (bf2.l_vec[2]-1) * _overlap(a1, bf1.l_vec, bf1.xyz, a2, bf2.l_vec+np.array([0,0,-2]), bf2.xyz)))
    return res

def _nuc_attraction(alpha : float, L1 : np.array, xyz1 : np.array,
                    beta : float, L2 : np.array, xyz2 : np.array,
                    xyza : np.array) -> float:
    p = alpha + beta
    xyzp = (xyz1 * alpha + xyz2 * beta) / p

    res = 0.0

    for t in range(L1[0]+L2[0]+1):
        for u in range(L1[1]+L2[1]+1):
            for v in range(L1[2]+L2[2]+1):
                res += (E(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta)
                        * E(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta)
                        * E(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta)
                        * R(t, u, v, 0, p, xyzp, xyza))
    return res * 2.0 * np.pi / p

def potential_1e_item(bf1 : ContractedGaussianFunction,
                      bf2 : ContractedGaussianFunction,
                      at : Atom) -> float:
    res = 0.0
    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            res += c1 * c2 * _nuc_attraction(a1, bf1.l_vec, bf1.xyz,
                                             a2, bf2.l_vec, bf2.xyz,
                                             at.xyz)
    return res

def _electron_repulsion(alpha : float, L1 : np.array, xyz1 : np.array,
                        beta : float, L2 : np.array, xyz2 : np.array,
                        gamma : float, L3 : np.array, xyz3 : np.array,
                        delta : float, L4 : np.array, xyz4 : np.array):
    
    p = alpha + beta
    q = gamma + delta
    xyzp = (alpha * xyz1 + beta * xyz2) / p 
    xyzq = (gamma * xyz3 + delta * xyz4) / q

    res = 0.0

    for t in range(L1[0]+L2[0]+1):
        for u in range(L1[1]+L2[1]+1):
            for v in range(L1[2]+L2[2]+1):
                for tau in range(L3[0]+L4[0]+1):
                    for nu in range(L3[1]+L4[1]+1):
                        for phi in range(L3[2]+L4[2]+1):
                            res += (E(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta)
                                    * E(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta)
                                    * E(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta)
                                    * E(L3[0], L4[0], tau, xyz3[0]-xyz4[0], gamma, delta)
                                    * E(L3[1], L4[1], nu, xyz3[1]-xyz4[1], gamma, delta)
                                    * E(L3[2], L4[2], phi, xyz3[2]-xyz4[2], gamma, delta)
                                    * R(t+tau, u+nu, v+phi, 0, p*q/(p+q), xyzp, xyzq)
                                    * (-1) ** (tau + nu + phi))
    
    return res * 2.0 * np.pi ** 2.5 / (p * q * np.sqrt(p + q))

def potential_2e_item(bf1 : ContractedGaussianFunction,
                      bf2 : ContractedGaussianFunction,
                      bf3 : ContractedGaussianFunction,
                      bf4 : ContractedGaussianFunction):
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            for c3, a3 in zip(bf3.coeffs, bf3.alphas):
                for c4, a4 in zip(bf4.coeffs, bf4.alphas):
                    res += c1 * c2 * c3 * c4 * _electron_repulsion(a1, bf1.l_vec, bf1.xyz,
                                                                   a2, bf2.l_vec, bf2.xyz,
                                                                   a3, bf3.l_vec, bf3.xyz,
                                                                   a4, bf4.l_vec, bf4.xyz)

    return res

def _overlap_der(alpha : float, L1 : np.array, xyz1 : np.array,
                 beta : float, L2 : np.array, xyz2 : np.array,
                 center : int, dim : int) -> float:

    if (dim == 0):
        return (E_der(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta, center)
                * E(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta)
                * E(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta)
                * (np.pi / (alpha + beta)) ** 1.50)
    elif (dim == 1):
        return (E(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta)
                * E_der(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta, center)
                * E(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta)
                * (np.pi / (alpha + beta)) ** 1.50)
    elif (dim == 2):
        return (E(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta)
                * E(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta)
                * E_der(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta, center)
                * (np.pi / (alpha + beta)) ** 1.50)

def overlap_der_item(bf1 : ContractedGaussianFunction,
                     bf2 : ContractedGaussianFunction,
                     center : int, dim : int) -> float:
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            res += c1 * c2 * _overlap_der(a1, bf1.l_vec, bf1.xyz,
                                          a2, bf2.l_vec, bf2.xyz,
                                          center, dim)
    return res

def _kinetic_der(alpha : float, L1 : np.array, xyz1 : np.array,
                 beta : float, L2 : np.array, xyz2 : np.array,
                 center : int, dim : int):
    
    if (dim == 0):
        kin1 = (beta * (2 * L2[0] + 1) * E_der(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta, center)
                - 2.0 * beta ** 2 * E_der(L1[0], L2[0]+2, 0, xyz1[0]-xyz2[0], alpha, beta, center)
                - 0.5 * L2[0] * (L2[0] - 1) * E_der(L1[0], L2[0]-2, 0, xyz1[0]-xyz2[0], alpha, beta, center))
    else:
        kin1 = (beta * (2 * L2[0] + 1) * E(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta)
                - 2.0 * beta ** 2 * E(L1[0], L2[0]+2, 0, xyz1[0]-xyz2[0], alpha, beta)
                - 0.5 * L2[0] * (L2[0] - 1) * E(L1[0], L2[0]-2, 0, xyz1[0]-xyz2[0], alpha, beta))
    if (dim == 1):
        kin2 = (beta * (2 * L2[1] + 1) * E_der(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta, center)
                - 2.0 * beta ** 2 * E_der(L1[1], L2[1]+2, 0, xyz1[1]-xyz2[1], alpha, beta, center)
                - 0.5 * L2[1] * (L2[1] - 1) * E_der(L1[1], L2[1]-2, 0, xyz1[1]-xyz2[1], alpha, beta, center))
    else:
        kin2 = (beta * (2 * L2[1] + 1) * E(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta)
                - 2.0 * beta ** 2 * E(L1[1], L2[1]+2, 0, xyz1[1]-xyz2[1], alpha, beta)
                - 0.5 * L2[1] * (L2[1] - 1) * E(L1[1], L2[1]-2, 0, xyz1[1]-xyz2[1], alpha, beta))
    if (dim == 2):
        kin3 = (beta * (2 * L2[2] + 1) * E_der(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta, center)
                - 2.0 * beta ** 2 * E_der(L1[2], L2[2]+2, 0, xyz1[2]-xyz2[2], alpha, beta, center)
                - 0.5 * L2[2] * (L2[2] - 1) * E_der(L1[2], L2[2]-2, 0, xyz1[2]-xyz2[2], alpha, beta, center))
    else:
        kin3 = (beta * (2 * L2[2] + 1) * E(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta)
                - 2.0 * beta ** 2 * E(L1[2], L2[2]+2, 0, xyz1[2]-xyz2[2], alpha, beta)
                - 0.5 * L2[2] * (L2[2] - 1) * E(L1[2], L2[2]-2, 0, xyz1[2]-xyz2[2], alpha, beta))
        
    if (dim == 0):
        tmp = E_der(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta, center)
        kin2 *= tmp
        kin3 *= tmp
    else:
        tmp = E(L1[0], L2[0], 0, xyz1[0]-xyz2[0], alpha, beta)
        kin2 *= tmp
        kin3 *= tmp
    if (dim == 1):
        tmp = E_der(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta, center)
        kin1 *= tmp
        kin3 *= tmp
    else:
        tmp = E(L1[1], L2[1], 0, xyz1[1]-xyz2[1], alpha, beta)
        kin1 *= tmp
        kin3 *= tmp
    if (dim == 2):
        tmp = E_der(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta, center)
        kin1 *= tmp
        kin2 *= tmp
    else:
        tmp = E(L1[2], L2[2], 0, xyz1[2]-xyz2[2], alpha, beta)
        kin1 *= tmp
        kin2 *= tmp

    return (kin1 + kin2 + kin3) * (np.pi / (alpha + beta)) ** 1.50

def kinetic_der_item(bf1 : ContractedGaussianFunction,
                     bf2 : ContractedGaussianFunction,
                     center : int, dim : int) -> float:
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            res += c1 * c2 * _kinetic_der(a1, bf1.l_vec, bf1.xyz,
                                          a2, bf2.l_vec, bf2.xyz,
                                          center, dim)
    
    return res

def _nuc_attraction_der(alpha : float, L1 : np.array, xyz1 : np.array,
                        beta : float, L2 : np.array, xyz2 : np.array,
                        xyza : np.array, dim : int) -> float:
    # center is always = 0
    p = alpha + beta
    xyzp = (xyz1 * alpha + xyz2 * beta) / p

    res = 0.0

    if (dim == 0):
        for t in range(L1[0]+L2[0]+1+1+1):
            for u in range(L1[1]+L2[1]+1+1):
                for v in range(L1[2]+L2[2]+1+1):
                    res += (E_der(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta, 0)
                            * E(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta)
                            * E(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta)
                            * R(t, u, v, 0, p, xyzp, xyza))
    elif (dim == 1):
        for t in range(L1[0]+L2[0]+1+1):
            for u in range(L1[1]+L2[1]+1+1+1):
                for v in range(L1[2]+L2[2]+1+1):
                    res += (E(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta)
                            * E_der(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta, 0)
                            * E(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta)
                            * R(t, u, v, 0, p, xyzp, xyza))
    elif (dim == 2):
        for t in range(L1[0]+L2[0]+1+1):
            for u in range(L1[1]+L2[1]+1+1):
                for v in range(L1[2]+L2[2]+1+1+1):
                    res += (E(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta)
                            * E(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta)
                            * E_der(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta, 0)
                            * R(t, u, v, 0, p, xyzp, xyza))
                    
    return res * 2.0 * np.pi / p
                    
def potential_1e_der_item(bf1 : ContractedGaussianFunction,
                          bf2 : ContractedGaussianFunction,
                          at : Atom, center : int, dim : int) -> float:
    if (center == 0):
        return potential_1e_der_item_un(bf1, bf2, at, dim)
    elif (center == 1):
        return potential_1e_der_item_un(bf2, bf1, at, dim) 

def potential_1e_der_item_un(bf1 : ContractedGaussianFunction,
                            bf2 : ContractedGaussianFunction,
                            at : Atom, dim : int) -> float:
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2. alphas):
            res += c1 * c2 * _nuc_attraction_der(a1, bf1.l_vec, bf1.xyz,
                                                 a2, bf2.l_vec, bf2.xyz,
                                                 at.xyz, dim)
    return res

def _electron_repulsion_der(alpha : float, L1 : np.array, xyz1 : np.array,
                               beta : float, L2 : np.array, xyz2 : np.array,
                               gamma : float, L3 : np.array, xyz3 : np.array,
                               delta : float, L4 : np.array, xyz4 : np.array,
                               dim : int) -> float:
    # center is always = 0
    
    p = alpha + beta
    q = gamma + delta
    xyzp = (alpha * xyz1 + beta * xyz2) / p 
    xyzq = (gamma * xyz3 + delta * xyz4) / q

    res = 0.0
    if (dim == 0):
        for t in range(L1[0]+L2[0]+1+1):
            for u in range(L1[1]+L2[1]+1):
                for v in range(L1[2]+L2[2]+1):
                    for tau in range(L3[0]+L4[0]+1):
                        for nu in range(L3[1]+L4[1]+1):
                            for phi in range(L3[2]+L4[2]+1):
                                res += (E_der(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta, 0)
                                        * E(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta)
                                        * E(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta)
                                        * E(L3[0], L4[0], tau, xyz3[0]-xyz4[0], gamma, delta)
                                        * E(L3[1], L4[1], nu, xyz3[1]-xyz4[1], gamma, delta)
                                        * E(L3[2], L4[2], phi, xyz3[2]-xyz4[2], gamma, delta)
                                        * R(t+tau, u+nu, v+phi, 0, p*q/(p+q), xyzp, xyzq)
                                        * (-1) ** (tau + nu + phi))
    if (dim == 1):
        for t in range(L1[0]+L2[0]+1):
            for u in range(L1[1]+L2[1]+1+1):
                for v in range(L1[2]+L2[2]+1):
                    for tau in range(L3[0]+L4[0]+1):
                        for nu in range(L3[1]+L4[1]+1):
                            for phi in range(L3[2]+L4[2]+1):
                                res += (E(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta)
                                        * E_der(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta, 0)
                                        * E(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta)
                                        * E(L3[0], L4[0], tau, xyz3[0]-xyz4[0], gamma, delta)
                                        * E(L3[1], L4[1], nu, xyz3[1]-xyz4[1], gamma, delta)
                                        * E(L3[2], L4[2], phi, xyz3[2]-xyz4[2], gamma, delta)
                                        * R(t+tau, u+nu, v+phi, 0, p*q/(p+q), xyzp, xyzq)
                                        * (-1) ** (tau + nu + phi))
    if (dim == 2):
        for t in range(L1[0]+L2[0]+1):
            for u in range(L1[1]+L2[1]+1):
                for v in range(L1[2]+L2[2]+1+1):
                    for tau in range(L3[0]+L4[0]+1):
                        for nu in range(L3[1]+L4[1]+1):
                            for phi in range(L3[2]+L4[2]+1):
                                res += (E(L1[0], L2[0], t, xyz1[0]-xyz2[0], alpha, beta)
                                        * E(L1[1], L2[1], u, xyz1[1]-xyz2[1], alpha, beta)
                                        * E_der(L1[2], L2[2], v, xyz1[2]-xyz2[2], alpha, beta, 0)
                                        * E(L3[0], L4[0], tau, xyz3[0]-xyz4[0], gamma, delta)
                                        * E(L3[1], L4[1], nu, xyz3[1]-xyz4[1], gamma, delta)
                                        * E(L3[2], L4[2], phi, xyz3[2]-xyz4[2], gamma, delta)
                                        * R(t+tau, u+nu, v+phi, 0, p*q/(p+q), xyzp, xyzq)
                                        * (-1) ** (tau + nu + phi))
    
    return res * 2.0 * np.pi ** 2.5 / (p * q * np.sqrt(p + q))

def potential_2e_der_item_un(bf1 : ContractedGaussianFunction,
                             bf2 : ContractedGaussianFunction,
                             bf3 : ContractedGaussianFunction,
                             bf4 : ContractedGaussianFunction,
                             dim : int) -> float:
    
    res = 0.0

    for c1, a1 in zip(bf1.coeffs, bf1.alphas):
        for c2, a2 in zip(bf2.coeffs, bf2.alphas):
            for c3, a3 in zip(bf3.coeffs, bf3.alphas):
                for c4, a4 in zip(bf4.coeffs, bf4.alphas):
                    res += c1 * c2 * c3 * c4 *_electron_repulsion_der(a1, bf1.l_vec, bf1.xyz,
                                                                      a2, bf2.l_vec, bf2.xyz,
                                                                      a3, bf3.l_vec, bf3.xyz,
                                                                      a4, bf4.l_vec, bf4.xyz,
                                                                      dim)
    return res

def potential_2e_der_item(bf1 : ContractedGaussianFunction,
                          bf2 : ContractedGaussianFunction,
                          bf3 : ContractedGaussianFunction,
                          bf4 : ContractedGaussianFunction,
                          center : int, dim : int) -> float:
    if (center == 0):
        return potential_2e_der_item_un(bf1, bf2, bf3, bf4, dim)
    elif (center == 1):
        return potential_2e_der_item_un(bf2, bf1, bf4, bf3, dim)
    elif (center == 2):
        return potential_2e_der_item_un(bf3, bf4, bf1, bf2, dim)
    elif (center == 3):
        return potential_2e_der_item_un(bf4, bf3, bf1, bf2, dim)
