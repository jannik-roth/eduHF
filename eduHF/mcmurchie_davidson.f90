function boys(T, n)
    implicit none
    real*8 :: boys
    real*8, intent(in) :: T
    integer, intent(in) :: n

    real*8 :: work(n+1)
    integer :: i
    real*8 :: PI = 3.14159265358979323846264338327950288419

    if ( abs(T) .le. 1.0e-11 ) then 
        boys = 1.0 / (2.0 * n + 1)
    else if ( n == 0 ) then
        boys = 0.5 * sqrt(PI/T) * erf(sqrt(T))
    else
        work(1) = 0.5 * sqrt(PI/T) * erf(sqrt(T))
        do i = 2, n+1
            work(i) = ( (2.0*(i-1.0) - 1.0) * work(i-1) - exp(-T) ) / (2.0*T)
        end do
        boys = work(n+1)
    end if
end function

recursive function e(n, m, T, AB, exp1, exp2) result(res)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: T
    real*8, intent(in)  :: AB
    real*8, intent(in)  :: exp1
    real*8, intent(in)  :: exp2

    real*8 :: p, q, res

    p = exp1 + exp2
    q = exp1 * exp2 / p

    if ( (T .lt. 0) .or. (T .gt. (n+m))) then
        res = 0.0
    else if ((T .eq. 0 .and. n .eq. 0) .and. m .eq. 0) then
        res = exp(-q * AB**2)
    else if (n .eq. 0) then
        res = 1.0/(2.0*p) * e(n, m-1, T-1, AB, exp1, exp2) &
            + exp1*AB/p * e(n, m-1, T, AB, exp1, exp2) &
            + (T+1.0) * e(n, m-1, T+1, AB, exp1, exp2)
    else 
        res = 1.0/(2.0*p) * e(n-1, m, T-1, AB, exp1, exp2) &
            - exp2*AB/p * e(n-1, m, T, AB, exp1, exp2) &
            + (T+1.0) * e(n-1, m, T+1, AB, exp1, exp2)
    end if
end function e

recursive function r(t, u, v, n, p, X_PC, Y_PC, Z_PC, R_pc2) result(res)
    implicit none
    integer, intent(in) :: t, u, v, n
    real*8, intent(in) :: p, X_PC, Y_PC, Z_PC, R_pc2

    real*8 :: boys
    real*8 :: res
    real*8 :: tmp
    tmp = 0.0
    if ((t .eq. 0 .and. u .eq. 0) .and. v .eq. 0) then
        res = (-2.0*p)**n * boys(p*R_pc2, n)
    else if (t .eq. 0 .and. u .eq. 0) then
        if (v .gt. 1) then
            tmp = (v-1)*r(t, u, v-2, n+1, p, X_PC, Y_PC, Z_PC, R_pc2)
        end if
        res = tmp + Z_PC*r(t, u, v-1, n+1, p, X_PC, Y_PC, Z_PC, R_pc2)
    else if (t .eq. 0) then
        if (u .gt. 1) then
            tmp = (u-1)*r(t, u-2, v, n+1, p, X_PC, Y_PC, Z_PC, R_pc2)
        end if
        res = tmp + Y_PC*r(t, u-1, v, n+1, p, X_PC, Y_PC, Z_PC, R_pc2)
    else 
        if (t .gt. 1) then
            tmp = (t-1)*r(t-2, u, v, n+1, p, X_PC, Y_PC, Z_PC, R_pc2)
        end if
        res = tmp + X_PC*r(t-1, u, v, n+1, p, X_PC, Y_PC, Z_PC, R_pc2)
    end if
end function r

function overlap(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2) result(S)
    implicit none
    integer, intent(in) :: ng1
    real*8, dimension(ng1), intent(in) :: coeffs1
    real*8, dimension(ng1), intent(in) :: exps1
    integer, dimension(3), intent(in) :: l1
    real*8, dimension(3), intent(in) :: xyz1
    integer, intent(in) :: ng2
    real*8, dimension(ng2), intent(in) :: coeffs2
    real*8, dimension(ng2), intent(in) :: exps2
    integer, dimension(3), intent(in) :: l2
    real*8, dimension(3), intent(in) :: xyz2

    real*8 S, overlap_element
    integer i, j

    S = 0.0
    do i=1, ng1
        do j=1, ng2
            S  = S + coeffs1(i) * coeffs2(j) &
               * overlap_element(exps1(i), l1, xyz1, exps2(j), l2, xyz2)
        end do
    end do

end function overlap

function overlap_element(exp1, l1, xyz1, exp2, l2, xyz2) result(S_inner)
    implicit none
    real*8, intent(in) :: exp1
    integer, dimension(3), intent(in) :: l1
    real*8, dimension(3),intent(in) :: xyz1
    real*8, intent(in) :: exp2
    integer, dimension(3), intent(in) :: l2
    real*8, dimension(3), intent(in) :: xyz2

    real*8 :: S_inner, e
    real*8 :: PI = 3.14159265358979323846264338327950288419

    S_inner = e(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2) &
            * e(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2) &
            * e(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2) &
            * (PI / (exp1 + exp2))**1.50

end function overlap_element

function kinetic(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2) result(T)
    implicit none
    integer, intent(in) :: ng1
    real*8, dimension(ng1), intent(in) :: coeffs1
    real*8, dimension(ng1), intent(in) :: exps1
    integer, dimension(3), intent(in) :: l1
    real*8, dimension(3), intent(in) :: xyz1
    integer, intent(in) :: ng2
    real*8, dimension(ng2), intent(in) :: coeffs2
    real*8, dimension(ng2), intent(in) :: exps2
    integer, dimension(3), intent(in) :: l2
    real*8, dimension(3), intent(in) :: xyz2

    real*8 :: T, tmp1, tmp2, tmp3, overlap_element
    integer :: i, j

    T = 0.0
    do i=1, ng1
        do j=1, ng2
            tmp1 = exps2(j) * (2 * (l2(1) + l2(2) + l2(3)) + 3) * overlap_element(exps1(i), l1, xyz1, exps2(j), l2, xyz2)
            tmp2 = -2 * exps2(j) ** 2 * (overlap_element(exps1(i), l1, xyz1, exps2(j), l2 + (/2, 0, 0/), xyz2) &
                                         + overlap_element(exps1(i), l1, xyz1, exps2(j), l2 + (/0, 2, 0/), xyz2) &
                                         + overlap_element(exps1(i), l1, xyz1, exps2(j), l2 + (/0, 0, 2/), xyz2))
            tmp3 = -0.5 * (l2(1) * (l2(1) - 1) * overlap_element(exps1(i), l1, xyz1, exps2(j), l2 + (/-2, 0, 0/), xyz2) &
                           + l2(2) * (l2(2) - 1) * overlap_element(exps1(i), l1, xyz1, exps2(j), l2 + (/0, -2, 0/), xyz2) &
                           + l2(3) * (l2(3) - 1) * overlap_element(exps1(i), l1, xyz1, exps2(j), l2 + (/0, 0, -2/), xyz2))
            T = T + coeffs1(i)*coeffs2(j)*(tmp1 + tmp2 + tmp3)
        end do 
    end do

end function kinetic

function potential_1e(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, xyza) result(V)
    implicit none
    integer, intent(in) :: ng1
    real*8, dimension(ng1), intent(in) :: coeffs1
    real*8, dimension(ng1), intent(in) :: exps1
    integer, dimension(3), intent(in) :: l1
    real*8, dimension(3), intent(in) :: xyz1
    integer, intent(in) :: ng2
    real*8, dimension(ng2), intent(in) :: coeffs2
    real*8, dimension(ng2), intent(in) :: exps2
    integer, dimension(3), intent(in) :: l2
    real*8, dimension(3), intent(in) :: xyz2
    real*8, dimension(3), intent(in) :: xyza

    integer :: i, j
    real*8 :: V, nuc_attraction

    V = 0.0

    do i=1, ng1
        do j=1, ng2
            V = V + coeffs1(i)*coeffs2(j)*nuc_attraction(exps1(i), l1, xyz1, exps2(j), l2,xyz2, xyza)
        end do
    end do

end function potential_1e

function nuc_attraction(exp1, l1, xyz1, exp2, l2, xyz2, xyza) result(nuc_att)
    implicit none
    real*8, intent(in) :: exp1, exp2
    integer, intent(in) :: l1(3), l2(3)
    real*8, intent(in) :: xyz1(3), xyz2(3), xyza(3)

    real*8 :: p
    real*8, dimension(3) :: xyzp
    real*8 :: R_pc2
    real*8, parameter :: PI = 3.14159265358979323846264338327950288419

    real*8 :: e, r, nuc_att
    integer :: t, u , v

    p = exp1 + exp2
    xyzp = (xyz1 * exp1 + xyz2 * exp2) / p
    R_pc2 = (xyzp(1) - xyza(1))**2 + (xyzp(2) - xyza(2))**2 + (xyzp(3)-xyza(3))**2
    nuc_att = 0.0


    do t=0, (l1(1)+l2(1))
        do u=0, (l1(2)+l2(2))
            do v=0, (l1(3)+l2(3))
                nuc_att = nuc_att + e(l1(1), l2(1), t, (xyz1(1)-xyz2(1)), exp1, exp2) &
                            * e(l1(2), l2(2), u, (xyz1(2)-xyz2(2)), exp1, exp2) &
                            * e(l1(3), l2(3), v, (xyz1(3)-xyz2(3)), exp1, exp2) &
                            * r(t, u, v, 0, p, xyzp(1)-xyza(1), xyzp(2)-xyza(2), xyzp(3)-xyza(3), R_pc2)
            end do
        end do
    end do
    nuc_att = nuc_att * 2.0 * PI / p

end function nuc_attraction

function potential_2e(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, &
                      ng3, coeffs3, exps3, l3, xyz3, ng4, coeffs4, exps4, l4, xyz4) result(eri)
    implicit none
    integer, intent(in) :: ng1, ng2, ng3, ng4
    real*8, intent(in) :: coeffs1(ng1), coeffs2(ng2), coeffs3(ng3), coeffs4(ng4)
    real*8, intent(in) :: exps1(ng1), exps2(ng2), exps3(ng3), exps4(ng4)
    integer, dimension(3), intent(in) :: l1, l2, l3, l4
    real*8, dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4

    real*8 :: eri, electron_repulsion
    integer :: i, j, k, l

    eri = 0.0

    do i=1, ng1
        do j=1, ng2
            do k=1, ng3
                do l=1, ng4
                    eri = eri + coeffs1(i) * coeffs2(j) * coeffs3(k) * coeffs4(l) &
                                * electron_repulsion(exps1(i), l1, xyz1, exps2(j), l2, xyz2, &
                                                     exps3(k), l3, xyz3, exps4(l), l4, xyz4)
                end do
            end do
        end do
    end do

end function potential_2e

function electron_repulsion(exp1, l1, xyz1, exp2, l2, xyz2, exp3, l3, xyz3, exp4, l4, xyz4) result(el_rep)
    implicit none
    real*8, intent(in) :: exp1, exp2, exp3, exp4
    integer, dimension(3), intent(in) :: l1, l2, l3, l4
    real*8, dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4

    real*8 :: p, q, el_rep, R_pq2, e, r
    real*8, dimension(3) :: xyzp, xyzq
    integer :: t, u, v, tau, nu, phi
    real*8, parameter :: PI = 3.14159265358979323846264338327950288419

    p = exp1 +exp2
    q = exp3 + exp4
    xyzp = (exp1 * xyz1 + exp2 * xyz2) / p
    xyzq = (exp3 * xyz3 + exp4 * xyz4) / q
    R_pq2 = (xyzp(1) - xyzq(1))**2 + (xyzp(2) - xyzq(2))**2 + (xyzp(3) - xyzq(3))**2

    el_rep = 0.0

    do t=0, (l1(1)+l2(1))
        do u=0, (l1(2)+l2(2))
            do v=0, (l1(3)+l2(3))
                do tau=0, (l3(1)+l4(1))
                    do nu=0, (l3(2)+l4(2))
                        do phi=0, (l3(3)+l4(3))
                            el_rep = el_rep + e(l1(1), l2(1), t, (xyz1(1) - xyz2(1)), exp1, exp2) &
                                              * e(l1(2), l2(2), u, (xyz1(2) - xyz2(2)), exp1, exp2) &
                                              * e(l1(3), l2(3), v, (xyz1(3) - xyz2(3)), exp1, exp2) &
                                              * e(l3(1), l4(1), tau, (xyz3(1) - xyz4(1)), exp3, exp4) &
                                              * e(l3(2), l4(2), nu, (xyz3(2) - xyz4(2)), exp3, exp4) &
                                              * e(l3(3), l4(3), phi, (xyz3(3) - xyz4(3)), exp3, exp4) &
                                              * r(t+tau, u+nu, v+phi, 0, p*q/(p+q), &
                                                  xyzp(1)-xyzq(1), xyzp(2)-xyzq(2), xyzp(3)-xyzq(3), R_pq2) &
                                              * (-1)**(tau+nu+phi)
                        end do
                    end do
                end do
            end do
        end do
    end do

    el_rep = el_rep * 2 * PI ** (2.50) / (p * q * sqrt(p + q))
end function electron_repulsion