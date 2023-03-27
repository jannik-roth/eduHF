module globals
    implicit none
    real*8 :: PI = 3.14159265358979323846264338327950288419
end module globals


function boys(T, n)
    use globals
    implicit none
    real*8 :: boys
    real*8, intent(in) :: T
    integer, intent(in) :: n

    real*8 :: work(n+1)
    integer :: i

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

function e_der(n, m, T, AB, exp1, exp2, center) result(res)
    implicit none
    integer, intent(in) :: n, m, T, center
    real*8, intent(in) :: AB, exp1, exp2

    real*8 :: res, e
    res = 0.0

    if (center .eq. 0) then
        res = 2.d0 * exp1 * e(n+1, m, T, AB, exp1, exp2)
        if (n .gt. 0) then
            res = res - n * e (n-1, m, T, AB, exp1, exp2)
        end if
    else if (center .eq. 1) then
        res = 2.d0 * exp2 * e(n, m+1, T, AB, exp1, exp2)
        if (m .gt. 0) then
            res = res - m * e(n, m-1, T, AB, exp1, exp2)
        end if
    end if

end function e_der

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

function overlap_der(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, center, dim) result(S_der)
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
    integer, intent(in) :: center, dim

    real*8 :: S_der, overlap_der_element
    integer :: i, j

    S_der = 0.0
    do i=1, ng1
        do j=1, ng2
            S_der  = S_der + coeffs1(i) * coeffs2(j) &
               * overlap_der_element(exps1(i), l1, xyz1, exps2(j), l2, xyz2, center, dim)
        end do
    end do

end function overlap_der

function overlap_element(exp1, l1, xyz1, exp2, l2, xyz2) result(S_inner)
    use globals
    implicit none
    real*8, intent(in) :: exp1
    integer, dimension(3), intent(in) :: l1
    real*8, dimension(3),intent(in) :: xyz1
    real*8, intent(in) :: exp2
    integer, dimension(3), intent(in) :: l2
    real*8, dimension(3), intent(in) :: xyz2

    real*8 :: S_inner, e

    S_inner = e(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2) &
            * e(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2) &
            * e(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2) &
            * (PI / (exp1 + exp2))**1.50

end function overlap_element

function overlap_der_element(exp1, l1, xyz1, exp2, l2, xyz2, center, dim) result(res)
    use globals
    implicit none
    real*8, intent(in) :: exp1, exp2
    integer, dimension(3), intent(in) :: l1, l2
    real*8, dimension(3), intent(in) :: xyz1, xyz2
    integer, intent(in) :: center, dim

    real*8 :: res, tmp1, tmp2, tmp3, e, e_der
    tmp1 = 0.0
    tmp2 = 0.0
    tmp3 = 0.0
    
    if (dim .eq. 0) then
        tmp1 = e_der(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2, center)
        tmp2 = e(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2)
        tmp3 = e(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2)
    else if (dim .eq. 1) then
        tmp1 = e(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2)
        tmp2 = e_der(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2, center)
        tmp3 = e(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2)
    else if (dim .eq. 2) then
        tmp1 = e(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2)
        tmp2 = e(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2)
        tmp3 = e_der(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2, center)
    end if
    res = tmp1 * tmp2 * tmp3 * (PI / (exp1 + exp2)) ** 1.50
end function

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

function kinetic_der_element(exp1, l1, xyz1, exp2, l2, xyz2, center, dim) result(res)
    use globals
    implicit none
    real*8, intent(in) :: exp1, exp2
    integer, dimension(3), intent(in) :: l1, l2
    real*8, dimension(3), intent(in) :: xyz1, xyz2
    integer, intent(in) :: center, dim

    real*8 :: res, kin1, kin2, kin3, e_der, e, tmp

    ! derivative wrt to the x coordinate
    if (dim .eq. 0) then
        kin1 = exp2 * (2.0*l2(1)+1.0) * e_der(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2, center) &
               - 2.0 * exp2**2 * e_der(l1(1), l2(1)+2, 0, xyz1(1)-xyz2(1), exp1, exp2, center) &
               - 0.5 * l2(1) * (l2(1)-1.0) * e_der(l1(1), l2(1)-2, 0, xyz1(1)-xyz2(1), exp1, exp2, center)
    else
        kin1 = exp2 * (2.0*l2(1)+1.0) * e(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2) &
               - 2.0 * exp2**2 * e(l1(1), l2(1)+2, 0, xyz1(1)-xyz2(1), exp1, exp2) &
               - 0.5 * l2(1) * (l2(1)-1.0) * e(l1(1), l2(1)-2, 0, xyz1(1)-xyz2(1), exp1, exp2)
    end if

    ! derivative wrt to the y coordinate
    if (dim .eq. 1) then
        kin2 = exp2 * (2.0*l2(2)+1.0) * e_der(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2, center) &
               - 2.0 * exp2**2 * e_der(l1(2), l2(2)+2, 0, xyz1(2)-xyz2(2), exp1, exp2, center) &
               - 0.5 * l2(2) * (l2(2)-1.0) * e_der(l1(2), l2(2)-2, 0, xyz1(2)-xyz2(2), exp1, exp2, center)
    else
        kin2 = exp2 * (2.0*l2(2)+1.0) * e(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2) &
               - 2.0 * exp2**2 * e(l1(2), l2(2)+2, 0, xyz1(2)-xyz2(2), exp1, exp2) &
               - 0.5 * l2(2) * (l2(2)-1.0) * e(l1(2), l2(2)-2, 0, xyz1(2)-xyz2(2), exp1, exp2)
    end if

    ! derivative wrt to the z coordinate
    if (dim .eq. 2) then
        kin3 = exp2 * (2.0*l2(3)+1.0) * e_der(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2, center) &
               - 2.0 * exp2**2 * e_der(l1(3), l2(3)+2, 0, xyz1(3)-xyz2(3), exp1, exp2, center) &
               - 0.5 * l2(3) * (l2(3)-1.0) * e_der(l1(3), l2(3)-2, 0, xyz1(3)-xyz2(3), exp1, exp2, center)
    else
        kin3 = exp2 * (2.0*l2(3)+1.0) * e(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2) &
               - 2.0 * exp2**2 * e(l1(3), l2(3)+2, 0, xyz1(3)-xyz2(3), exp1, exp2) &
               - 0.5 * l2(3) * (l2(3)-1.0) * e(l1(3), l2(3)-2, 0, xyz1(3)-xyz2(3), exp1, exp2)
    end if

    if (dim .eq. 0) then
        tmp = e_der(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2, center)
        kin2 = kin2 * tmp
        kin3 = kin3 * tmp
    else
        tmp = e(l1(1), l2(1), 0, xyz1(1)-xyz2(1), exp1, exp2)
        kin2 = kin2 * tmp
        kin3 = kin3 * tmp
    end if

    if (dim .eq. 1) then
        tmp = e_der(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2, center)
        kin1 = kin1 * tmp
        kin3 = kin3 * tmp
    else
        tmp = e(l1(2), l2(2), 0, xyz1(2)-xyz2(2), exp1, exp2)
        kin1 = kin1 * tmp
        kin3 = kin3 * tmp
    end if

    if (dim .eq. 2) then
        tmp = e_der(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2, center)
        kin1 = kin1 * tmp
        kin2 = kin2 * tmp
    else
        tmp = e(l1(3), l2(3), 0, xyz1(3)-xyz2(3), exp1, exp2)
        kin1 = kin1 * tmp
        kin2 = kin2 * tmp
    end if
  
    res = (kin1 + kin2 + kin3) * (PI / (exp1 + exp2))**1.50

end function kinetic_der_element

function kinetic_der(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, center, dim) result(T_der)
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
    integer, intent(in) :: center, dim

    real*8 :: T_der, kinetic_der_element
    integer :: i, j

    T_der = 0.0
    do i=1, ng1
        do j=1, ng2
            T_der  = T_der + coeffs1(i) * coeffs2(j) &
                     * kinetic_der_element(exps1(i), l1, xyz1, exps2(j), l2, xyz2, center, dim)
        end do
    end do

end function kinetic_der

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
    use globals
    implicit none
    real*8, intent(in) :: exp1, exp2
    integer, intent(in) :: l1(3), l2(3)
    real*8, intent(in) :: xyz1(3), xyz2(3), xyza(3)

    real*8 :: p
    real*8, dimension(3) :: xyzp
    real*8 :: R_pc2

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

function nuc_attraction_der(exp1, l1, xyz1, exp2, l2, xyz2, xyza, dim) result(res)
    ! Always differentiate the first center !
    use globals
    implicit none
    real*8, intent(in) :: exp1, exp2
    integer, intent(in) :: l1(3), l2(3), dim
    real*8, intent(in) :: xyz1(3), xyz2(3), xyza(3)

    real*8 :: p
    real*8, dimension(3) :: xyzp
    real*8 :: R_pc2

    real*8 :: e, e_der, r, res
    integer :: t, u , v

    p = exp1 + exp2
    xyzp = (xyz1 * exp1 + xyz2 * exp2) / p
    R_pc2 = (xyzp(1) - xyza(1))**2 + (xyzp(2) - xyza(2))**2 + (xyzp(3)-xyza(3))**2
    
    res = 0.0
    ! derivative wrt to x
    if (dim .eq. 0) then
        do t=0, (l1(1)+l2(1)+1+1)
            do u=0, (l1(2)+l2(2)+1)
                do v=0, (l1(3)+l2(3)+1)
                    res = res + e_der(l1(1), l2(1), t, xyz1(1)-xyz2(1), exp1, exp2, 0) &
                                * e(l1(2), l2(2), u, xyz1(2)-xyz2(2), exp1, exp2) &
                                * e(l1(3), l2(3), v, xyz1(3)-xyz2(3), exp1, exp2) &
                                * r(t, u, v, 0, p, xyzp(1)-xyza(1), xyzp(2)-xyza(2), xyzp(3)-xyza(3), R_pc2)
                end do
            end do
        end do
    ! derivative wrt to y
    else if (dim .eq. 1) then
        do t=0, (l1(1)+l2(1)+1)
            do u=0, (l1(2)+l2(2)+1+1)
                do v=0, (l1(3)+l2(3)+1)
                    res = res + e(l1(1), l2(1), t, xyz1(1)-xyz2(1), exp1, exp2) &
                                * e_der(l1(2), l2(2), u, xyz1(2)-xyz2(2), exp1, exp2, 0) &
                                * e(l1(3), l2(3), v, xyz1(3)-xyz2(3), exp1, exp2) &
                                * r(t, u, v, 0, p, xyzp(1)-xyza(1), xyzp(2)-xyza(2), xyzp(3)-xyza(3), R_pc2)
                end do
            end do
        end do
    ! derivative wrt to z
    else if (dim .eq. 2) then
        do t=0, (l1(1)+l2(1)+1)
            do u=0, (l1(2)+l2(2)+1)
                do v=0, (l1(3)+l2(3)+1+1)
                    res = res + e(l1(1), l2(1), t, xyz1(1)-xyz2(1), exp1, exp2) &
                                * e(l1(2), l2(2), u, xyz1(2)-xyz2(2), exp1, exp2) &
                                * e_der(l1(3), l2(3), v, xyz1(3)-xyz2(3), exp1, exp2, 0) &
                                * r(t, u, v, 0, p, xyzp(1)-xyza(1), xyzp(2)-xyza(2), xyzp(3)-xyza(3), R_pc2)
                end do
            end do
        end do
    end if

    res = res * 2.0 * PI / p
end function nuc_attraction_der

function potential_1e_der_ordered(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, xyza, dim) result(V_der)
    ! always differentiate the first function
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
    integer, intent(in) :: dim

    integer :: i, j
    real*8 :: V_der, nuc_attraction_der

    V_der = 0.0

    do i=1, ng1
        do j=1, ng2
            V_der = V_der + coeffs1(i)*coeffs2(j)*nuc_attraction_der(exps1(i), l1, xyz1, exps2(j), l2,xyz2, xyza, dim)
        end do
    end do

end function potential_1e_der_ordered

function potential_1e_der(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, xyza, center, dim) result(V_der)
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
    integer, intent(in) :: center, dim

    real*8 :: V_der, potential_1e_der_ordered

    V_der = 0.0
    if (center .eq. 0) then
        V_der = potential_1e_der_ordered(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, xyza, dim)
    else if (center .eq. 1) then
        V_der = potential_1e_der_ordered(ng2, coeffs2, exps2, l2, xyz2, ng1, coeffs1, exps1, l1, xyz1, xyza, dim)
    end if

end function potential_1e_der

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
    use globals
    implicit none
    real*8, intent(in) :: exp1, exp2, exp3, exp4
    integer, dimension(3), intent(in) :: l1, l2, l3, l4
    real*8, dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4

    real*8 :: p, q, el_rep, R_pq2, e, r
    real*8, dimension(3) :: xyzp, xyzq
    integer :: t, u, v, tau, nu, phi

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

function electron_repulsion_der(exp1, l1, xyz1, exp2, l2, xyz2, exp3, l3, xyz3, exp4, l4, xyz4, dim) result(res)
    use globals
    ! always differentiate the first bf
    implicit none
    real*8, intent(in) :: exp1, exp2, exp3, exp4
    integer, dimension(3), intent(in) :: l1, l2, l3, l4
    real*8, dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4
    integer, intent(in) :: dim

    real*8 :: p, q, R_pq2, e, e_der, r, res
    real*8, dimension(3) :: xyzp, xyzq
    integer :: t, u, v, tau, nu, phi

    p = exp1 +exp2
    q = exp3 + exp4
    xyzp = (exp1 * xyz1 + exp2 * xyz2) / p
    xyzq = (exp3 * xyz3 + exp4 * xyz4) / q
    R_pq2 = (xyzp(1) - xyzq(1))**2 + (xyzp(2) - xyzq(2))**2 + (xyzp(3) - xyzq(3))**2

    res = 0.0
    ! differentiate wrt to x
    if (dim .eq. 0) then
        do t=0, (l1(1)+l2(1)+1)
            do u=0, (l1(2)+l2(2))
                do v=0, (l1(3)+l2(3))
                    do tau=0, (l3(1)+l4(1))
                        do nu=0, (l3(2)+l4(2))
                            do phi=0, (l3(3)+l4(3))
                                res = res + e_der(l1(1), l2(1), t, xyz1(1)-xyz2(1), exp1, exp2, 0) &
                                            * e(l1(2), l2(2), u, xyz1(2)-xyz2(2), exp1, exp2) &
                                            * e(l1(3), l2(3), v, xyz1(3)-xyz2(3), exp1, exp2) &
                                            * e(l3(1), l4(1), tau, xyz3(1)-xyz4(1), exp3, exp4) &
                                            * e(l3(2), l4(2), nu, xyz3(2)-xyz4(2), exp3, exp4) &
                                            * e(l3(3), l4(3), phi, xyz3(3)-xyz4(3), exp3, exp4) &
                                            * r(t+tau, u+nu, v+phi, 0, p*q/(p+q), &
                                                xyzp(1)-xyzq(1), xyzp(2)-xyzq(2), xyzp(3)-xyzq(3), R_pq2) &
                                            * (-1.0)**(tau + nu + phi)
                                
                            end do
                        end do
                    end do
                end do
            end do
        end do
    ! differentiate wrt to y
    else if (dim .eq. 1) then
        do t=0, (l1(1)+l2(1))
            do u=0, (l1(2)+l2(2)+1)
                do v=0, (l1(3)+l2(3))
                    do tau=0, (l3(1)+l4(1))
                        do nu=0, (l3(2)+l4(2))
                            do phi=0, (l3(3)+l4(3))
                                res = res + e(l1(1), l2(1), t, xyz1(1)-xyz2(1), exp1, exp2) &
                                            * e_der(l1(2), l2(2), u, xyz1(2)-xyz2(2), exp1, exp2, 0) &
                                            * e(l1(3), l2(3), v, xyz1(3)-xyz2(3), exp1, exp2) &
                                            * e(l3(1), l4(1), tau, xyz3(1)-xyz4(1), exp3, exp4) &
                                            * e(l3(2), l4(2), nu, xyz3(2)-xyz4(2), exp3, exp4) &
                                            * e(l3(3), l4(3), phi, xyz3(3)-xyz4(3), exp3, exp4) &
                                            * r(t+tau, u+nu, v+phi, 0, p*q/(p+q), &
                                                xyzp(1)-xyzq(1), xyzp(2)-xyzq(2), xyzp(3)-xyzq(3), R_pq2) &
                                            * (-1.0)**(tau + nu + phi)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    ! differentiate wrt to z
    else if (dim .eq. 2) then
        do t=0, (l1(1)+l2(1))
            do u=0, (l1(2)+l2(2))
                do v=0, (l1(3)+l2(3)+1)
                    do tau=0, (l3(1)+l4(1))
                        do nu=0, (l3(2)+l4(2))
                            do phi=0, (l3(3)+l4(3))
                                res = res + e(l1(1), l2(1), t, xyz1(1)-xyz2(1), exp1, exp2) &
                                            * e(l1(2), l2(2), u, xyz1(2)-xyz2(2), exp1, exp2) &
                                            * e_der(l1(3), l2(3), v, xyz1(3)-xyz2(3), exp1, exp2, 0) &
                                            * e(l3(1), l4(1), tau, xyz3(1)-xyz4(1), exp3, exp4) &
                                            * e(l3(2), l4(2), nu, xyz3(2)-xyz4(2), exp3, exp4) &
                                            * e(l3(3), l4(3), phi, xyz3(3)-xyz4(3), exp3, exp4) &
                                            * r(t+tau, u+nu, v+phi, 0, p*q/(p+q), &
                                                xyzp(1)-xyzq(1), xyzp(2)-xyzq(2), xyzp(3)-xyzq(3), R_pq2) &
                                            * (-1.0)**(tau + nu + phi)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if

    res = res * 2.0 * PI**2.50 / (p * q * sqrt(p + q))

end function electron_repulsion_der

function potential_2e_der_ordered(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, &
                                  ng3, coeffs3, exps3, l3, xyz3, ng4, coeffs4, exps4, l4, xyz4, &
                                  dim) result(res)
    implicit none
    integer, intent(in) :: ng1, ng2, ng3, ng4
    real*8, intent(in) :: coeffs1(ng1), coeffs2(ng2), coeffs3(ng3), coeffs4(ng4)
    real*8, intent(in) :: exps1(ng1), exps2(ng2), exps3(ng3), exps4(ng4)
    integer, dimension(3), intent(in) :: l1, l2, l3, l4
    real*8, dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4
    integer, intent(in) :: dim

    real*8 :: res, electron_repulsion_der
    integer :: i, j, k, l

    res = 0.0

    do i=1, ng1
        do j=1, ng2
            do k=1, ng3
                do l=1, ng4
                    res = res + coeffs1(i) * coeffs2(j) * coeffs3(k) * coeffs4(l) &
                                * electron_repulsion_der(exps1(i), l1, xyz1, exps2(j), l2, xyz2, &
                                                         exps3(k), l3, xyz3, exps4(l), l4, xyz4, &
                                                         dim)
                end do
            end do
        end do
    end do

end function potential_2e_der_ordered

function potential_2e_der(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, &
                          ng3, coeffs3, exps3, l3, xyz3, ng4, coeffs4, exps4, l4, xyz4, &
                          center, dim) result(res)
    implicit none
    integer, intent(in) :: ng1, ng2, ng3, ng4
    real*8, intent(in) :: coeffs1(ng1), coeffs2(ng2), coeffs3(ng3), coeffs4(ng4)
    real*8, intent(in) :: exps1(ng1), exps2(ng2), exps3(ng3), exps4(ng4)
    integer, dimension(3), intent(in) :: l1, l2, l3, l4
    real*8, dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4
    integer, intent(in) :: center, dim

    real*8 :: res, potential_2e_der_ordered

    res = 0.0

    if (center .eq. 0) then
        res = potential_2e_der_ordered(ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, &
                                       ng3, coeffs3, exps3, l3, xyz3, ng4, coeffs4, exps4, l4, xyz4, &
                                       dim)
    else if (center .eq. 1) then 
        res = potential_2e_der_ordered(ng2, coeffs2, exps2, l2, xyz2, ng1, coeffs1, exps1, l1, xyz1, &
                                       ng4, coeffs4, exps4, l4, xyz4, ng3, coeffs3, exps3, l3, xyz3, &
                                       dim)
    else if (center .eq. 2) then 
        res = potential_2e_der_ordered(ng3, coeffs3, exps3, l3, xyz3, ng4, coeffs4, exps4, l4, xyz4, &
                                       ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, &
                                       dim)
    else if (center .eq. 3) then
        res = potential_2e_der_ordered(ng4, coeffs4, exps4, l4, xyz4, ng3, coeffs3, exps3, l3, xyz3, &
                                       ng1, coeffs1, exps1, l1, xyz1, ng2, coeffs2, exps2, l2, xyz2, &
                                       dim)
    end if
end function potential_2e_der