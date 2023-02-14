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

recursive function e(n, m, T, AB, alpha, beta) result(res)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: T
    real*8, intent(in)  :: AB
    real*8, intent(in)  :: alpha
    real*8, intent(in)  :: beta

    real*8 :: p, q, res

    p = alpha + beta
    q = alpha * beta / p

    if ( (T .lt. 0) .or. (T .gt. (n+m))) then
        res = 0.0
    else if ((T .eq. 0 .and. n .eq. 0) .and. m .eq. 0) then
        res = exp(-q * AB**2)
    else if (n .eq. 0) then
        res = 1.0/(2.0*p) * e(n, m-1, T-1, AB, alpha, beta) &
            + alpha*AB/p * e(n, m-1, T, AB, alpha, beta) &
            + (T+1.0) * e(n, m-1, T+1, AB, alpha, beta)
    else 
        res = 1.0/(2.0*p) * e(n-1, m, T-1, AB, alpha, beta) &
            - beta*AB/p * e(n-1, m, T, AB, alpha, beta) &
            + (T+1.0) * e(n-1, m, T+1, AB, alpha, beta)
    end if
end function e

recursive function r(t, u, v, n, p, X_PC, Y_PC, Z_PC, R_pc2) result(res)
    implicit none
    integer, intent(in) :: t, u, v, n
    real*8, intent(in) :: p, X_PC, Y_PC, Z_PC, R_pc2

    real*8 :: boys
    real*8 :: res
    real*8 :: tmp = 0.0
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

function overlap(ng1, coeffs1, exps1, l1, center1, ng2, coeffs2, exps2, l2, center2) result(S)
    implicit none
    integer, intent(in) :: ng1
    real*8, dimension(ng1), intent(in) :: coeffs1
    real*8, dimension(ng1), intent(in) :: exps1
    integer, dimension(3), intent(in) :: l1
    real*8, dimension(3), intent(in) :: center1
    integer, intent(in) :: ng2
    real*8, dimension(ng2), intent(in) :: coeffs2
    real*8, dimension(ng2), intent(in) :: exps2
    integer, dimension(3), intent(in) :: l2
    real*8, dimension(3), intent(in) :: center2

    real*8 S, overlap_element
    integer i, j

    S = 0.0
    do i=1, ng1
        do j=1, ng2
            S  = S + coeffs1(i) * coeffs2(j) &
               * overlap_element(exps1(i), l1, center1, exps2(j), l2, center2)
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

    real*8 S_inner, e
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

function nuc_attraction(alpha, l1, xyz1, beta, l2, xyz2, xyza) result(nuc_att)
    implicit none
    real*8, intent(in) :: alpha, beta
    integer, intent(in) :: l1(3), l2(3)
    real*8, intent(in) :: xyz1(3), xyz2(3), xyza(3)

    real*8 :: p
    real*8, dimension(3) :: xyzp
    real*8 :: R_pc2
    real*8, parameter :: PI = 3.14159265358979323846264338327950288419

    real*8 :: e, r, nuc_att
    integer :: t, u , v

    p = alpha + beta
    xyzp = (xyz1 * alpha + xyz2 * beta) / p
    R_pc2 = (xyzp(1) - xyza(1))**2 + (xyzp(2) - xyza(2))**2 + (xyzp(3)-xyza(3))**2
    nuc_att = 0.0


    do t=0, (l1(1)+l2(1))
        do u=0, (l1(2)+l2(2))
            do v=0, (l1(3)+l2(3))
                nuc_att = nuc_att + e(l1(1), l2(1), t, (xyz1(1)-xyz2(1)), alpha, beta) &
                            * e(l1(2), l2(2), u, (xyz1(2)-xyz2(2)), alpha, beta) &
                            * e(l1(3), l2(3), v, (xyz1(3)-xyz2(3)), alpha, beta) &
                            * r(t, u, v, 0, p, xyzp(1)-xyza(1), xyzp(2)-xyza(2), xyzp(3)-xyza(3), R_pc2)
            end do
        end do
    end do
    nuc_att = nuc_att * 2.0 * PI / p

end function nuc_attraction

! function testin(nob, vec1) result(tmp)
!     implicit none
!     integer*4, intent(in) :: nob
!     real*8, dimension(nob), intent(in) :: vec1


!     real*8 :: tmp
!     real*8, dimension(nob) :: testing
!     testing = vec1 + (/0, 0, 10, 0, 0/)
!     print *, testing
    
! end function