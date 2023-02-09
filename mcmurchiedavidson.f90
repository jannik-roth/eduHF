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