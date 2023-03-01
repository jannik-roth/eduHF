subroutine slater_exp(exps, coeffs, zeta, ng, nl_quant)
    implicit none
    real*8, intent(inout)    :: exps(ng)
    real*8, intent(inout)    :: coeffs(ng)
    real, intent(in)     :: zeta
    integer, intent(in)      :: ng
    character(len = 2), intent(in) :: nl_quant
    integer :: l, nl_to_l
    real, parameter :: double_fac(8) = (/1.0,1.0,3.0,15.0,105.0,945.0,10395.0,135135.0/) ! OEIS A001147
    real, parameter :: two_over_pi = 0.636619772367582

    l = nl_to_l(nl_quant) ! get l quantum number s = 0, p = 1, ...
    
    select case(ng)
    case(1)
        call sto1g(exps, coeffs, zeta, nl_quant)
    case(2)
        call sto2g(exps, coeffs, zeta, nl_quant)
    case(3)
        call sto3g(exps, coeffs, zeta, nl_quant)    
    end select

    ! need to normalize
    coeffs = coeffs * (two_over_pi * exps)**(0.75) * sqrt(4*exps)**(l / sqrt(double_fac(l+1)))
end subroutine

subroutine sto1g(exps, coeffs, zeta, nl_quant)
    ! Taken from: The Journal of Chemical Physics 52, 431 (1970)
    ! http://dx.doi.org/10.1063/1.1672702
    implicit none
    real, intent(in)         :: zeta
    character(2), intent(in) :: nl_quant
    real*8, dimension(1), intent(inout)    :: exps
    real*8, dimension(1), intent(inout)    :: coeffs

    select case(nl_quant)
    case("1s"); exps(1) = 2.709498091e-1
    case("2s"); exps(1) = 1.012151084e-1
    case("3s"); exps(1) = 5.296881757e-2
    case("4s"); exps(1) = 3.264600274e-2
    case("5s"); exps(1) = 2.216912938e-2
    case("2p"); exps(1) = 1.759666885e-1
    case("3p"); exps(1) = 9.113614253e-2
    case("4p"); exps(1) = 5.578350235e-2
    case("5p"); exps(1) = 3.769845216e-2
    case("3d"); exps(1) = 1.302270363e-1
    case("4d"); exps(1) = 7.941656339e-2
    case("5d"); exps(1) = 5.352200793e-2
    case("4f"); exps(1) = 1.033434062e-1
    case("5f"); exps(1) = 6.952785407e-2
    case("5g"); exps(1) = 8.565417784e-2
    end select
    exps(1) = exps(1) * zeta**2
    coeffs(1) = 1.0
end subroutine

subroutine sto2g(exps, coeffs, zeta, nl_quant)
    ! Taken from: The Journal of Chemical Physics 52, 431 (1970)
    ! http://dx.doi.org/10.1063/1.1672702
    implicit none
    real, intent(in)         :: zeta
    character(2), intent(in) :: nl_quant
    real*8, dimension(2), intent(inout)    :: exps
    real*8, dimension(2), intent(inout)    :: coeffs

    select case(nl_quant)
    case("1s")
        exps = (/8.518186635e-1, 1.516232927e-1/)
        coeffs = (/4.301284983e-1, 6.789135305e-1/)
    case("2s")
        exps = (/1.292278611e-1, 4.908584205e-2/)
        coeffs = (/7.470867124e-1,2.855980556e-1/)
    case("3s")
        exps = (/6.694095822e-1, 5.837135094e-2/)
        coeffs = (/-1.529645716e-1,1.051370110e+0/)
    case("4s")
        exps = (/2.441785453e-1, 4.051097664e-2/)
        coeffs = (/-3.046656896e-1,1.146877294e+0/)
    case("5s")
        exps = (/1.213425654e-1, 3.133152144e-2/)
        coeffs = (/-5.114756049e-1,1.307377277e+0/)
    case("2p")
        exps = (/4.323908358e-1, 1.069139065e-1/)
        coeffs = (/4.522627513e-1,6.713122642e-1/)
    case("3p")
        exps = (/1.458620964e-1, 5.664210742e-2/)
        coeffs = (/5.349653114e-1,5.299607212e-1/)
    case("4p")
        exps = (/6.190052680e-2, 2.648418407e-2/)
        coeffs = (/8.743116767e-1,1.513640107e-1/)
    case("5p")
        exps = (/2.691294191e-1, 3.980805011e-2/)
        coeffs = (/-1.034227010e-1,1.033376378e+0/)
    case("3d")
        exps = (/2.777427345e-1, 8.336507714e-2/)
        coeffs = (/4.666137923e-1,6.644706516e-1/)
    case("4d")
        exps = (/1.330958892e-1, 5.272119659e-2/)
        coeffs = (/4.932764167e-1,5.918727866e-1/)
    case("5d")
        exps = (/6.906014388e-2, 3.399457777e-2/)
        coeffs = (/6.539405185e-1,3.948945302e-1/)
    case("4f")
        exps = (/2.006693538e-1, 6.865384900e-2/)
        coeffs = (/4.769346276e-1,6.587383976e-1/)
    case("5f")
        exps = (/1.156094555e-1, 4.778940916e-2/)
        coeffs = (/4.856637346e-1,6.125980914e-1/)
    case("5g")
        exps = (/1.554531559e-1, 5.854079811e-2/)
        coeffs = (/4.848298074e-1,6.539381621e-1/)
    end select
    exps = exps * zeta**2
end subroutine

subroutine sto3g(exps, coeffs, zeta, nl_quant)
    ! Taken from: The Journal of Chemical Physics 52, 431 (1970)
    ! http://dx.doi.org/10.1063/1.1672702
    implicit none
    real, intent(in)         :: zeta
    character(2), intent(in) :: nl_quant
    real*8, dimension(3), intent(inout)    :: exps
    real*8, dimension(3), intent(inout)    :: coeffs

    select case(nl_quant)
    case("1s")
        exps = (/2.227660584e+0, 4.057711562e-1, 1.098175104e-1/)
        coeffs = (/1.543289673e-1, 5.353281423e-1, 4.446345422e-1/)
    case("2s")
        exps = (/2.581578398e+0, 1.567622104e-1, 6.018332272e-2/)
        coeffs = (/-5.994474934e-2, 5.960385398e-1, 4.581786291e-1/)
    case("3s")
        exps = (/5.641487709e-1, 6.924421391e-2, 3.269529097e-2/)
        coeffs = (/-1.782577972e-1, 8.612761663e-1, 2.261841969e-1/)
    case("4s")
        exps = (/2.267938753e-1, 4.448178019e-2, 2.195294664e-2/)
        coeffs = (/-3.349048323e-1, 1.056744667e+0, 1.256661680e-1/)
    case("5s")
        exps = (/1.080198458e-1, 4.408119382e-2, 2.610811810e-2/)
        coeffs = (/-6.617401158e-1, 7.467595004e-1, 7.146490945e-1/)
    case("2p")
        exps = (/9.192379002e-1, 2.359194503e-1, 8.009805746e-2/)
        coeffs = (/1.623948553e-1, 5.661708862e-1, 4.223071752e-1/)
    case("3p")
        exps = (/2.692880368e+0, 1.489359592e-1, 5.739585040e-2/)
        coeffs = (/-1.061945788e-2, 5.218564264e-1, 5.450015143e-1/)
    case("4p")
        exps = (/4.859692220e-1, 7.430216918e-2, 3.653340923e-2/)
        coeffs = (/-6.147823411e-2, 6.604172234e-1, 3.932639495e-1/)
    case("5p")
        exps = (/2.127482317e-1, 4.729648620e-2, 2.604865324e-2/)
        coeffs = (/-1.389529695e-1, 8.076691064e-1, 2.726029342e-1/)
    case("3d")
        exps = (/5.229112225e-1, 1.639595876e-1, 6.386630021e-2/)
        coeffs = (/1.686596060e-1, 5.847984817e-1, 4.056779523e-1/)
    case("4d")
        exps = (/1.777717219e-1, 8.040647350e-2, 3.949855551e-2/)
        coeffs = (/2.308552718e-1, 6.042409177e-1, 2.595768926e-1/)
    case("5d")
        exps = (/4.913352950e-1, 7.329090601e-2, 3.594209290e-2/)
        coeffs = (/-2.010175008e-2, 5.899310608e-1, 4.658445960e-1/)
    case("4f")
        exps = (/3.483826963e-1, 1.249380537e-1, 5.349995725e-2/)
        coeffs = (/1.737856685e-1, 5.973380628e-1, 3.929395614e-1/)
    case("5f")
        exps = (/1.649233885e-1, 7.487066646e-2, 3.135787219e-2/)
        coeffs = (/1.909729355e-1, 6.146060459e-1, 3.059611271e-1/)
    case("5g")
        exps = (/2.545432122e-1, 1.006544376e-1, 4.624463922e-2/)
        coeffs = (/1.780980905e-1, 6.063757846e-1, 3.828552923e-1/)
    end select
    exps = exps * zeta**2
end subroutine

function nl_to_l(nl_quant) result(l)
    implicit none
    integer :: l
    character(2), intent(in) :: nl_quant
    character(1) :: snd
    l = 0

    snd = nl_quant(2:)
    select case(snd)
    case("s"); l = 0
    case("p"); l = 1
    case("d"); l = 2
    case("f"); l = 3
    case("g"); l = 4
    end select
end function