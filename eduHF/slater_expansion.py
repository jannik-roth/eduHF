import numpy as np
from scipy.special import factorial2

def slater_expansion(zeta : float, nl_quant : str, ng : int) -> tuple[np.array, np.array]:
    match ng:
        case 1:
            coeffs, alphas = sto1g(zeta, nl_quant)
        case 2:
            coeffs, alphas = sto2g(zeta, nl_quant)
        case 3:
            coeffs, alphas = sto3g(zeta, nl_quant)
        case 4:
            coeffs, alphas = sto4g(zeta, nl_quant)
        case 5:
            coeffs, alphas = sto5g(zeta, nl_quant)
        case 6:
            coeffs, alphas = sto6g(zeta, nl_quant)

    l = _letter_to_l(nl_quant[1])
    coeffs = np.multiply(coeffs, 
                         (2.0 * alphas/ np.pi) ** 0.75 * np.sqrt(4.0 * alphas)**(l / np.sqrt(factorial2(2*l-1))))

    return coeffs, alphas

def sto1g(zeta : float, nl_quant : str) -> tuple[np.array, np.array]:
    # Taken from: The Journal of Chemical Physics 52, 431 (1970)
    # http://dx.doi.org/10.1063/1.1672702
    #alphas = np.zeros(1)
    match nl_quant:
        case("1s"): 
            alphas = np.array([2.709498091e-1])
        case("2s"):
            alphas = np.array([1.012151084e-1])
        case("3s"):
            alphas = np.array([5.296881757e-2])
        case("4s"):
            alphas = np.array([3.264600274e-2])
        case("5s"):
            alphas = np.array([2.216912938e-2])
        case("2p"):
            alphas = np.array([1.759666885e-1])
        case("3p"):
            alphas = np.array([9.113614253e-2])
        case("4p"):
            alphas = np.array([5.578350235e-2])
        case("5p"):
            alphas = np.array([3.769845216e-2])
        case("3d"):
            alphas = np.array([1.302270363e-1])
        case("4d"):
            alphas = np.array([7.941656339e-2])
        case("5d"):
            alphas = np.array([5.352200793e-2])
        case("4f"):
            alphas = np.array([1.033434062e-1])
        case("5f"):
            alphas = np.array([6.952785407e-2])
        case("5g"):
            alphas = np.array([8.565417784e-2])
    alphas *= zeta ** 2
    coeffs = np.array([1.0])

    return coeffs, alphas

def sto2g(zeta : float, nl_quant : str) -> tuple[np.array, np.array]:
    # Taken from: The Journal of Chemical Physics 52, 431 (1970)
    # http://dx.doi.org/10.1063/1.1672702
    #alphas = np.zeros(1)
    match nl_quant:
        case("1s"):
            alphas = np.array([8.518186635e-1, 1.516232927e-1])
            coeffs = np.array([4.301284983e-1, 6.789135305e-1])
        case("2s"):
            alphas = np.array([1.292278611e-1, 4.908584205e-2])
            coeffs = np.array([7.470867124e-1,2.855980556e-1])
        case("3s"):
            alphas = np.array([6.694095822e-1, 5.837135094e-2])
            coeffs = np.array([-1.529645716e-1,1.051370110e+0])
        case("4s"):
            alphas = np.array([2.441785453e-1, 4.051097664e-2])
            coeffs = np.array([-3.046656896e-1,1.146877294e+0])
        case("5s"):
            alphas = np.array([1.213425654e-1, 3.133152144e-2])
            coeffs = np.array([-5.114756049e-1,1.307377277e+0])
        case("2p"):
            alphas = np.array([4.323908358e-1, 1.069139065e-1])
            coeffs = np.array([4.522627513e-1,6.713122642e-1])
        case("3p"):
            alphas = np.array([1.458620964e-1, 5.664210742e-2])
            coeffs = np.array([5.349653114e-1,5.299607212e-1])
        case("4p"):
            alphas = np.array([6.190052680e-2, 2.648418407e-2])
            coeffs = np.array([8.743116767e-1,1.513640107e-1])
        case("5p"):
            alphas = np.array([2.691294191e-1, 3.980805011e-2])
            coeffs = np.array([-1.034227010e-1,1.033376378e+0])
        case("3d"):
            alphas = np.array([2.777427345e-1, 8.336507714e-2])
            coeffs = np.array([4.666137923e-1,6.644706516e-1])
        case("4d"):
            alphas = np.array([1.330958892e-1, 5.272119659e-2])
            coeffs = np.array([4.932764167e-1,5.918727866e-1])
        case("5d"):
            alphas = np.array([6.906014388e-2, 3.399457777e-2])
            coeffs = np.array([6.539405185e-1,3.948945302e-1])
        case("4f"):
            alphas = np.array([2.006693538e-1, 6.865384900e-2])
            coeffs = np.array([4.769346276e-1,6.587383976e-1])
        case("5f"):
            alphas = np.array([1.156094555e-1, 4.778940916e-2])
            coeffs = np.array([4.856637346e-1,6.125980914e-1])
        case("5g"):
            alphas = np.array([1.554531559e-1, 5.854079811e-2])
            coeffs = np.array([4.848298074e-1,6.539381621e-1])
    alphas *= zeta ** 2

    return coeffs, alphas

def sto3g(zeta : float, nl_quant : str) -> tuple[np.array, np.array]:
    # Taken from: The Journal of Chemical Physics 52, 431 (1970)
    # http://dx.doi.org/10.1063/1.1672702
    #alphas = np.zeros(1)
    match nl_quant:
        case("1s"):
            alphas = np.array([2.227660584e+0, 4.057711562e-1, 1.098175104e-1])
            coeffs = np.array([1.543289673e-1, 5.353281423e-1, 4.446345422e-1])
        case("2s"):
            alphas = np.array([2.581578398e+0, 1.567622104e-1, 6.018332272e-2])
            coeffs = np.array([-5.994474934e-2, 5.960385398e-1, 4.581786291e-1])
        case("3s"):
            alphas = np.array([5.641487709e-1, 6.924421391e-2, 3.269529097e-2])
            coeffs = np.array([-1.782577972e-1, 8.612761663e-1, 2.261841969e-1])
        case("4s"):
            alphas = np.array([2.267938753e-1, 4.448178019e-2, 2.195294664e-2])
            coeffs = np.array([-3.349048323e-1, 1.056744667e+0, 1.256661680e-1])
        case("5s"):
            alphas = np.array([1.080198458e-1, 4.408119382e-2, 2.610811810e-2])
            coeffs = np.array([-6.617401158e-1, 7.467595004e-1, 7.146490945e-1])
        case("2p"):
            alphas = np.array([9.192379002e-1, 2.359194503e-1, 8.009805746e-2])
            coeffs = np.array([1.623948553e-1, 5.661708862e-1, 4.223071752e-1])
        case("3p"):
            alphas = np.array([2.692880368e+0, 1.489359592e-1, 5.739585040e-2])
            coeffs = np.array([-1.061945788e-2, 5.218564264e-1, 5.450015143e-1])
        case("4p"):
            alphas = np.array([4.859692220e-1, 7.430216918e-2, 3.653340923e-2])
            coeffs = np.array([-6.147823411e-2, 6.604172234e-1, 3.932639495e-1])
        case("5p"):
            alphas = np.array([2.127482317e-1, 4.729648620e-2, 2.604865324e-2])
            coeffs = np.array([-1.389529695e-1, 8.076691064e-1, 2.726029342e-1])
        case("3d"):
            alphas = np.array([5.229112225e-1, 1.639595876e-1, 6.386630021e-2])
            coeffs = np.array([1.686596060e-1, 5.847984817e-1, 4.056779523e-1])
        case("4d"):
            alphas = np.array([1.777717219e-1, 8.040647350e-2, 3.949855551e-2])
            coeffs = np.array([2.308552718e-1, 6.042409177e-1, 2.595768926e-1])
        case("5d"):
            alphas = np.array([4.913352950e-1, 7.329090601e-2, 3.594209290e-2])
            coeffs = np.array([-2.010175008e-2, 5.899310608e-1, 4.658445960e-1])
        case("4f"):
            alphas = np.array([3.483826963e-1, 1.249380537e-1, 5.349995725e-2])
            coeffs = np.array([1.737856685e-1, 5.973380628e-1, 3.929395614e-1])
        case("5f"):
            alphas = np.array([1.649233885e-1, 7.487066646e-2, 3.135787219e-2])
            coeffs = np.array([1.909729355e-1, 6.146060459e-1, 3.059611271e-1])
        case("5g"):
            alphas = np.array([2.545432122e-1, 1.006544376e-1, 4.624463922e-2])
            coeffs = np.array([1.780980905e-1, 6.063757846e-1, 3.828552923e-1])
    alphas *= zeta ** 2

    return coeffs, alphas

def sto4g(zeta : float, nl_quant : str) -> tuple[np.array, np.array]:
    # Taken from: The Journal of Chemical Physics 52, 431 (1970)
    # http://dx.doi.org/10.1063/1.1672702
    match nl_quant:
        case("1s"):
            alphas = np.array([5.216844534e+0, 9.546182760e-1, 2.652034102e-1, 8.801862774e-2])
            coeffs = np.array([5.675242080e-2, 2.601413550e-1, 5.328461143e-1, 2.916254405e-1])
        case("2s"):
            alphas = np.array([1.161525551e+1, 2.000243111e+0, 1.607280687e-1, 6.125744532e-2])
            coeffs = np.array([-1.198411747e-2,-5.472052539e-2, 5.805587176e-1, 4.770079976e-1])
        case("3s"):
            alphas = np.array([1.513265591e+0, 4.262497508e-1, 7.643320863e-2, 3.760545063e-2])
            coeffs = np.array([-3.295496352e-2,-1.724516959e-1, 7.518511194e-1, 3.589627317e-1])
        case("4s"):
            alphas = np.array([3.242212833e-1, 1.663217177e-1, 5.081097451e-2, 2.829066600e-2])
            coeffs = np.array([-1.120682822e-1,-2.845426863e-1, 8.909873788e-1, 3.517811205e-1])
        case("5s"):
            alphas = np.array([8.602284252e-1, 1.189050200e-1, 3.446076176e-2, 1.974798796e-2])
            coeffs = np.array([1.103657561e-2,-5.606519023e-1, 1.179429987e+0, 1.734974376e-1])
        case("2p"):
            alphas = np.array([1.798260992e+0, 4.662622228e-1, 1.643718620e-1, 6.543927065e-2])
            coeffs = np.array([5.713170255e-2, 2.857455515e-1, 5.517873105e-1, 2.632314924e-1])
        case("3p"):
            alphas = np.array([1.853180239e+0, 1.915075719e-1, 8.655487938e-2, 4.184253862e-2])
            coeffs = np.array([-1.434249391e-2, 2.755177589e-1, 5.846750879e-1, 2.144986514e-1])
        case("4p"):
            alphas = np.array([1.492607880e+0, 4.327619272e-1, 7.553156064e-2, 3.706272183e-2])
            coeffs = np.array([-6.035216774e-3,-6.013310874e-2, 6.451518200e-1, 4.117923820e-1])
        case("5p"):
            alphas = np.array([3.962838833e-1, 1.838858552e-1, 4.943555157e-2, 2.750222273e-2])
            coeffs = np.array([-1.801459207e-2,-1.360777372e-1, 7.533973719e-1, 3.409304859e-1])
        case("3d"):
            alphas = np.array([9.185846715e-1, 2.920461109e-1, 1.187568890e-1, 5.286755896e-2])
            coeffs = np.array([5.799057705e-2, 3.045581349e-1, 5.601358038e-1, 2.432423313e-1])
        case("4d"):
            alphas = np.array([1.995825422e+0, 1.823461280e-1, 8.197240896e-2, 4.000634951e-2])
            coeffs = np.array([-2.816702620e-3, 2.177095871e-1, 6.058047348e-1, 2.717811257e-1])
        case("5d"):
            alphas = np.array([4.230617826e-1, 8.293863702e-2, 4.590326388e-2, 2.628744797e-2])
            coeffs = np.array([-2.421626009e-2, 3.937644956e-1, 5.489520286e-1, 1.190436963e-1])
        case("4f"):
            alphas = np.array([5.691670217e-1, 2.074585819e-1, 9.298346885e-2, 4.473508853e-2])
            coeffs = np.array([5.902730589e-2, 3.191828952e-1, 5.639423893e-1, 2.284796537e-1])
        case("5f"):
            alphas = np.array([2.017831152e-1, 1.001952178e-1, 5.441006630e-2, 3.037569283e-2])
            coeffs = np.array([9.174268830e-2, 4.023496947e-1, 4.937432100e-1, 1.254001522e-1])
        case("5g"):
            alphas = np.array([3.945205573e-1, 1.588100623e-1, 7.646521729e-1, 3.898703611e-2])
            coeffs = np.array([6.010484250e-2, 3.309738329e-1, 5.655207585e-1, 2.171122608e-1])
    alphas *= zeta ** 2

    return coeffs, alphas

def sto5g(zeta : float, nl_quant : str) -> tuple[np.array, np.array]:
    # Taken from: The Journal of Chemical Physics 52, 431 (1970)
    # http://dx.doi.org/10.1063/1.1672702
    match nl_quant:
        case("1s"):
            alphas = np.array([1.130563696e+1, 2.071728178e+0, 5.786484833e-1, 1.975724573e-1, 7.445271746e-2])
            coeffs = np.array([2.214055312e-2, 1.135411520e-1, 3.318161484e-1, 4.825700713e-1, 1.935721966e-1])
        case("2s"):
            alphas = np.array([8.984956862e+0, 1.673710636e+0, 1.944726668e-1, 8.806345634e-2, 4.249068522e-2])
            coeffs = np.array([-1.596349096e-2,-5.685884883e-2, 3.698265599e-1, 5.480512593e-1, 1.472634893e-1])
        case("3s"):
            alphas = np.array([4.275877914e+0, 1.132409433e+0, 4.016256968e-1, 7.732370620e-2, 3.800708627e-2])
            coeffs = np.array([-3.920358850e-3,-4.168430506e-2,-1.637440990e-1, 7.419373723e-1, 3.724364929e-1])
        case("4s"):
            alphas = np.array([2.980263783e+0, 3.792228833e-1, 1.789717224e-1, 5.002110360e-2, 2.789361681e-2])
            coeffs = np.array([1.513948997e-3,-7.316801518e-2,-3.143703799e-1, 9.032615169e-1, 3.294210848e-1])
        case("5s"):
            alphas = np.array([7.403763257e-1, 1.367990863e-1, 9.135301779e-2, 3.726907315e-2, 2.241490836e-2])
            coeffs = np.array([1.375523371e-2,-3.097344179e-1,-3.199192259e-1, 1.084547038e+0, 3.345288361e-1])
        case("2p"):
            alphas = np.array([3.320386533e+0, 8.643257633e-1, 3.079819284e-1, 1.273309895e-1, 5.606243164e-2])
            coeffs = np.array([2.079051117e-2, 1.235472099e-1, 3.667738986e-1, 4.834930290e-1, 1.653444074e-1])
        case("3p"):
            alphas = np.array([6.466803859e+0, 1.555914802e+0, 1.955925255e-1, 8.809647701e-2, 4.234835707e-2])
            coeffs = np.array([-2.329023747e-3,-1.357395221e-2, 2.632185383e-1, 5.880427024e-1, 2.242794445e-1])
        case("4p"):
            alphas = np.array([1.091977298e+0, 3.719985051e-1, 8.590019352e-2, 4.786503860e-2, 2.730479990e-2])
            coeffs = np.array([-1.143929558e-2,-6.322651538e-2, 4.398907721e-1, 5.245859166e-1, 1.017072253e-1])
        case("5p"):
            alphas = np.array([3.422168934e-1, 1.665099900e-1, 5.443732013e-2, 3.367775277e-2, 2.091949042e-2])
            coeffs = np.array([-3.113958289e-2,-1.374007017e-1, 5.573881018e-1, 4.855428100e-1, 6.605423564e-2])
        case("3d"):
            alphas = np.array([1.539033958e+0, 4.922090297e-1, 2.029756928e-1, 9.424112917e-2, 4.569058269e-2])
            coeffs = np.array([2.020869128e-2, 1.321157923e-1, 3.911240346e-1, 4.779609701e-1, 1.463662294e-1])
        case("4d"):
            alphas = np.array([1.522122079e+0, 2.173041823e-1, 1.084876577e-1, 5.836797641e-2, 3.206682246e-2])
            coeffs = np.array([-3.673711876e-3, 1.167122499e-1, 4.216476416e-1, 4.547673415e-1, 1.037803318e-1])
        case("5d"):
            alphas = np.array([9.702946470e-1, 3.603270196e-1, 8.668717752e-2, 4.833708379e-2, 2.751899341e-2])
            coeffs = np.array([-3.231527611e-3,-2.434931372e-2, 3.440817054e-1, 5.693674316e-1, 1.511340183e-1])
        case("4f"):
            alphas = np.array([8.925960415e-1, 3.277589120e-1, 1.492869962e-1, 7.506099109e-2, 3.892475795e-2])
            coeffs = np.array([1.999839052e-2, 1.395427440e-1, 4.091508237e-1, 4.708252119e-1, 1.328082566e-1])
        case("5f"):
            alphas = np.array([1.670735676e+0, 2.072477219e-1, 1.024709357e-1, 5.531913898e-2, 3.072866652e-2])
            coeffs = np.array([-7.301193568e-4, 8.414991343e-2, 3.923683153e-1, 5.040033146e-1, 1.328979300e-1])
        case("5g"):
            alphas = np.array([5.895429375e-1, 2.393343780e-1, 1.172646904e-1, 6.254074479e-2, 3.411243214e-2])
            coeffs = np.array([1.998085812e-2, 1.460384050e-1, 4.230565459e-1, 4.635699665e-1, 1.226411691e-1])
    alphas *= zeta ** 2

    return coeffs, alphas

def sto6g(zeta : float, nl_quant : str) -> tuple[np.array, np.array]:
    # Taken from: The Journal of Chemical Physics 52, 431 (1970)
    # http://dx.doi.org/10.1063/1.1672702
    match nl_quant:
        case("1s"):
            alphas = np.array([2.310303149e+1, 4.235915534e+0, 1.185056519e+0, 4.070988982e-1, 1.580884151e-1, 6.510953954e-2])
            coeffs = np.array([9.163596280e-3, 4.936149294e-2, 1.685383049e-1, 3.705627997e-1, 4.164915298e-1, 1.303340841e-1])
        case("2s"):
            alphas = np.array([2.768496241e+1, 5.077140627e+0, 1.426786050e+0, 2.040335729e-1, 9.260298399e-2, 4.416183978e-2])
            coeffs = np.array([-4.151277819e-3,-2.067024148e-2,-5.150303337e-2, 3.346271174e-1, 5.621061301e-1, 1.712994697e-1])
        case("3s"):
            alphas = np.array([3.273031938e+0, 9.200611311e-1, 3.593349765e-1, 8.636686991e-2, 4.797373812e-2, 2.724741144e-2])
            coeffs = np.array([-6.775596947e-3,-5.639325779e-2,-1.587856086e-1, 5.534527651e-1, 5.015351020e-1, 7.223633674e-2])
        case("4s"):
            alphas = np.array([3.232838646e+0, 3.605788802e-1, 1.717905487e-1, 5.277666487e-2, 3.163400284e-2, 1.874093091e-2])
            coeffs = np.array([1.374817488e-3,-8.666390043e-2,-3.130627309e-1, 7.812787397e-1, 4.389247988e-1, 2.487178756e-2])
        case("5s"):
            alphas = np.array([1.410128298e+0, 5.077878915e-1, 1.847926858e-1, 1.061070594e-1, 3.669584901e-2, 2.213558430e-2])
            coeffs = np.array([2.695439582e-3, 1.850157487e-2,-9.588628125e-2, -5.200673560e-1, 1.087619490e+0, 3.103964343e-1])
        case("2p"):
            alphas = np.array([5.868285913e+0, 1.530329631e+0, 5.475665231e-1, 2.288932733e-1, 1.046655969e-1, 4.948220127e-2])
            coeffs = np.array([7.924233646e-3, 5.144104825e-2, 1.898400060e-1, 4.049863191e-1, 4.012362861e-1, 1.051855189e-1])
        case("3p"):
            alphas = np.array([5.077973607e+0, 1.340786940e+0, 2.248434849e-1, 1.131741848e-1, 6.076408893e-2, 3.315424265e-2])
            coeffs = np.array([-3.329929840e-3,-1.419488340e-2, 1.639395770e-1, 4.485358256e-1, 3.908813050e-1, 7.411456232e-2])
        case("4p"):
            alphas = np.array([2.389722618e+0, 7.960947826e-1, 3.415541380e-1, 8.847434525e-2, 4.958248334e-2, 2.816929784e-2])
            coeffs = np.array([-1.665913575e-3,-1.657464971e-2,-5.958513378e-2, 4.053115554e-1, 5.433958189e-1, 1.204970491e-1])
        case("5p"):
            alphas = np.array([3.778623374e+0, 3.499121109e-1, 1.683175469e-1, 5.404070736e-2, 3.328911801e-2, 2.063815019e-2])
            coeffs = np.array([1.163246387e-4,-2.920771322e-2,-1.381051233e-1, 5.706134877e-1, 4.768808140e-1, 6.021665516e-2])
        case("3d"):
            alphas = np.array([2.488296923e+0, 7.981487853e-1, 3.311327490e-1, 1.559114463e-1, 7.817734732e-2, 4.058484363e-2])
            coeffs = np.array([7.283828112e-3, 5.386799363e-2, 2.072139149e-1, 4.266269092e-1, 3.843100204e-1, 8.902827546e-2])
        case("4d"):
            alphas = np.array([4.634239420e+0, 1.341648295e+0, 2.209593028e-1, 1.101467943e-1, 5.904190370e-2, 3.232628887e-2])
            coeffs = np.array([-4.749842876e-4,-3.566777891e-3, 1.108670481e-1, 4.159646930e-1, 4.621672517e-1, 1.081250196e-1])
        case("5d"):
            alphas = np.array([8.820520428e-1, 3.410838409e-1, 9.204308840e-2, 5.472831774e-2, 3.391202830e-2, 2.108227374e-2])
            coeffs = np.array([-4.097311019e-3,-2.508271857e-2, 2.648458555e-1, 5.097437054e-1, 2.654483467e-1, 2.623132212e-2])
        case("4f"):
            alphas = np.array([1.357718039e+0, 5.004907278e-1, 2.296565064e-1, 1.173146814e-1, 6.350097171e-2, 3.474556673e-2])
            coeffs = np.array([6.930234381e-3, 5.634263745e-2, 2.217065797e-1, 4.411388883e-1, 3.688112625e-1, 7.787514504e-2])
        case("5f"):
            alphas = np.array([1.334096840e+0, 2.372312347e-1, 1.269485144e-1, 7.290318381e-2, 4.351355997e-2, 2.598071843e-2])
            coeffs = np.array([-9.486751531e-4, 4.624275998e-2, 2.373699784e-1, 4.589112231e-1, 3.205010548e-1, 5.077063693e-2])
        case("5g"):
            alphas = np.array([8.574668996e-1, 3.497184772e-1, 1.727917060e-1,  9.373643151e-2, 5.340032759e-2, 3.051364464e-2])
            coeffs = np.array([6.729778096e-3, 5.874145170e-2, 2.339955227e-1, 4.512983737e-1, 3.552053926e-1, 6.974153145e-2])
    alphas *= zeta ** 2

    return coeffs, alphas

def _letter_to_l(letter):
    match letter.lower():
        case 's':
            return 0
        case 'p':
            return 1
        case 'd':
            return 2
        case 'f':
            return 3
        case 'g':
            return 4