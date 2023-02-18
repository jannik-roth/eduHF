cd eduHF
python3 -m numpy.f2py -c slater_expansion.f90 -m slater_expansion
python3 -m numpy.f2py -c mcmurchie_davidson.f90 -m mcmurchie_davidson
cd ..