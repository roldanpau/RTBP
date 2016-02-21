# Try to compute periodic orbit for various energy levels (specified in input
# file energies.dat). Leave result in output file energies.res.

for i in `cat energies.dat`; do echo "0.95387536e-3 $i 1.302533194 0.0 1.123473908" | portbp >>energies.res; done
