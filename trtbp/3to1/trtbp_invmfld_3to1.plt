# 3:1 periodic orbit of the R3BP in ROTATING coordinates, corresponding to
# eccentricity e=0.8, and a trajectory on its invariant unstable manifold.

#set out "trtbp_invmfld_3to1.eps"
#set term post eps

set size square

plot \
"trtbp_3to1_ell07_invmfld.res" using 2:3 w l lc 2 t "trajectory in inv. mfld.", \
"L1.dat" t "L1"

#"trtbp_3to1_ell08.res" using 2:3 w l lc 1 t "e=0.8", \
#set term pop
#set out
