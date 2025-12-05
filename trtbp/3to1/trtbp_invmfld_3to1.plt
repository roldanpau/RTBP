# 3:1 periodic orbit of the R3BP in ROTATING coordinates, corresponding to
# eccentricity e=0.8, and a trajectory on its invariant unstable manifold.

set xlabel "$x$"
set ylabel "$y$"
set size square

set term tikz
set out "trtbp_invmfld_3to1.tex"

plot \
"trtbp_3to1_ell08.res" using 2:3 w l lw 3 not, \
"trtbp_3to1_ell08_invmfld.res" using 2:3 w l dt 2 not

#"L1.dat" t "L1"

#"trtbp_3to1_ell08.res" using 2:3 w l lc 1 t "e=0.8", \

set out
set term pop
