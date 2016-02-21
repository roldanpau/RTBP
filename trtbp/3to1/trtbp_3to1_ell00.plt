# 3:1 periodic orbits of the R3BP in ROTATING coordinates, corresponding to
# different eccentricities (from e=0.0 until e=0.85).

set out "trtbp_3to1_ell00.eps"
set term post eps

set size square

set xlabel "x"
set ylabel "y"
plot [-.6:1] [-.8:.8]"trtbp_3to1_ell00.res" using 2:3 w l t "H=-1.733, C=3.467"

set term pop
set out
