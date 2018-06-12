#set out "intersecs.eps"
#set term postscript enhanced monochrome eps

set key top left horizontal

set xlabel "H"
set ylabel "L"

plot [:-1.4894] \
"intersecs_unst_SECg_br1.res" u 1:6 w l t "z1p", \
"intersecs_unst_SECg_br2.res" u 1:6 w l t "z1m", \
"intersecs_unst_SECg2_br1.res" u 1:6 w l t "z2p", \
"intersecs_unst_SECg2_br2.res" u 1:6 w l t "z2m"

set term pop
set out 
