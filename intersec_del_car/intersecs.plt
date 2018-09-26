#set out "intersecs.bak2.eps"
#set term postscript enhanced monochrome eps

set key top left horizontal

set xlabel "H"
set ylabel "L"


#"../approxint_del_car/approxints_unst_SECg_br1.res" u 1:6 w l t "z1p", \
#"../approxint_del_car/approxints_unst_SECg_br2.res" u 1:6 w l t "z1m", \
#"../approxint_del_car/approxints_unst_SECg2_br1.res.bak" u 1:6 w l t "z2p", \
#"../approxint_del_car/approxints_unst_SECg2_br2.res.bak.2" u 1:6 w l t "z2m", \

#plot [:-1.4894] \

plot \
"intersecs_unst_SECg_br1.res" u 1:6 w l t "z1p", \
"intersecs_unst_SECg_br2.res" u 1:6 w l t "z1m", \
"intersecs_unst_SECg2_br1.res" u 1:6 w lp t "z2p", \
"../approxint_del_car/approxints_unst_SECg2_br1.res" u 1:6 w l t "z2p", \
"intersecs_unst_SECg2_br2.res" u 1:6 w lp t "z2m"

set term pop
set out 
