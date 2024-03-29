#set out "intersecs.eps"
#set term postscript enhanced color eps

set key top left horizontal spacing 4

set xlabel "H"
set ylabel "L"


#"../approxint_del_car/approxints_unst_SECg_br1.res" u 1:6 w l t "z1p", \
#"../approxint_del_car/approxints_unst_SECg_br2.res" u 1:6 w l t "z1m", \
#"../approxint_del_car/approxints_unst_SECg2_br1.res" u 1:6 w l t "z2p", \
#"../approxint_del_car/approxints_unst_SECg2_br2.res" u 1:6 w l t "z2m", \

#plot [:-1.4894] \

plot \
"intersecs_unst_SECg_br1.res" u 1:6 w l t "z1p", \
"intersecs_unst_SECg_br2.res" u 1:6 w l t "z1m", \
"intersecs_unst_SECg2_br1.res" u 1:6 w l t "z2p", \
"intersecs_unst_SECg2_br2.res" u 1:6 w l t "z2m"

set term pop
set out 
