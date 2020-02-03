#set out "approxints.eps"
#set term eps 

set xlabel "H"
set ylabel "x"

set key right center

plot "approxints_unst_br1.res" u 1:5 w l, \
"approxints_unst_br2.res" u 1:5 w l, \
"approxints_st_br1.res" u 1:5 w l, \
"approxints_st_br2.res" u 1:5 w l

set term pop
set out 
