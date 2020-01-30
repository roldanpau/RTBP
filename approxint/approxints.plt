#set out "approxints.eps"
#set term eps 

plot "approxints_unst_br1.res" u 1:5, \
"approxints_unst_br2.res" u 1:5

set term pop
set out 
