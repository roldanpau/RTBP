set out "trtbp_unstmfld_H\-1.4514.eps"
set term post enhanced color eps

set size square
f(x)=0.0
plot "trtbp_unstmfld_H-1.4514_p1_br2.res" u 2:3 w l, \
f(x) w l lt 3 lc 2 not

set term pop
set out
