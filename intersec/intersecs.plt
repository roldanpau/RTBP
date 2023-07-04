#set out "intersecs.eps"
#set term post eps color

set key center right

set xlabel "H"
set ylabel "x"

plot \
	"intersecs_unst_br1.res" u 1:7 w l, \
	"intersecs_unst_br2.res_nokinks" u 1:7 w l, \
    "intersecs_st_br1.res" u 1:7 w l, \
    "intersecs_st_br2.res_nokinks" u 1:7 w l

#set term pop
#set out
