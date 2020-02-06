set key center right

plot \
	"intersecs_unst_br1.res" u 1:5 w l, \
	"intersecs_unst_br2.res" u 1:5 w l, \
    "intersecs_st_br1.res" u 1:5 w l, \
    "intersecs_st_br2.res" u 1:5 w l
