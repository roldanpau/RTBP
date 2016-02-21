set out "homoclinic.tex"
set term epslatex color
set title "Homoclinic orbit"
set xlabel "$x$"
set ylabel "$p_x$"

plot "unstmfld.res" w dots linewidth 5 notitle, \
	"stmfld.res" w dots linewidth 5 notitle, \
	"homoclinic.res"
set term pop
set out
#pause -1

set out "homoclinic_zoom.tex"
set term epslatex color
set title "Homoclinic orbit (zoom)"
set xlabel "$x$"
set ylabel "$p_x$"

plot [1.95:] [:-0.39] "unstmfld.res" ps 0.5 notitle, \
	"stmfld.res" ps 0.5 notitle, \
	"homoclinic.res"
set term pop
set out
#pause -1

