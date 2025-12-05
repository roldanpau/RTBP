# Family of resonant periodic orbits in rotating coordinates

set xlabel "$x$"
set ylabel "$y$"
set size square

set term tikz

set style line 1 lc rgb 'yellow' pt 7 ps 2  # circle
set style line 2 lc rgb 'brown' pt 7 ps 2  # circle

do for [e=0:8] {
    pltfile=sprintf("trtbp_3to1_ell0%d.res",e)
    outfile=sprintf("trtbp_3to1_ell0%d.tex",e)
    set output outfile
    plot [-1.1:1.1] [-1:1] \
    pltfile using 2:3 w l not, \
    "<echo '-0.001 0'" with points ls 1 not, \
    "<echo '0.999 0'" with points ls 2 not

#    "L1.dat" not
    set out
}

set term pop
