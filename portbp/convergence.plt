set out "convergence.tex"
set term epslatex color
set title "Convergence of Newton-like method"
set xlabel "$d_{k}$"
set ylabel "$d_{k+1}$"

set log xy
set xrange [1.e-14:1]
set yrange [1.e-14:1]                                                                    
set size square
plot "portbp.res.paste" using \
(sqrt(($1)*($1)+($2)*($2))):(sqrt(($3)*($3)+($4)*($4))) notitle
set term pop
set out
