set term tikz
set out "omega_neg_test.tex"

set xlabel "$N$"
set ylabel "$\\abs{\\alpha_1^+(N) - \\alpha_1^+(N-1)}$" noenhanced
plot "omega_neg_test_2.dat" u 1:(abs($2)) w lp not
#t "Convergence of the integral as a function of N"
set out
set term pop
