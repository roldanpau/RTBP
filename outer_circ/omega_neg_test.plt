set term eps 
set out "omega_neg_test.eps"
set xlabel "N"
set ylabel "difference |int_0^{N} - int_0^{N-1}|"
plot "omega_neg_test_2.dat" u 1:(abs($2)) w lp not
#t "Convergence of the integral as a function of N"
