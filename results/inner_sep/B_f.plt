# input:
# H, \Im(B^+), \Re(B_{in}^f), \Im(B_{in}^f), \omega_pos

set xlabel "H" 
unset ylabel

# remember \omega_out = \omega_+ - \omega_- = 2\omega_+

# remember that, since Marcel changed signs in the last version of the paper,
# our computed value of \omega_+ is now the negative of Marcel's \omega_+.

#set out "omega_f.eps"
#set term post eps

# plot Re(B^f), Im(B^f)
plot \
"./B_f.res" u 1:2 w l t "re", \
"./B_f.res" u 1:3 w l t "im"

#set term pop
#set out
