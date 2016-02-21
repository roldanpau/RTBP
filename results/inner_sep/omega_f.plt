# input: H \omega_+ \omega_in

set xlabel "H" 
unset ylabel

# remember \omega_out = \omega_+ - \omega_- = 2\omega_+
# and \omega = \omega_out + \omega_in

# remember that, since Marcel changed signs in the last version of the paper,
# our computed value of \omega_+ is now the negative of Marcel's \omega_+.

#set out "omega_f.eps"
#set term post eps

# plot omega_+, omega_in, omega
plot \
"./omega_f.dat" u 1:(-$2) w l t "wpos", \
"./omega_f.dat" u 1:($3) w l t "win", \
"./omega_f.res" u 1:2 w l t "w"

#set term pop
#set out
