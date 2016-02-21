# input: H \omega_+ \omega_in

set xlabel "H" 
unset ylabel

set key t l

# remember \omega_out = \omega_+ - \omega_- = 2\omega_+
# and \omega = \omega_out + \omega_in

# remember that, since Marcel changed signs in the last version of the paper,
# our computed value of \omega_+ is now the negative of Marcel's \omega_+.

#set out "omega_fb.eps"
#set term post eps

# plot omega_f, omega_b

plot [-1.81:-1.56]\
"./inner_sep/omega_f.dat" u 1:(-$2) w l t "wposf", \
"./inner_sep/omega_f.dat" u 1:($3) w l t "winf", \
"./inner_sep/omega_f.res" w l t "wf", \
"./omega_b.dat" u 1:(-$2) w l t "wposb", \
"./omega_b.dat" u 1:($3) w l t "winb", \
"./omega_b.res" w l t "wb"

#set term pop
#set out
