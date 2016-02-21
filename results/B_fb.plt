# input:
# H, \Im(B^+), \Re(B_{in}^f), \Im(B_{in}^f), \omega_pos
# H, \Im(B^+), \Re(B_{in}^b), \Im(B_{in}^b), \omega_{in}^b

set xlabel "H" 
unset ylabel

set key t l spacing 3.0

# remember \omega_out = \omega_+ - \omega_- = 2\omega_+

# remember that, since Marcel changed signs in the last version of the paper,
# our computed value of \omega_+ is now the negative of Marcel's \omega_+.

#set out "B_fb.eps"
#set term post eps color

# plot Re(B^f), Im(B^f), Re(B^b), Im(B^b)
#plot [-1.81:-1.56] \
#"./inner_sep/B_f.dat" u 1:($3*cos(-0.00095387536*2.0*$5) - $4*sin(-0.00095387536*2.0*$5)) w l t "reBf", \
#"./inner_sep/B_f.dat" u 1:(2.0*0.00095387536*$2 + $4*cos(-0.00095387536*2.0*$5) + $3*sin(-0.00095387536*2.0*$5)) w l t "imBf", \
#"./B_b.dat" u 1:($3 - 2.0*0.00095387536*$2*sin(0.00095387536*$5)) w l t "reBb", \
#"./B_b.dat" u 1:($4 + 2.0*0.00095387536*$2*cos(0.00095387536*$5)) w l t "imBb"

plot [-1.81:-1.56] \
"./inner_sep/B_f.res" u 1:2 w l t "reBfXXXXXX", \
"./inner_sep/B_f.res" u 1:3 w l t "imBfXXXXXX", \
"./B_b.res" u 1:2 w l t "reBbXXXXXX", \
"./B_b.res" u 1:3 w l t "imBbXXXXXX"


#set term pop
#set out
