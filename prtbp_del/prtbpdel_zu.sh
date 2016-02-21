echo "0.95387536e-3 SEC1 1" >prtbpdel_zu.dat
cut -d ' ' -f 2-5 ../cardel/cardels_2d_zu.res >>prtbpdel_zu.dat
prtbpdel <prtbpdel_zu.dat >prtbpdel_zu.res
rm prtbpdel_zu.dat
