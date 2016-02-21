echo "0.95387536e-3 SEC1 1" >prtbpdel_p1.dat
cut -d ' ' -f 2-5 ../cardel/cardels.res >>prtbpdel_p1.dat
prtbpdel <prtbpdel_p1.dat >prtbpdel_p1.res
rm prtbpdel_p1.dat
