echo "0.95387536e-3 SEC1 2" >prtbpdel_p2.dat
cut -d ' ' -f 2-5 ../cardel/cardels.res >>prtbpdel_p2.dat
prtbpdel <prtbpdel_p2.dat >prtbpdel_p2.res
rm prtbpdel_p2.dat
