echo "0.95387536e-3 SEC2 2" >prtbpdel_zf.dat
cut -d ' ' -f 2-5 ../cardel/cardels_2d_zf.res >>prtbpdel_zf.dat
prtbpdel <prtbpdel_zf.dat >prtbpdel_zf.res
rm prtbpdel_zf.dat
