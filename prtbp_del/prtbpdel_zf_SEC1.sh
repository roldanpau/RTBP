echo "0.95387536e-3 SEC1 1" >prtbpdel_zf.dat
cut -d ' ' -f 2-5 ../cardel/cardels_2d_zf.res >>prtbpdel_zf.dat
prtbpdel <prtbpdel_zf.dat >prtbpdel_zf_SEC1.res
rm prtbpdel_zf.dat
