#!/bin/bash
PROG=approxintdel
DAT=approxintdel_p1.dat
RES=approxintdel_p1.res

cut -d ' ' -f 1-4 sec1sec2s_p1.res > temp1.txt   # e, H, p[2]
cut -d ' ' -f 7-8 hyperdels_p1.res > temp2.txt   # v_s[2]
cut -d ' ' -f 4 hyperdels_p1.res > temp3.txt     # rho_s

# Header: mu, k, stable flag, axis line $G=a$
echo "0.95387536e-3 3 1 0.0" > $DAT

paste -d ' ' temp1.txt temp2.txt temp3.txt >> $DAT

rm temp1.txt temp2.txt temp3.txt
