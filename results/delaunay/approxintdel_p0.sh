#!/bin/bash
PROG=approxintdel
DAT=approxintdel_p0.dat
RES=approxintdel_p0.res

cut -d ' ' -f 1-4 sec1sec2s_p0.res > temp1.txt      # e, H, p[2]
cut -d ' ' -f 5-6 hyperdels_p0.res > temp2.txt      # v_u[2]
cut -d ' ' -f 3 hyperdels_p0.res > temp3.txt        # rho_u

# Header: mu, k, stable flag, axis line $G=a$
echo "0.95387536e-3 3 0 0.0" > $DAT

paste -d ' ' temp1.txt temp2.txt temp3.txt >> $DAT

rm temp1.txt temp2.txt temp3.txt
