#!/bin/bash

# Intersection of unst mfld of p0 with vertical axis.

# This should give exactly the same result as intersec_p1.sh, i.e.
# intersection of st mfld of p1 with vertical axis.

PROG=intersecdel
DAT=intersecdel_p0.dat
RES=intersecdel_p0.res

cut -d ' ' -f 1-4 sec1sec2s_p0.res > temp1.txt      # e, H, p[2]
cut -d ' ' -f 5-6 hyperdels_p0.res > temp2.txt      # v_u[2]
cut -d ' ' -f 3 hyperdels_p0.res > temp3.txt        # rho_u
cut -d ' ' -f 3-5 approxintdel_p0.res > temp4.txt   # n, h_1, h_2

# Header: mu, stable flag, axis line $G=a$
echo "0.95387536e-3 0 0.0" > $DAT

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt >> $DAT

rm temp*.txt
