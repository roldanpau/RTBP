#!/bin/bash
PROJ=/home/roldan/research/splitting
PROG=sec1sec2del_inv
DAT=sec1sec2s.dat
RES=sec1sec2s_p1.res

if [ -f $DAT ]; then
   echo $DAT exists!
   exit 1
fi

echo "0.95387536e-3" > $DAT

# e,H,g,G
cut -d ' ' -f 1,2,6,7 $PROJ/portbp_del/porbitsdel.res >> $DAT
$PROG < $DAT  > $RES
