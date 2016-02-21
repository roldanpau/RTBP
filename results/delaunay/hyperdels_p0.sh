#!/bin/bash
PROG=hyperdels
DAT=hyperdels_p0.dat
RES=hyperdels_p0.res

if [ -f $DAT ]; then
   echo $DAT exists!
   exit 1
fi

echo "0.95387536e-3 SEC2 3" > $DAT
cut -d ' ' -f 1-4 sec1sec2s_p0.res >> $DAT	# e, H, (g,G)
$PROG < $DAT  > $RES
