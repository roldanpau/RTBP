echo "0.95387536e-3" >cardels.dat
cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2.res >>cardels.dat
cardels_2d <cardels.dat >cardels.res
rm cardels.dat
