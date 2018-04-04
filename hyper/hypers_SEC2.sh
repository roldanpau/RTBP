echo "0.95387536e-3 SEC2 2" >hypers.dat
cut -d ' ' -f 1-3 ../sec1sec2/sec1sec2_2d.res >>hypers.dat
hypers <hypers.dat >hypers.res

rm hypers.dat 
