echo "0.95387536e-3 SEC1 2" >hypers_SEC1.dat
cut -d ' ' -f 1-3 ../portbp/porbits.res >>hypers_SEC1.dat
./hypers <hypers_SEC1.dat >hypers_SEC1.res

rm hypers_SEC1.dat 
