echo "0.95387536e-3 SEC2 2" >hypers_apo_SEC2.dat
cut -d ' ' -f 1-3 ../portbp_apo/porbits_apo.res >>hypers_apo_SEC2.dat
./hypers <hypers_apo_SEC2.dat >hypers_apo_SEC2.res

rm hypers_apo_SEC2.dat 
