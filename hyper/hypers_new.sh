echo "0.95387536e-3 SEC2 2" >hypers_ell00-ell06_new.dat
cut -d ' ' -f 1,3 ../portbp/porbits_ell00-ell06_new.res >>temp1.dat 
awk '{print $1 " " $2 " 0"}' temp1.dat >>hypers_ell00-ell06_new.dat
./hypers <hypers_ell00-ell06_new.dat >hypers_ell00-ell06_new.res

#rm temp1.dat hypers_ell00-ell06_new.dat
