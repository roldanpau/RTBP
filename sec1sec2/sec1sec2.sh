echo "0.95387536e-3" >sec1sec2.dat
cut -d ' ' -f 1,3,5 ../portbp/porbits.res >>sec1sec2.dat
sec1sec2 <sec1sec2.dat >sec1sec2.res
