echo "0.95387536e-3" >sec1sec2_2d.dat
cut -d ' ' -f 1,3,5 ../portbp/porbits.res >>sec1sec2_2d.dat
sec1sec2_2d <sec1sec2_2d.dat >sec1sec2_2d.res
rm sec1sec2_2d.dat
