echo "0.95387536e-3" >sec1sec2s.dat
cut -d ' ' -f 1,3,4 ../portbp/porbits.res >>sec1sec2s.dat
sec1sec2 <sec1sec2s.dat >sec1sec2s_p3.res
