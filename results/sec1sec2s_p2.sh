echo "0.95387536e-3" >sec1sec2s.dat
cut -d ' ' -f 1,3,4 ../portbp/porbits.res >>sec1sec2s.dat
sec1sec2_inv <sec1sec2s.dat >sec1sec2s_p2.res
