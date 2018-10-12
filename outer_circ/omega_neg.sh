# Put together the four families of $\omega_neg$ in one file

resfile=omega_neg.res

cut -d ' ' -f 1 ../portbp/porbits.res > temp1    # H

cut -d ' ' -f 1 omega_neg_SECg_br1.res > temp2
cut -d ' ' -f 1 omega_neg_SECg_br2.res > temp3
cut -d ' ' -f 1 omega_neg_SECg2_br1.res > temp4
cut -d ' ' -f 1 omega_neg_SECg2_br2.res > temp5

paste -d ' ' temp1 temp2 temp3 temp4 temp5 > $resfile
rm temp1 temp2 temp3 temp4 temp5
