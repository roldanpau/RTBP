# Put together the four families of $\omega_neg$ in one file

resfile=omega_neg.res

cut -d ' ' -f 1 ../portbp/porbits.res > temp1    # H

cut -d ' ' -f 1 omega_neg_unst_br1.res > temp2
cut -d ' ' -f 1 omega_neg_unst_br2.res > temp3
cut -d ' ' -f 1 omega_pos_st_br1.res > temp4
cut -d ' ' -f 1 omega_pos_st_br2.res > temp5

# Select only the first 117 lines, corresponding to energies H<-1.4854.
# The reason for this is because, for larger energies, the p.o. is not
# in almost resonance anymore, and g=0/PI is not a surface of section.
paste -d ' ' temp1 temp2 temp3 temp4 temp5 | sed -n '1,117p' > $resfile
#rm temp1 temp2 temp3 temp4 temp5
