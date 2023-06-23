# Compute integral A_ín(I) related to outer map of elliptic problem.

#lines=`wc -l ../outer_circ/omega_neg_SECg_br1.res | cut -d ' ' -f 1`
lines=117

echo "0.95387536e-3" >inner_ell_stoch.dat	# mu

# Select only the first 117 lines, corresponding to energies H<-1.4854.
# The reason for this is because, for larger energies, g=0/PI is not a surface
# of section for the flow (either for po or hom.).

cut -d ' ' -f 1 ../portbp/porbits.res | head -n $lines > temp1 	# H
cut -d ' ' -f 1-4 ../cardel/periodicorbits_del.res | head -n $lines \
	> temp2  # periodic point p (in Delaunay coords).

paste -d ' ' temp1 temp2 >> inner_ell_stoch.dat
rm temp1 temp2 

nohup inner_ell_stoch < inner_ell_stoch.dat > inner_ell_stoch.res \
    2> inner_ell_stoch.err &
