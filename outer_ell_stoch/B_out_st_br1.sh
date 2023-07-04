# Compute function $B_{out}^1$ associated to outer map of elliptic problem.
#
# Actually, we only compute the integral called $B^+$ in my notes. This is
# enough, because 
#    B_{out}^j = -\mu B^+ + \mu C^+$,
# where
#    re(C^+)=re(B^+), and im(C^+)=-im(B^+).

NAMEROOT=B_out_st_br1

DATFILE=$NAMEROOT.dat
RESFILE=$NAMEROOT.res
ERRFILE=$NAMEROOT.err

# Select only the first 117 lines, corresponding to energies H<-1.4854.
# The reason for this is because, for larger energies, g=0/PI is not a surface
# of section for the flow (either for po or hom.).

#lines=`wc -l ../outer_circ/omega_neg_SECg_br1.res | cut -d ' ' -f 1`
lines=117

cut -d ' ' -f 1 ../portbp/porbits.res >temp1 	# H
cut -d ' ' -f 1-4 ../prtbp_del/periodicpoints_SECg2.res \
    > temp2  # periodic point p (in Delaunay coords, with g=\pi).

# zs: preimage of homoclinic point z (in Delaunay coords, with g=\pi)
cut -d ' ' -f 1-4 ../prtbp_del/intersecs_preimgs_st_br1_SECg2.res > temp3

# z_s: preimage of homoclinic point z
# We could use 
#   cut -d ' ' -f 13-16 ../intersec_del_car/intersecs_st_SECg_br1.res > temp3
# but we decide to simply use the symmetry z_s = z_u, but reversing the sign of
# l and g. (This is done in the program).

# t: integration time to reach homoclinic point z from z_s
cut -d ' ' -f 6 ../intersec/intersecs_st_br1.res > temp4

echo "0.95387536e-3 1" > $DATFILE	# mu, STABLE

# H, p, z_u, \omega_pos, N
paste -d ' ' temp1 temp2 temp3 ../outer_circ/omega_pos_st_br1.res \
	temp4 | sed -n '1,117p' >> $DATFILE
	
rm temp*

nohup ./outer_ell_stoch <$DATFILE >$RESFILE 2>$ERRFILE &
