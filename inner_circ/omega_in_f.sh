echo "0.95387536e-3 1" >omega_in_f.dat
cut -d ' ' -f 1-4 ../prtbp_del/prtbpdel_p2.res >>omega_in_f.dat
inner_circ <omega_in_f.dat >omega_in_f.res
rm omega_in_f.dat

# Compute partial shift $\omega_in^f$ corresponding to INNER separatrix
# (bInner flag == 1). This will be used later for outer map.
# In particular, compute integral along the periodic orbit with initial
# condition p2=(l,L,g,G).
