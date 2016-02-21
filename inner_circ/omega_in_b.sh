echo "0.95387536e-3 0" >omega_in_b.dat
cut -d ' ' -f 1-4 ../prtbp_del/prtbpdel_p1.res >>omega_in_b.dat
inner_circ <omega_in_b.dat >omega_in_b.res
rm omega_in_b.dat

# Compute partial shift $\omega_in^b$ corresponding to OUTER separatrix
# (bInner flag == 0). This will be used later for outer map.
# In particular, compute integral along the periodic orbit with initial
# condition p1=(l,L,g,G).
