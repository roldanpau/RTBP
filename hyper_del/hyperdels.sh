echo "0.95387536e-3 SEC1 3" >hyperdels.dat

# e, H, g, G
cut -d ' ' -f 1,2,6,7 ../portbp_del/porbitsdel.res >>hyperdels.dat
hyperdels <hyperdels.dat >hyperdels.res
