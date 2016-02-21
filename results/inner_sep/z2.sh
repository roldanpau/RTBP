# Iterate homoclinic point z2^u until it hits z2.

echo "0.95387536e-3 14" >z2.dat
cut -d ' ' -f 1-4 z2_u.res >>z2.dat	# z2_u
prtbpdel <z2.dat >z2.res
