NAMEROOT=periodicpoints_SECg2

datfile=$NAMEROOT.dat
resfile=$NAMEROOT.res
errfile=$NAMEROOT.err

echo "0.95387536e-3 SECg2 1" > $datfile	# mu, sec, num of iterates

# periodic point p (in Delaunay coords).
cut -d ' ' -f 1-4 ../cardel/periodicorbits_del.res >> $datfile

./prtbpdel < $datfile > $resfile
rm $datfile
