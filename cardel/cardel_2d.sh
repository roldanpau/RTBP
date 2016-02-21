echo "0.95387536e-3 -1.7194" >unstmfld_H-1.7194_p1.dat
cat ../invmfld/unstmfld_H-1.7194_p1.res >>unstmfld_H-1.7194_p1.dat
cardel_2d <unstmfld_H-1.7194_p1.dat >unstmfld_H-1.7194_p1.res

echo "0.95387536e-3 -1.7194" >stmfld_H-1.7194_p1.dat
cat ../invmfld/stmfld_H-1.7194_p1.res >>stmfld_H-1.7194_p1.dat
cardel_2d <stmfld_H-1.7194_p1.dat >stmfld_H-1.7194_p1.res

echo "0.95387536e-3 -1.7194" >stmfld_H-1.7194_p2.dat
cat ../invmfld/unstmfld_H-1.7194_p1.res >temp1.dat
awk 'NF>0 { neg = -$2; print $1, neg }' temp1.dat >>stmfld_H-1.7194_p2.dat
cardel_2d <stmfld_H-1.7194_p2.dat >stmfld_H-1.7194_p2.res

echo "0.95387536e-3 -1.7194" >unstmfld_H-1.7194_p2.dat
cat ../invmfld/stmfld_H-1.7194_p1.res >temp1.dat
awk 'NF>0 { neg = -$2; print $1, neg }' temp1.dat >>unstmfld_H-1.7194_p2.dat
cardel_2d <unstmfld_H-1.7194_p2.dat >unstmfld_H-1.7194_p2.res

rm \
   temp1.dat \
   unstmfld_H-1.7194_p1.dat stmfld_H-1.7194_p1.dat \
   stmfld_H-1.7194_p2.dat unstmfld_H-1.7194_p2.dat
