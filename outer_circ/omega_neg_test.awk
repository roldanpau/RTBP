# awk -f omega_neg_test.awk omega_neg_test_1.dat >omega_neg_test_2.dat
NR==1{l=$2;next}{s=$2-l; printf("%d %e\n", 68-NR+1, s); l=$2}
