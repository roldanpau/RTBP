# Join intersection points computed for +h_1, -h_1, +h_2, -h_2 in one single
# file.

# H, G value corresponding to g=0+h
cut -d ' ' -f 1,8 intersecs_unst_SEC2_br2_h1p.res > temp1.txt   
# G value corresponding to g=0-h
cut -d ' ' -f 8 intersecs_unst_SEC2_br2_h1m.res > temp2.txt   
# G value corresponding to g=0+2h
cut -d ' ' -f 8 intersecs_unst_SEC2_br2_h2p.res > temp3.txt   
# G value corresponding to g=0-2h
cut -d ' ' -f 8 intersecs_unst_SEC2_br2_h2m.res > temp4.txt   

paste -d ' ' temp1.txt temp2.txt temp3.txt temp4.txt |tac > \
joint_results_SEC2_br2.res

# Compute the splitting angle using an AWK script.
awk -f splitting.awk joint_results_SEC2_br2.res > splitting_SEC2_br2.res
