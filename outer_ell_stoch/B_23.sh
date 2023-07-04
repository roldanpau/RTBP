# Paste result files B_2.res and B_3.res into one single file B_23.res

paste -d ' ' B_2.res B_3.res > B_23.res

# now simply plot the file B_23.res with gnuplot using B_23.plt
