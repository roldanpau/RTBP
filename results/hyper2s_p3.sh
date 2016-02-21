echo "0.95387536e-3 6" >hyper2s_p3.dat
cut -d ' ' -f 1-3 sec1sec2s_p3.res >>hyper2s_p3.dat     # H,p[2]
hyper2s <hyper2s_p3.dat >hyper2s_p3.res
