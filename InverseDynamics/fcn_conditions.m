function out = fcn_conditions(vr)

cond_1 = vr<0;
true_1 = 1-vr;
cond_2 = vr<2;
true_2 = ((vr-2)^2)/4;
true_3 = 0;
out = if_else(cond_1,true_1,if_else(cond_2,true_2,true_3));