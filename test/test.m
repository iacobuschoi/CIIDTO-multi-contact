close all; clc; clear;

addpath('InverseDynamics');
q1 = randn(6,1);
q2 = randn(6,1);
q3 = randn(6,1);

num_test = 1e3;
test_start = tic;
test_time_log = [];


for ii=1:num_test
    [u_dgc_m, u_dgc_n, u_dgc_p, gv_n_dgc_m, gv_n_dgc_n] = getInverseDynamicsDerivatives(q1,q2,q3); 
end
test_end = toc(test_start);
test_end = test_end/num_test;
test_time_log = [test_time_log test_end];


test_start = tic;
for ii=1:num_test
    % [u_dgc_m, u_dgc_n, u_dgc_p, gv_n_dgc_m, gv_n_dgc_n] = getInverseDynamicsDerivatives(q1,q2,q3); 
    udqm = getUdqm(q1,q2,q3);
end
test_end = toc(test_start);
test_end = test_end/num_test;
test_time_log = [test_time_log test_end];


test_start = tic;
for ii=1:num_test
    udqn = getUdqn(q1,q2,q3);
end
test_end = toc(test_start);
test_end = test_end/num_test;
test_time_log = [test_time_log test_end];


test_start = tic;
for ii=1:num_test
    % [u_dgc_m, u_dgc_n, u_dgc_p, gv_n_dgc_m, gv_n_dgc_n] = getInverseDynamicsDerivatives(q1,q2,q3); 
    udqp = getUdqp(q1,q2,q3);
end
test_end = toc(test_start);
test_end = test_end/num_test;
test_time_log = [test_time_log test_end];


test_start = tic;
for ii=1:num_test
    % [u_dgc_m, u_dgc_n, u_dgc_p, gv_n_dgc_m, gv_n_dgc_n] = getInverseDynamicsDerivatives(q1,q2,q3); 
    vdqm = getVdqm();
end
test_end = toc(test_start);
test_end = test_end/num_test;
test_time_log = [test_time_log test_end];


test_start = tic;
for ii=1:num_test
    % [u_dgc_m, u_dgc_n, u_dgc_p, gv_n_dgc_m, gv_n_dgc_n] = getInverseDynamicsDerivatives(q1,q2,q3); 
    vdqn = getVdqn();
end
test_end = toc(test_start);
test_end = test_end/num_test;
test_time_log = [test_time_log test_end]
