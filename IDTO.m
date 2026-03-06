addpath('contact');
addpath('InverseDynamics');

close all; clc; clear;

p = getParams();
dt = p.dt_MPC;
nq = p.nq;
N  = p.N;

q0 = [0; 0;
      0.3; 0.5;
      0.8; 1];
des= [1;1];
v0 = zeros(6,1);

qd = repmat([zeros(4,1);des], N, 1);
vd = zeros(nq*N, 1);

q  = genPath();

% [l,q,h] = getCostDerivatives(q,qd,vd,q0,v0,p);
%%
% Solver skeleton: damped Newton + backtracking Armijo line search
% Assumes you have:
%   [l,g,h] = getCostDerivatives(q,qd,vd,q0,v0,p);
%   l = getCost(q,qd,vd,q0,v0,p);   % or reuse l from derivatives
%
% Key upgrades vs your snippet:
% 1) Regularize Hessian to ensure descent direction
% 2) Guard for non-finite cost
% 3) Termination criteria
% 4) Step acceptance + lambda update (LM-style)
% 5) Optional clamp on step length

maxIter   = 50;
tolGrad   = 1e-6;
tolStep   = 1e-10;

rho_ls    = 0.6;     % backtracking shrink
c1        = 1e-1;    % Armijo
alpha_min = 1e-10;

lambda    = 1e-6;    % Hessian regularization
lambdaUp  = 10;
lambdaDn  = 0.3;

cost = @(qq) getCost(qq,qd,vd,q0,v0,p);

for it = 1:maxIter
    disp(q(1+(N-1)*nq:N*nq).');
    [l,g,H] = getCostDerivatives(q,qd,vd,q0,v0,p);
    ng = norm(g);

    fprintf('it=%02d  cost=%.6e  ||g||=%.3e  lambda=%.3e\n', it, l, ng, lambda);

    if ~isfinite(l) || any(~isfinite(g)) || any(~isfinite(H(:)))
        error('Non-finite cost/derivatives');
    end

    if ng < tolGrad
        break;
    end

    accepted = false;
    for trial = 1:12
        Hreg = H + lambda*eye(size(H));

        pp = -Hreg \ g;

        if or(~isfinite(pp),(g.'*pp >= 0))
            lambda = lambda * lambdaUp;
            continue;
        end

        alpha = 1.0;
        while alpha > alpha_min
            q_new = q + alpha*pp;
            l_new = cost(q_new);

            if and(isfinite(l_new), (l_new <= l + c1*alpha*(g.'*pp)))
                q = q_new;
                lambda = max(1e-12, lambda*lambdaDn);
                accepted = true;
                break;
            end

            alpha = alpha * rho_ls;
        end

        if accepted
            break;
        else
            lambda = lambda * lambdaUp;
        end
    end

    if ~accepted
        warning('Step rejected (lambda too large / alpha too small). Stopping.');
        break;
    end

    if norm(alpha*pp) < tolStep
        break;
    end
end

disp('q_N =');
disp(q(1+(N-1)*nq:N*nq).');
%%

function [l,g,h] = getCostDerivatives(q,qd,vd,q0,v0,p)
    %%%% Params
    dt = p.dt_MPC;
    nq = p.nq;
    N  = p.N;
    R  = p.R;
    R0 = p.R0;
    Qq = p.Qq;
    Qv = p.Qv;
    Qqf= p.Qqf;
    Qvf= p.Qvf;
    v  = zeros(nq*N,1);
    u0 = 0;
    u  = zeros(nq*(N-1),1);
    
    %%%% Calculat Velocity
    v(1:nq) = (q(1:nq) - q0)/dt;
    for i = 2:N
        v(1+(i-1)*nq:i*nq) = (q(1+(i-1)*nq:i*nq) - q(1+(i-2)*nq:(i-1)*nq))/dt;
    end
    
    %%%% Inverse Dynamics
    u0 = getU(q0 - v0*dt, q0, q(1:nq));
    for i = 1:N-1
        if i==1
            q1 = q0;
            q2 = q(1:nq);
            q3 = q(nq+1:2*nq);
        else
            q1 = q(1+(i-2)*nq:(i-1)*nq);
            q2 = q(1+(i-1)*nq:i*nq);
            q3 = q(1+i*nq:(i+1)*nq);
        end
        u(1+(i-1)*nq:i*nq) = getU(q1,q2,q3);
    end

    %%%% Cost
    l = 0.5*u0'*R0*u0;
    for i=1:N
        dq = q(1+(i-1)*nq:i*nq) - qd(1+(i-1)*nq:i*nq);
        dv = v(1+(i-1)*nq:i*nq) - vd(1+(i-1)*nq:i*nq);
        if i==N
            l = l + 0.5*(dq'*Qqf*dq + dv'*Qvf*dv);
        else
            du = u(1+(i-1)*nq:i*nq);
            l = l + 0.5*(dq'*Qq*dq + dv'*Qv*dv + du'*R*du);
        end
    end
    
    % cost gredient
    g = zeros(nq*N,1);
    for k=1:N
        g(1+(k-1)*nq:k*nq) = ldq(k-1,k) + ldq(k,k) + ldq(k+1,k);
    end

    % cost hessian
    h = zeros(nq*N,nq*N);
    for k=3:N
        h(1+(k-3)*nq:(k-2)*nq,1+(k-1)*nq:k*nq) = ldqdq(k,k-2);
        h(1+(k-1)*nq:k*nq,1+(k-3)*nq:(k-2)*nq) = ldqdq(k-2,k);
    end
    for k=2:N
        h(1+(k-2)*nq:(k-1)*nq,1+(k-1)*nq:k*nq) = ldqdq(k,k-1);
        h(1+(k-1)*nq:k*nq,1+(k-2)*nq:(k-1)*nq) = ldqdq(k-1,k);
    end
    for k=1:N
        h(1+(k-1)*nq:k*nq,1+(k-1)*nq:k*nq) = ldqdq(k,k);
    end
    


    %%%% Functions
    % Cost gradient
    function x = ldq(i,k)
        %%%% dl_i/dq_j output: (6,1)
        if i==k-1
            if i==0
                du = u0;
                x = du'*R0*udq(k-1,k);
            else
                du = u(1+(i-1)*nq:i*nq);
                x = du'*R*udq(k-1,k);
            end
        elseif i==k
            if k < N
                dq = q(1+(i-1)*nq:i*nq) - qd(1+(i-1)*nq:i*nq);
                dv = v(1+(i-1)*nq:i*nq) - vd(1+(i-1)*nq:i*nq);
                du = u(1+(i-1)*nq:i*nq);
                x  = dq'*Qq + dv'*Qv*vdq(i,k) + du'*R*udq(i,k);
            else
                dq = q(1+(i-1)*nq:i*nq) - qd(1+(i-1)*nq:i*nq);
                dv = v(1+(i-1)*nq:i*nq) - vd(1+(i-1)*nq:i*nq);
                x  = dq'*Qqf + dv'*Qvf*vdq(i,k);
            end
        else
            if k < N-1
                dv = v(1+(i-1)*nq:i*nq) - vd(1+(i-1)*nq:i*nq);
                du = u(1+(i-1)*nq:i*nq);
                x  = dv'*Qv*vdq(i,k) + du'*R*udq(i,k);
            elseif k == N-1
                dv = v(1+(i-1)*nq:i*nq) - vd(1+(i-1)*nq:i*nq);
                x  = dv'*Qvf*vdq(i,k);
            else
                x = zeros(1,6);
            end
        end
        x = x';
    end
    
    % Cost Hessian
    function x = ldqdq(k,j)
        %%%% H_k,j
        if j == k-2     % k = 3~N
            x = udq(k-1,k-2)'*R*udq(k-1,k);
        elseif j == k-1 % k = 2~N
            if k < N 
                x = udq(k-1,k-1)' * R  * udq(k-1,k) +...
                    vdq(k,k-1)'   * Qv * vdq(k,k) +...
                    udq(k,k-1)'   * R  * udq(k,k);
            else
                x = udq(k-1,k-1)' * R  * udq(k-1,k) +...
                    vdq(k,k-1)'   * Qvf* vdq(k,k);
            end
        elseif j == k   % k = 1~N
            if k == 1
                x = udq(k-1,k)' * R0 * udq(k-1,k) +...
                                  Qq              +...
                    vdq(k,k)'   * Qv * vdq(k,k)   +...
                    udq(k,k)'   * R  * udq(k,k)   +...
                    vdq(k+1,k)' * Qv * vdq(k+1,k) +...
                    udq(k+1,k)' * R  * udq(k+1,k);
            elseif k < N-1
                x = udq(k-1,k)' * R  * udq(k-1,k) +...
                                  Qq              +...
                    vdq(k,k)'   * Qv * vdq(k,k)   +...
                    udq(k,k)'   * R  * udq(k,k)   +...
                    vdq(k+1,k)' * Qv * vdq(k+1,k) +...
                    udq(k+1,k)' * R  * udq(k+1,k);
            elseif k == N-1
                x = udq(k-1,k)' * R  * udq(k-1,k) +...
                                  Qq              +...
                    vdq(k,k)'   * Qv * vdq(k,k)   +...
                    udq(k,k)'   * R  * udq(k,k)   +...
                    vdq(k+1,k)' * Qvf* vdq(k+1,k);
            else % k == N
                x = udq(k-1,k)' * R  * udq(k-1,k) +...
                                  Qqf             +...
                    vdq(k,k)'   * Qvf* vdq(k,k);
            end
        else % j > k
            x = ldqdq(j,k)';
        end
    end

    % Inverse-Dynamics Derivatives
    function x = vdq(i,j)
        %%%% dv_i/dq_j
        if i==j
            x = getVdqn();
        else
            x = getVdqm();
        end
    end

    function x = udq(i,j)
        %%%% du_i/dq_j
        if i==0
            q1 = q0 - v0*dt;
            q2 = q0;
            q3 = q(1:nq);
        elseif i==1
            q1 = q0;
            q2 = q(1:nq);
            q3 = q(nq+1:2*nq);
        else
            q1 = q(1:nq);
            q2 = q(nq+1:2*nq);
            q3 = q(2*nq+1:3*nq);
        end
        if j==i-1
            x = getUdqm(q1,q2,q3);
        elseif j==i
            x = getUdqn(q1,q2,q3);
        else
            x = getUdqp(q1,q2,q3);
        end
    end
end

function l = getCost(q,qd,vd,q0,v0,p)
    %%%% Params
    dt = p.dt_MPC;
    nq = p.nq;
    N  = p.N;
    R  = p.R;
    R0 = p.R0;
    Qq = p.Qq;
    Qv = p.Qv;
    Qqf= p.Qqf;
    Qvf= p.Qvf;
    v  = zeros(nq*N,1);
    u0 = 0;
    u  = zeros(nq*(N-1),1);
    
    %%%% Calculat Velocity
    v(1:nq) = (q(1:nq) - q0)/dt;
    for i = 2:N
        v(1+(i-1)*nq:i*nq) = (q(1+(i-1)*nq:i*nq) - q(1+(i-2)*nq:(i-1)*nq))/dt;
    end
    
    %%%% Inverse Dynamics
    u0 = getU(q0 - v0*dt, q0, q(1:nq));
    for i = 1:N-1
        if i==1
            q1 = q0;
            q2 = q(1:nq);
            q3 = q(nq+1:2*nq);
        else
            q1 = q(1+(i-2)*nq:(i-1)*nq);
            q2 = q(1+(i-1)*nq:i*nq);
            q3 = q(1+i*nq:(i+1)*nq);
        end
        u(1+(i-1)*nq:i*nq) = getU(q1,q2,q3);
    end

    %%%% Cost
    l = 0.5*u0'*R0*u0;
    for i=1:N
        dq = q(1+(i-1)*nq:i*nq) - qd(1+(i-1)*nq:i*nq);
        dv = v(1+(i-1)*nq:i*nq) - vd(1+(i-1)*nq:i*nq);
        if i==N
            l = l + 0.5*(dq'*Qqf*dq + dv'*Qvf*dv);
        else
            du = u(1+(i-1)*nq:i*nq);
            l = l + 0.5*(dq'*Qq*dq + dv'*Qv*dv + du'*R*du);
        end
    end
end

function [q_sol, f_val, info] = solve_ipopt(x,q_init, Qd, Vd, p)
    %%%% ================= Problem size =================
    nq = p.nq;
    N  = p.N;
    n  = nq * N;
    
    %%%% ================= Bounds =================
    lbx = -inf(n,1);
    ubx =  inf(n,1);
    
    %%%% ================= Constraint size =================
    ceq0 = getEqualityConstraint(x, q_init);
    nc   = length(ceq0);
    
    cl = zeros(nc,1);
    cu = zeros(nc,1);
    
    %%%% ================= IPOPT problem =================
    problem.objective   = @objective;
    problem.gradient    = @gradient;
    problem.constraints = @constraints;
    problem.jacobian    = @jacobian;
    problem.hessian     = @hessian;
    
    problem.jacobianstructure = @jacobianstructure;
    problem.hessianstructure  = @hessianstructure;
    
    problem.x0  = q_init;
    problem.x_L = lbx;
    problem.x_U = ubx;
    problem.cl  = cl;
    problem.cu  = cu;
    
    %%%% ================= IPOPT options =================
    options.ipopt.max_iter = 1000;
    options.ipopt.tol = 1e-6;
    options.ipopt.acceptable_tol = 1e-4;
    options.ipopt.print_level = 0;
    options.ipopt.sb = 'yes';
    options.ipopt.mu_strategy = 'adaptive';
    
    %%%% ================= Solve =================
    [q_sol, info] = ipopt(q_init, problem, options);
    f_val = getCost(q_sol, x, Qd, Vd, p);
    
    %%%% =================================================
    %%%% =============== CALLBACKS =======================
    %%%% =================================================
    function f = objective(Q)
        f = getCost(Q, x, Qd, Vd, p);
    end
    function g = gradient(Q)
        g = getDerivativeCost(Q, x, Qd, Vd, p);
    end
    function c = constraints(Q)
        [c,~] = getEqualityConstraint(x, Q);
    end
    function J = jacobian(Q)
        [~,J] = getEqualityConstraint(x, Q);
        J = sparse(double(J));
    end
    function H = hessian(Q, sigma, lambda)
        H = sigma * getHessianCost(Q, x, p);
        H = sparse(tril(real(double(H))));
    end


    %%%% =============== SPARSITY =====================
    function S = hessianstructure()
        S = sparse(n,n);

        % 3-block banded symmetric Hessian
        for i = 0:2
            for j = 0:N-i-1
                idx1 = j*nq + (1:nq);
                idx2 = (j+i)*nq + (1:nq);

                S(idx2, idx1) = 1;
                S(idx1, idx2) = 1;   % symmetry
            end
        end
        S = tril(S);
        S = sparse(S);
    end
    function S = jacobianstructure()
        % Safe default (can be tightened)
        S = sparse(ones(nc, n));
    end
end
