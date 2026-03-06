%% build_getContactForce.m
close all; clear; clc;
addpath('../');
p = getParams();

mu        = p.mu;
sigma     = p.sigma;
k         = p.k;
v_d       = p.v_d;
v_s       = p.v_s;

mu_ground = p.mu_ground;
v_gs      = p.v_gs;

m1 = p.m_vec(1); m2 = p.m_vec(2); m3 = p.m_vec(3);
g  = p.g;

r1 = p.r_vec(1); r2 = p.r_vec(2); r3 = p.r_vec(3);

gc = sym('gc',[6,1],'real'); % [x1;y1;x2;y2;x3;y3]
gv = sym('gv',[6,1],'real'); % [vx1;vy1;vx2;vy2;vx3;vy3]

% positions
p1 = gc(1:2);
p2 = gc(3:4);
p3 = gc(5:6);

% velocities
v1 = gv(1:2);
v2 = gv(3:4);
v3 = gv(5:6);

% helper: smooth normal force magnitude (same structure as your c())
c_from_phi = @(phi) sigma*k*log(1+exp(-phi/sigma));

% piecewise d(vr) exactly as your code
d_from_vr = @(vr) piecewise( ...
    vr<0, 1-vr, ...
    vr<2, ((vr-2)^2)/4, ...
    0);

% contact pair force (returns world force acting on "i"; force on "j" is -Fi)
pairForce = @(pi, pj, vi, vj, Ri, Rj) localPairForce(pi,pj,vi,vj,Ri,Rj,mu,v_d,v_s,c_from_phi,d_from_vr);

% accumulate ball-ball contact forces (world)
F1 = sym(zeros(2,1)); F2 = sym(zeros(2,1)); F3 = sym(zeros(2,1));

F12 = pairForce(p1,p2,v1,v2,r1,r2); % force on 1 due to 2
F23 = pairForce(p2,p3,v2,v3,r2,r3); % force on 2 due to 3
F13 = pairForce(p1,p3,v1,v3,r1,r3); % force on 1 due to 3

% apply Newton's 3rd law
F1 = F1 + F12 + F13;
F2 = F2 - F12 + F23;
F3 = F3 - F23 - F13;

gf_contact = [F1; F2; F3];

%%%% --------- Ground friction (symbolic, smooth) ----------
% Your original code:
% if norm(v) > eps: f = -v/norm(v) * m*g*mu_ground else 0
% Replace with smooth version:
% f = -(m*g*mu_ground) * v / sqrt(v_gs^2 + ||v||^2)
groundFriction = @(v, m) -(m*g*mu_ground) * v / sqrt(v_gs^2 + v.'*v);

Fg1 = groundFriction(v1, m1);
Fg2 = groundFriction(v2, m2);
Fg3 = groundFriction(v3, m3);

gf_ground = [Fg1; Fg2; Fg3];

%%%% total generalized force output
gf = gf_contact + gf_ground;
% gf = simplify(gf);  % optional; can be slow

matlabFunction( ...
    gf, 'File', 'getContactForce', ...
    'Outputs', {'gf'}, ...
    'Vars', {gc, gv});

disp('Generated getContactForce.m (includes ball-ball contact + ground friction)');

%% -------- local function used by symbolic construction ----------
function Fi = localPairForce(pi,pj,vi,vj,Ri,Rj,mu,v_d,v_s,c_from_phi,d_from_vr)
    % gap phi = distance - (Ri+Rj)
    dp = pi - pj;
    dist = sqrt(dp.'*dp + 1e-12); % avoid singularity
    n = dp / dist;                % normal (from j->i)
    t = [-n(2); n(1)];            % tangent

    vrel = vi - vj;
    vt = t.'*vrel;
    vn = n.'*vrel;

    phi = dist - (Ri + Rj);

    c  = c_from_phi(phi);
    vr = vn / v_d;
    d  = d_from_vr(vr);

    % contact-frame forces (same structure as your model)
    ft = -mu * vt * c * d / sqrt(v_s^2 + vt^2);
    fn =  c * d;

    % world force on i
    Fi = t*ft + n*fn;
end