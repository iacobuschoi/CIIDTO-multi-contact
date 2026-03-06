addpath('../contact');
addpath('../InverseDynamics/');

p = getParams();
m1 = p.m_vec(1);
m2 = p.m_vec(2);
m3 = p.m_vec(3);
r1 = p.r_vec(1);
r2 = p.r_vec(2);
r3 = p.r_vec(3);

M  = diag([m1 m1 m2 m2 m3 m3]);

gc = [0; 0;
      0.3; 0.5;
      0.81; 1];
gv = [0; 0; 0; 0; 0; 0];

dt = p.dt_MPC;
T  = p.horizon;
N  = p.N;
theta = linspace(0, 2*pi, 80);

figure; clf;
axis equal; grid on; hold on;

B = zeros(6,2);
B(1:2,1:2) = eye(2);
q_log = [];

for i = 1:N
    q_log = [q_log;gc];
    if i==1
        u = [0.4;1.2]*m1/dt;
    else
        u = [0;0];
    end
    gf = getContactForce(gc, gv);
    ga = M \ (gf + B*u);

    gv = gv + ga * dt;
    gc = gc + gv * dt;

    %%%% plot
    x1 = gc(1); y1 = gc(2);
    x2 = gc(3); y2 = gc(4);
    x3 = gc(5); y3 = gc(6);

    c1x = x1 + r1*cos(theta);  c1y = y1 + r1*sin(theta);
    c2x = x2 + r2*cos(theta);  c2y = y2 + r2*sin(theta);
    c3x = x3 + r3*cos(theta);  c3y = y3 + r3*sin(theta);

    cla;
    plot(c1x, c1y, 'LineWidth', 2); hold on;
    plot(c2x, c2y, 'LineWidth', 2);
    plot(c3x, c3y, 'LineWidth', 2);

    plot([x1 x2 x3], [y1 y2 y3], '.', 'MarkerSize', 12);

    axis equal; grid on;
    xlim([-1 2]);
    ylim([-1 2]);

    title(sprintf('t = %.3f s', i*dt));
    drawnow;
    pause(dt);
end