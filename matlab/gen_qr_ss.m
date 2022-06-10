%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script defines the state-space representation of a quadrotor       %
% linearized about the (0,0,0,...0,0) equilibrium, analyzes the system    %
% for controllability and observability, and finally designs an LQR       %
% state-feedback controller and observer to stabilize the system.         %
% Author: Zach Bortoff                                                    %
% Last Updated: June 10, 2022                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define symbolic system parameters, config variables, and state variables
% system params
syms M m L l g real

% state variables
syms x y z alpha beta gamma xdot ydot zdot alphadot betadot gammadot real

% state variable derivatives
syms xddot yddot zddot alphaddot betaddot gammaddot real

% equilibrium angular velocity
u0 = sqrt(1 / 4 * (M + 4 * m) * g);

% X = [x y z alpha beta gamma xdot ydot zdot alphadot betadot gammadot]';
A = [
     0 0 0 0     0    0     1    0    0    0        0       0      ;
     0 0 0 0     0    0     0    1    0    0        0       0      ;
     0 0 0 0     0    0     0    0    1    0        0       0      ;
     0 0 0 0     0    0     0    0    0    1        0       0      ;
     0 0 0 0     0    0     0    0    0    0        1       0      ;
     0 0 0 0     0    0     0    0    0    0        0       1      ;
     0 0 0 0     1/(M+4*m)*4*u0^2    0    0    0    0    0        0       0      ;
     0 0 0 0     0    1/(M+4*m)*4*u0^2     0    0    0    0        0       0      ;
     0 0 0 0     0    0     0    0    0    0        0       0      ;
     0 0 0 0     0    0     0    0    0    0        0       0      ;
     0 0 0 0     0    0     0    0    0    0        0       0      ;
     0 0 0 0     0    0     0    0    0    0        0       0      
];
B = [
     0  0  0  0 ;
     0  0  0  0 ;
     0  0  0  0 ;
     0  0  0  0 ;
     0  0  0  0 ;
     0  0  0  0 ;
     0  0  0  0 ;
     0  0  0  0 ;
     1/(M+4*m)*2*u0 1/(M+4*m)*-2*u0 1/(M+4*m)*2*u0 1/(M+4*m)*-2*u0
     m*l/((M+4*m)*L) m*l/((M+4*m)*L) m*l/((M+4*m)*L) m*l/((M+4*m)*L)
     1/m*2*u0 0 -1/m*2*u0 0
     0 1/m*2*u0 0 -1/m*2*u0
];

% only system outputs are (x,y,z,alpha); seems weird, but whatever
C = [
    1 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
];

D = zeros(6,4);

%% Test for controllability and observability
Mc = simplify([B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B A*A*A*A*A*A*B A*A*A*A*A*A*A*B A*A*A*A*A*A*A*A*B A*A*A*A*A*A*A*A*A*B A*A*A*A*A*A*A*A*A*A*B A*A*A*A*A*A*A*A*A*A*A*B]);
rank(double(simplify(subs(Mc, [M m L l g], [1 1 0.2 0.15 9.81])))) % rank is 12 -> therefore controllable

Mob = simplify([C; C*A; C*A*A; C*A*A*A; C*A*A*A*A; C*A*A*A*A*A; C*A*A*A*A*A*A; C*A*A*A*A*A*A*A; C*A*A*A*A*A*A*A*A; C*A*A*A*A*A*A*A*A*A; C*A*A*A*A*A*A*A*A*A*A; C*A*A*A*A*A*A*A*A*A*A*A; C*A*A*A*A*A*A*A*A*A*A*A*A]);
rank(double(simplify(subs(Mob, [M m L l g], [1 1 0.2 0.15 9.81])))) % rank is 12 -> therefore observable

A = double(subs(A, [M m L l g], [1 1 0.2 0.15 9.81]));
B = double(subs(B, [M m L l g], [1 1 0.2 0.15 9.81]));
C = double(subs(C, [M m L l g], [1 1 0.2 0.15 9.81]));

%% Controller Design with "place(...)" Command
K = place(A,B,[-1 -1.1 -1.2 -1.3 -1.4 -1.5 -1.6 -1.7 -1.8 -1.9 -1.91 -1.92]);

%% Controller Design with "lqr(...)" Command
Q = eye(12, 12); R = eye(4, 4);
Q(1,1) = 100;
Q(2,2) = 100;
Q(3,3) = 100;
Q(4,4) = 10;
Q(5,5) = 1e3;
Q(6,6) = 1e3;
[K,~,system_poles] = lqr(A, B, Q, R);  % internal dynamics pole placement
Mo = [Q; Q*A; Q*A*A; Q*A*A*A; Q*A*A*A*A; Q*A*A*A*A*A; Q*A*A*A*A*A*A; Q*A*A*A*A*A*A*A; Q*A*A*A*A*A*A*A*A; Q*A*A*A*A*A*A*A*A*A; Q*A*A*A*A*A*A*A*A*A*A; Q*A*A*A*A*A*A*A*A*A*A*A; Q*A*A*A*A*A*A*A*A*A*A*A*A]; % (A,Q) is observable => detectible
[H,~,observer_poles] = lqr(A', C', eye(12), eye(4,4)); % observer error dynamics pole placement

%% Augment the system with integral action
Atilde = [
    A zeros(12, 4);
    C zeros(4, 4)
];
Btilde = [
    B;
    zeros(4, 4)
];
rank(ctrb(Atilde, Btilde)) % rank = 16 => system is fully controllable with augmentation

% design a controller for the integral action state-feedback
[K, ~, poles] = lqr(Atilde, Btilde, eye(16), eye(4));
K1 = K(:, 1:12)
K2 = K(:, 13:16)

% design an observer by placing the poles at least 2x further into OLHP
H = place(A', C', [-0.8, -0.81, -0.13, -1.8, -1.81, -1.82, -1.83, -1.84, -1.85, -3, 3.01, -3.02])'

%% Perform basic simulation of linear system using 'lsim(...)' command
Acl = [A-B*K1 -B*K2 B*K1; C zeros(4,4) zeros(4,12); zeros(12,12) zeros(12,4) A-H*C];
Bcl = [zeros(12,4); -eye(4,4); zeros(12,4)];
Ccl = [C zeros(4,4) zeros(4, 12)];
Dcl = zeros(4,4);

cls_sf = ss(Acl, Bcl, Ccl, Dcl);
T = 0:0.001:15;

r_step = ones(1, length(T))';
r_zero = zeros(1, length(T))';
r_ramp = T' .* r_step;
yA = lsim(cls_sf, [r_ramp r_zero r_zero r_zero], T, zeros(28,1))
yB = lsim(cls_sf, [r_zero r_step r_zero r_zero], T, zeros(28,1))
yC = lsim(cls_sf, [r_zero r_zero r_step r_zero], T, zeros(28,1))
yD = lsim(cls_sf, [r_zero r_zero r_zero r_step], T, zeros(28,1))
figure;
subplot(2,2,1)
plot(T, yA, 'k'); hold on; plot(T, r_ramp, 'k--'); hold off;
subplot(2,2,2)
plot(T, yB, 'k'); hold on; plot(T, r_step, 'k--'); hold off;
subplot(2,2,3)
plot(T, yC, 'k'); hold on; plot(T, r_step, 'k--'); hold off;
subplot(2,2,4)
plot(T, yD, 'k'); hold on; plot(T, r_step, 'k--'); hold off;

% it's linear, so it doesn't matter, but perturb multiple things at once
yE = lsim(cls_sf, [r_step r_step r_step r_zero], T, zeros(28,1))
plot(T, yE)