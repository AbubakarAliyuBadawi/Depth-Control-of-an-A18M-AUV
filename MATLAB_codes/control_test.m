clear; clc; close all

%% Inputs
u0 = 2;
dt = 0.1;

Q = 1 * diag([100, 1e-6, 1e-5, 3, 1e-5]);
%Q = 1 * diag([0,0,5e5,0,1e-1]);
R = 1e3;

%% Computations

addpath('Tools/JsonLab')
AUVparams = loadjson('Conf/AUVParameters.json');

rho = 1026;
%nu = AUVparams.Environment.nu;
Sref = 0.385;
Lref = 5;

CMuw = 5;
CMuq = -6;

CZuw = -3;
CZuq = -4;

CZ0 = -0.01;
CM0 = -0.002;

m = 1974.26;
Zw_dot = 1845.9;
Zq_dot = 0;
Mv_dot = 0;
Iy = 4173.5;
Mq_dot = 3388.7;
zg = 3e-2;
zb = 0;

CZ = -0.2;
CM = -0.3;

gravity = 9.81;
%V = AUVparams.Mecanic.volume;

W = m * gravity;
B = W;

k = 0.5*rho*Sref;

%% Etats Equilibres

syms w q z int_z theta BARx

A_cal = u0 * [0; m*w];
B_cal = [0, m*(zg*q + u0); -m*(zg*q + u0), 0] * [w; q];
C_cal = k * [CZ0*abs(u0); Lref * CM0 * abs(u0)] * u0;
D_cal = k * [CZuw * u0, CZuq*Lref*u0; CMuw*Lref*u0, Lref^2*CMuq*u0]*[w; q];
E_cal = [(W - B)*cos(theta); -(zg*W - zb*B)*sin(theta)];
F_cal = k * [CZ*BARx*u0* abs(u0); Lref*CM*BARx*u0*abs(u0)];
M_cal = [(m+Zw_dot), Zq_dot; Mv_dot, (Iy + Mq_dot)];

z_dot = -u0*sin(theta) + w*cos(theta);

% The value of X_dot
f = [M_cal \ (A_cal + B_cal + C_cal + D_cal + E_cal + F_cal);
    z_dot; q; z];

% Computing Equilibrium Point
eqn = f(1:4) == 0;

%Solving the equation for the equilibrium point
S = vpasolve(eqn, [w, q, theta, BARx]);

% Equilibrium Points
w_e = double(S.w)
q_e = double(S.q)
theta_e = double(S.theta)
BAR_e = double(S.BARx)

%% A, B and C matrices

% Computing the Jacobians
A_mat = Jacobian(f, [w; q; z; theta; int_z]);
B_mat = Jacobian(f, BARx);

digits(4)

A = vpa(subs(A_mat, [theta; q; w; BARx], [theta_e; q_e; w_e; BAR_e]))
B = vpa(subs(B_mat, [theta; q; w; BARx], [theta_e; q_e; w_e; BAR_e]))

C = eye(5);

%% LQR Command

% Check for stability
%eig()

Qy = C'*Q*C;

K = lqrd(A,B,Qy,R,dt);
