clear; clc; close all

%% Inputs
u0 = 2;
dt = 0.1;

%Q = diag([1,3,5e6,10,1e-5]);
%R = 1e12;

% Q = diag([1e-6,1,100,10,0.1]);
% R = 1e9;


Q = diag([1e8,1e8,350,450,1e-5]);
R = 1e8;
%% Computations

addpath('Tools/JsonLab')
AUVparams = loadjson('Conf/AUVParameters.json');

rho = AUVparams.Environment.Rho;
nu = AUVparams.Environment.nu;
Sref = AUVparams.Mecanic.surface_reference;
Lref = AUVparams.Mecanic.length_reference;

CMuw = AUVparams.Hydrodynamic.CMuw;
CMuq = AUVparams.Hydrodynamic.CMuq;

CZuw = AUVparams.Hydrodynamic.CZuw;
CZuq = AUVparams.Hydrodynamic.CZuq;

CZ0 = AUVparams.Hydrodynamic.CZ0;
CM0 = AUVparams.Hydrodynamic.CM0;

m = AUVparams.Mecanic.mass;
Zw_dot = AUVparams.Mecanic.added_water_matrix(3, 3);
Zq_dot = AUVparams.Mecanic.added_water_matrix(3, 5);
Mw_dot = AUVparams.Mecanic.added_water_matrix(5, 3);
Iy = AUVparams.Mecanic.inertia_matrix(2, 2);
Mq_dot = AUVparams.Mecanic.added_water_matrix(5, 5);
zg = AUVparams.Mecanic.center_of_gravity_veh_ref(3);
zb = AUVparams.Mecanic.center_of_buoyancy_veh_ref(3);

CZ = AUVparams.XRearHelms.CZ;
CM = AUVparams.XRearHelms.CM;

gravity = 9.81;
V = AUVparams.Mecanic.volume;

W = m * gravity;
B = V*rho*gravity;

k = 0.5*rho*Sref;

%% Etats Equilibres

syms w q z int_z theta BAR

A_cal = [0, 0, m*(zg*q + u0); m*w, -m*(zg*q + u0), 0] * [u0; w; q];
B_cal = k * [CZ0*abs(u0), CZuw * u0, CZuq*Lref*u0; Lref * CM0 * abs(u0), CMuw*Lref*u0, Lref^2*CMuq*u0]*[u0; w; q];
C_cal = [(W - B)*cos(theta); -(zg*W - zb*B)*sin(theta)];
D_cal = k * [CZ*BAR*u0*abs(u0); Lref*CM*BAR*u0*abs(u0)];
M_cal = [(m+Zw_dot), Zq_dot; Mw_dot, (Iy + Mq_dot)];

%z_dot = -u0*sin(theta) + w*cos(theta);
z_dot = w;
%z_dot = -u0*theta + w;

% The value of X_dot
f = [M_cal \ (A_cal + B_cal + C_cal + D_cal);
    z_dot; q; z];

% Computing Equilibrium Point
eqn = f(1:4) == 0;

%Solving the equation for the equilibrium point
S = vpasolve(eqn, [w, q, theta, BAR]);

% Equilibrium Points
w_e = double(S.w);
q_e = double(S.q);
theta_e = double(S.theta);
BAR_e = double(S.BAR);

%% A, B and C matrices

% Computing the Jacobians
A_mat = Jacobian(f, [w; q; z; theta; int_z]);
B_mat = Jacobian(f, BAR);

A = double(subs(A_mat, [theta; q; w; BAR], [theta_e; q_e; w_e; BAR_e]))
B = double(subs(B_mat, [theta; q; w; BAR], [theta_e; q_e; w_e; BAR_e]))

C = eye(size(A, 2));

%% LQR Command

% Check for controllability (for discrete)
rank_controllability = rank(ctrb(A, B));

Qy = C'*Q*C;

K = lqrd(A,B,Qy,R,dt)

%% Equilibrium Points

equ_points = [w_e ; q_e; theta_e; BAR_e]