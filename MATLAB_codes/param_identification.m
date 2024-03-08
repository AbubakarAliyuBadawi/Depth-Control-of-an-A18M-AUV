clear, clc

% Load Data
load('Data/SpeedSteps_Noise_U.mat');
%load('FMINSEARCHBNH/FMINSEARCHBNH/test/fminsearchbnd.m');

%% Load Params
auv_param = jsondecode(fileread('Conf/AUVParameters.json'));

u = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data;

figure;
plot(vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Time,...
    vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data);
   
grid on; ylabel('Speed(m/s)'); xlabel('Time (s)');
title('U');

Fx = forces_values.HydrodynamicForces.Fx_N.Data(:);

%% Params
rho = auv_param.Environment.Rho;
nu = auv_param.Environment.nu;
Lref = auv_param.Mecanic.length_reference;
Sref = auv_param.Mecanic.surface_reference;

%% Finding ksh
Cxf = 0.07./(log10(u*Lref/nu)-2).^2;
eq = 0.5*rho*Sref.*Cxf.*abs(u).*u;
Ksh = lsqr(eq,Fx);

%% Plot Ksh estimation
figure; grid on; hold on;
plot(Fx);
plot(eq*Ksh);
legend('Data', 'Fit');
title('Ksh Fit');

%% 

load('Data/PitchSteps_NoiseOn_Pitch_W_Q.mat');

u = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data; 
t_u = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Time; 

figure; 

plot(t_u,u);
grid on; ylabel('speed(m/s)');xlabel('Time(s)');
title('U'); 


Fx = forces_values.HydrodynamicForces.Fx_N.Data(:); 

%Ax = b 
%x = lsqr(A,b) or A/b

%params
rho = auv_param.Environment.Rho; 
nu = auv_param.Environment.nu ; 
Lref = auv_param.Mecanic.length_reference;
Sref = auv_param.Mecanic.surface_reference;
CZ0 = auv_param.Hydrodynamic.CZ0; 
CM0 = auv_param.Hydrodynamic.CM0;

%acquiring CZuw,CZuq,
Fz = forces_values.HydrodynamicForces.Fz_N.Data(:);
w = vehicle_state.nu.rAUV_WaterSpeed.Wwater_ms.Data; 
q = vehicle_state.nu.AngularSpeed.Q_rads.Data; 

const = 0.5.*rho.*Sref;
A = [u.* w.*const Lref.*u.*q.*const];
b = Fz - CZ0.*abs(u).*u.*const; 
x = lsqr(A,b);
CZuw = x(1)
CZuq = x(2)

%plot CZuw,CZuq, estimation 
figure;grid on; hold on;
plot(Fz);
plot(A*x);
legend('data','fit')
title('CZuw,CZuq fit'); 


%acquiring CMuw,CMuq
My = forces_values.HydrodynamicForces.My_Nm.Data(:);
const = 0.5.*rho.*Sref.*Lref;
Ay = [u.* w.*const Lref.*u.*q.*const];
by = My - CM0.*abs(u).*u.*const; 
y = lsqr(Ay,by);        
CMuw = y(1)
CMuq = y(2)

%plot CMuw,CMuq, estimation 
figure;grid on; hold on;
plot(My);
plot(Ay*y);
legend('data','fit')
title('CMuw,CMuq fit'); 

%%%%
%%
load('Data/DieuDonne_AStep_NoiseOn_U_V_R_A.mat');

u = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data; 
t_u = vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Time; 

figure; 

plot(t_u,u);
grid on; ylabel('speed(m/s)');xlabel('Time(s)');
title('U'); 

%%
Fx = forces_values.HydrodynamicForces.Fx_N.Data(:); 

%Ax = b 
%x = lsqr(A,b) or A/b

%params
rho = auv_param.Environment.Rho; 
nu = auv_param.Environment.nu ; 
Lref = auv_param.Mecanic.length_reference;
Sref = auv_param.Mecanic.surface_reference;

%acquiring CYuv,CYur,
Fy = forces_values.HydrodynamicForces.Fy_N.Data(:);
v = vehicle_state.nu.rAUV_WaterSpeed.Vwater_ms.Data; 
r = vehicle_state.nu.AngularSpeed.R_rads.Data; 

const = 0.5.*rho.*Sref;
Ay= [u.*v.*const Lref.*u.*r.*const];
b = Fy; 
x0CY = [CZuw -CZuq];
funCY = @(xsolCY) sum((Fy -(u.*v.*const.*xsolCY(1) + Lref.*u.*r.*const.*xsolCY(2))).^2);
ub = [0.8*CZuw; -1.025*CZuq];
lb = [1.2*CZuw; -0.975*CZuq];

xsolCY = fminsearchbnd(funCY,x0CY,lb,ub); 

CYuv = xsolCY(1)
CYur = xsolCY(2)

figure;
grid on;
hold on;

plot(Fy, 'DisplayName', 'Reference Fy');
plot(Ay * xsolCY', 'DisplayName', 'Fitted Fy');

legend('show');
title('CYuv, CYur fit');


%acquiring CNuv,CNur
Mz = forces_values.HydrodynamicForces.Mz_Nm.Data(:);

const = 0.5.*rho.*Sref;
An = [u.*v.*Lref.*const Lref.^2.*u.*r.*const];
b = Mz;
x0CN = [-CMuw CMuq];
%ub = [-1.1*CMuw; 0.9*CMuq];
%lb = [-0.9*CMuw; 1.1*CMuq];
lb = [-CMuw-(0.1*CMuw)  CMuq-(0.1*-CMuq)];
ub = [-CMuw+(0.1*CMuw)  CMuq+(0.1*-CMuq)];
funCN = @(xsolCN) sum((Mz -(u.*v.*Lref.*const.*xsolCN(1) + Lref.^2.*u.*r.*const.*xsolCN(2))).^2);
xsolCN = fminsearchbnd(funCN,x0CN,lb,ub); 
CNuv = xsolCN(1)
CNur = xsolCN(2)


% Plotting for yaw moments (CNuv, CNur)
figure;
grid on;
hold on;

plot(Mz, 'DisplayName', 'Reference Mz');
plot(An * xsolCN', 'DisplayName', 'Fitted Mz');

legend('show');
title('CNuv, CNur fit');
%%%