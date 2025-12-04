clc; clear; close all;

%% Load in NACA data
load("output.mat")
x = out.indp.vec;
y = out.depValues;

poly = polyfit(x,y,4);
x1 = linspace(min(x),max(x));
y1 = polyval(poly,x1);

figure(1);
plot(x,y,"o",LineWidth=2,Color="red")
hold on;
plot(x1,y1,LineWidth=2,Color="black")
hold off;
title("NACA 0012 at Re = {5.9x10^6}")
xlabel("Coefficient of Lift ({c_l})")
ylabel("Coefficient of Drag ({c_d})")
grid on;
legend("Experimental Points","Polyfit",Location="northwest")
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

%% Task 2
% Get Cl across AoA
N = 100;
AoA = -10:0.5:10;

[x_vec_0012, y_vec_0012] = naca4_boundary('0012', N);
[x_vec_2412, y_vec_2412] = naca4_boundary('2412', N);
cl_AoA_0012 = zeros(length(AoA),1);
cl_AoA_2412 = zeros(length(AoA),1);
for i = 1:length(AoA)
    cl_AoA_0012(i) = Vortex_Panel(flip(x_vec_0012),flip(y_vec_0012),AoA(i));
    cl_AoA_2412(i) = Vortex_Panel(flip(x_vec_2412),flip(y_vec_2412),AoA(i));
end

% Fit Cd to Cl_AoA
cd_AoA = polyval(poly, cl_AoA_0012);

%Plot
figure(2);
hold on;
plot(AoA,cd_AoA,LineWidth=2,Color="red")
hold off;
title("Angle of Attack vs Profile Drag Coefficient")
legend("NACA 0012 at Re = {5.9x10^6}",Location="northwest")
xlabel("Angle of Attack (\alpha)")
ylabel("Profile Drag Coefficient ({C_{D,0}})")
grid on;
ax = gca; % axes handle
ax.YAxis.Exponent = 0;



%% Task 3

e = zeros(length(AoA),1);
cdi_AoA = zeros(length(AoA),1);
c_L = zeros(length(AoA),1);

b = 36; % span (ft)
c_r = 5.333; % root chord (ft)
c_t = 3.583; % tip chord (ft)

% Fit CL vs AoA TIP
pt = polyfit(AoA, cl_AoA_0012, 1);
slope_deg_t     = pt(1); % lift slope (per degree)
intercept_t     = pt(2);
% Convert to per radian for PLLT
a0_t        = slope_deg_t * (180/pi);
% Zero-lift AoA
alphaL0_tip_deg = -intercept_t / slope_deg_t;
aero_t          = deg2rad(alphaL0_tip_deg);

% Fit CL vs AoA ROOT
pr = polyfit(AoA, cl_AoA_2412, 1);
slope_deg_r     = pr(1); % lift slope (per degree)
intercept_r     = pr(2);
% Convert to per radian for PLLT
a0_r        = slope_deg_r * (180/pi);
% Zero-lift AoA
alphaL0_root_deg = -intercept_r / slope_deg_r;
aero_r          = deg2rad(alphaL0_root_deg);

for i = 1:length(AoA)
    alpha_tip = AoA(i);
    geo_t = deg2rad(alpha_tip);
    geo_r = deg2rad(alpha_tip + 2);  % +2 deg root twist
    [e(i), c_L(i), cdi_AoA(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
end

cd_tot = cdi_AoA + cd_AoA;

%Plot
figure(3);
hold on;
plot(AoA,cdi_AoA,LineWidth=2,Color="red")
hold off;
title("Angle of Attack vs Induced Drag Coefficient")
legend("NACA 0012 at Re = {5.9x10^6}",Location="northwest")
xlabel("Angle of Attack (\alpha)")
ylabel("Induced Drag Coefficient ({C_{D,i}})")
grid on;
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

figure(4);
hold on;
plot(AoA,cd_tot,LineWidth=2,Color="red")
hold off;
title("Angle of Attack vs Total Drag Coefficient")
legend("NACA 0012 at Re = {5.9x10^6}",Location="northwest")
xlabel("Angle of Attack (\alpha)")
ylabel("Total Drag Coefficient ({C_{D}})")
grid on;
ax = gca; % axes handle
ax.YAxis.Exponent = 0;



%% PLLT Function
function [e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    %halfspan points
    i = (1:N).';
    theta = i*pi/(2*N);
    eta = cos(theta); %nondimensional span
    % Lift slope
    a0 = a0_r + (a0_t - a0_r).*eta;
    %Chord
    c = c_r + (c_t - c_r).*eta;
    %Zero-lift angle
    alphaL0 = aero_r + (aero_t - aero_r).*eta;
    %Geometric AoA(twist + AC AoA)
    alphaG  = geo_r + (geo_t  - geo_r).*eta;

    %linear system for A_(2j-1)
    % Unknowns
    A_mat = zeros(N);       %coefficient matrix
    rhs = alphaG - alphaL0;  %right-hand side

    for ii = 1:N
        for jj = 1:N
            n = 2*jj - 1; %odd index
            sin_nth = sin(n*theta(ii));
            term1 = 4*b/(a0(ii)*c(ii))*sin_nth;
            term2 = n*sin_nth / sin(theta(ii));
            A_mat(ii,jj) = term1 + term2;
        end
    end

    %Solving for A_(2j-1)
    A = A_mat\rhs; %column vector length N
    %Wing lift & induced drag coefficients
    %S and AR
    S = 0.5*(c_r + c_t)*b; %trapezoidal wing area
    AR = b^2 / S;
    % A1 is n=1
    A1 = A(1);

    %Lift coefficient
    c_L = pi*AR*A1;

    %Induced drag coefficient(odd n only)
    c_Di = 0.0;
    for jj = 1:N
        n = 2*jj - 1;
        c_Di = c_Di + n*A(jj)^2;
    end
    c_Di = pi * AR * c_Di;

    %Span efficiency factor
    e = c_L^2 /(pi*AR*c_Di);
end