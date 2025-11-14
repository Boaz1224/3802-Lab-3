clc; clear; close all;

N_size = 50;

cl_test.NACA0006 = zeros(N_size,1);
cl_test.NACA0012 = zeros(N_size,1);
cl_test.NACA0018 = zeros(N_size,1);


%% CL versus Panel #

for N = 3:N_size
    [x_vec_0006, y_vec_0006] = naca4_boundary('0006', N);
    [x_vec_0012, y_vec_0012] = naca4_boundary('0012', N);
    [x_vec_0018, y_vec_0018] = naca4_boundary('0018', N);

%{
figure(1); 
plot(x_vec_0012, y_vec_0012,LineWidth=2); 
axis equal; 
grid on; 
xlabel('x'); 
ylabel('y'); 
title('NACA 2412 boundary (clockwise, LE \rightarrow TE)');
%}

    cl_test.NACA0006(N) = Vortex_Panel(flip(x_vec_0006),flip(y_vec_0006),5);
    cl_test.NACA0012(N) = Vortex_Panel(flip(x_vec_0012),flip(y_vec_0012),5);
    cl_test.NACA0018(N) = Vortex_Panel(flip(x_vec_0018),flip(y_vec_0018),5);

end

figure(2);
hold on;
%plot(1:N_size,cl_test.NACA0006,LineWidth=2);
l1 = plot(3:N_size,cl_test.NACA0012(3:end),LineWidth=2);
%plot(1:N_size,cl_test.NACA0018,LineWidth=2);
hold off;
title("Vortex Panel Method for Increasing Panels")
xlabel("Number of Panels")
ylabel("Coefficient of Lift (c_{l})")
yline(cl_test.NACA0012(end),LineWidth= 1.5,Color="black")
l2 = yline(cl_test.NACA0012(end)+0.01*cl_test.NACA0012(end),LineWidth= 1.5,Color="black",LineStyle="--");
yline(cl_test.NACA0012(end)-0.01*cl_test.NACA0012(end),LineWidth= 1.5,Color="black",LineStyle="--")

legend([l1,l2],{"NACA 0012","Long Term c_{l} +-1%"},Location="southeast")


%% CL versus Angle of Attack
N = 100;
AoA = -10:1:10;

for i = 1:length(AoA)
    [x_vec_0006, y_vec_0006] = naca4_boundary('0006', N);
    [x_vec_0012, y_vec_0012] = naca4_boundary('0012', N);
    [x_vec_0018, y_vec_0018] = naca4_boundary('0018', N);

%{
figure(1); 
plot(x_vec, y_vec,LineWidth=2); 
axis equal; 
grid on; 
xlabel('x'); 
ylabel('y'); 
title('NACA 2412 boundary (clockwise, LE \rightarrow TE)');
%}

    cl_test.NACA0006(i) = Vortex_Panel(flip(x_vec_0006),flip(y_vec_0006),AoA(i));
    cl_test.NACA0012(i) = Vortex_Panel(flip(x_vec_0012),flip(y_vec_0012),AoA(i));
    cl_test.NACA0018(i) = Vortex_Panel(flip(x_vec_0018),flip(y_vec_0018),AoA(i));

end

figure(3);
hold on;
plot(AoA,cl_test.NACA0006(1:length(AoA)),LineWidth=2);
plot(AoA,cl_test.NACA0012(1:length(AoA)),LineWidth=2);
plot(AoA,cl_test.NACA0018(1:length(AoA)),LineWidth=2);
hold off;
title("Vortex Panel Method Multiple Airfoils")
xlabel("Angle of Attack (\alpha)")
ylabel("Coefficient of Lift (c_{l})")
legend("NACA 0006","NACA 0012","NACA 0018",Location="southeast")

n06 = polyfit(AoA,cl_test.NACA0006(1:length(AoA)),1)
n12 = polyfit(AoA,cl_test.NACA0012(1:length(AoA)),1)
n18 = polyfit(AoA,cl_test.NACA0018(1:length(AoA)),1)
