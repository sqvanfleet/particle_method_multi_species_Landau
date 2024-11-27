clear all 
close all 
clc


data_dt_02 = ...
    load('multi_species_particle_2d_Coulomb_n_50_Example_9_gamma_-3_dv1_0.16_dv2_0.16_dt_0.02tmax50.mat');
data_dt_04 = ...
    load('multi_species_particle_2d_Coulomb_n_50_Example_9_gamma_-3_dv1_0.16_dv2_0.16_dt_0.04tmax50.mat');
data_dt_005 = ...
    load('multi_species_particle_2d_Coulomb_n_50_Example_9_gamma_-3_dv1_0.16_dv2_0.16_dt_0.005tmax50.mat');


Nt_02 = data_dt_02.Nt;
Nt_04 = data_dt_04.Nt;
Nt_005 = data_dt_005.Nt;


n_02 = data_dt_02.n;
n_04 = data_dt_04.n;
n_005 = data_dt_005.n;


dt_02 = data_dt_02.dt;
dt_04 = data_dt_04.dt;
dt_005 = data_dt_005.dt;


error_list_s1_02 = data_dt_02.error_list_s1;
error_list_s2_02 = data_dt_02.error_list_s2;
error_list_s1_04 = data_dt_04.error_list_s1;
error_list_s2_04 = data_dt_04.error_list_s2;
error_list_s1_005 = data_dt_005.error_list_s1;
error_list_s2_005 = data_dt_005.error_list_s2;

m = data_dt_02.m;

%Total Momentum
%x-direction
figure
plot(error_list_s1_02(:,1),m(1)*error_list_s1_02(:,2)...
    .*error_list_s1_02(:,4)+m(2)*error_list_s2_02(:,2).*...
    error_list_s2_02(:,4)-.25,'DisplayName',['dt = ',num2str(dt_02)])
hold on 
plot(error_list_s1_04(:,1),m(1)*error_list_s1_04(:,2)...
    .*error_list_s1_04(:,4)+m(2)*error_list_s2_04(:,2).*...
    error_list_s2_04(:,4)-.25,'DisplayName',['dt = ',num2str(dt_04)])
hold on 
plot(error_list_s1_005(:,1),m(1)*error_list_s1_005(:,2)...
    .*error_list_s1_005(:,4)+m(2)*error_list_s2_005(:,2).*...
    error_list_s2_005(:,4)-.25,'DisplayName',['dt = ',num2str(dt_005)])
xlabel('$t$','Interpreter','latex')
ylabel('$P-P_0$','Interpreter','latex')
title('\textbf{Total Momentum in $v_1$ Direction}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
exportgraphics(f,'Example 5 Total Momentum v1 direction dt test.eps')


%y-direction
figure
plot(error_list_s1_02(:,1),m(1)*error_list_s1_02(:,2)...
    .*error_list_s1_02(:,5)+m(2)*error_list_s2_02(:,2).*...
    error_list_s2_02(:,5)-.25,'DisplayName',['dt = ',num2str(dt_02)])
hold on 
plot(error_list_s1_04(:,1),m(1)*error_list_s1_04(:,2)...
    .*error_list_s1_04(:,5)+m(2)*error_list_s2_04(:,2).*...
    error_list_s2_04(:,5)-.25,'DisplayName',['dt = ',num2str(dt_04)])
hold on 
plot(error_list_s1_005(:,1),m(1)*error_list_s1_005(:,2)...
    .*error_list_s1_005(:,5)+m(2)*error_list_s2_005(:,2).*...
    error_list_s2_005(:,5)-.25,'DisplayName',['dt = ',num2str(dt_005)])
xlabel('$t$','Interpreter','latex')
ylabel('$P-P_0$','Interpreter','latex')
title('\textbf{Total Momentum in $v_1$ Direction}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
exportgraphics(f,'Example 5 Total Momentum v2 direction dt test.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data set 1
Vx1 = data_dt_04.Vrx1; Vy1 = data_dt_04.Vry1;
Vx2 = data_dt_04.Vrx2; Vy2 = data_dt_04.Vry2;

W1 = data_dt_04.W1; W2 = data_dt_04.W2;

E1 = m(1)*sum(W1.*(Vx1.^2 + Vy1.^2)) + m(2)*sum(W2.*(Vx2.^2 + Vy2.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data set 2
Vx1 = data_dt_02.Vrx1; Vy1 = data_dt_02.Vry1;
Vx2 = data_dt_02.Vrx2; Vy2 = data_dt_02.Vry2;

W1 = data_dt_02.W1; W2 = data_dt_04.W2;

E2 = m(1)*sum(W1.*(Vx1.^2 + Vy1.^2)) + m(2)*sum(W2.*(Vx2.^2 + Vy2.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data set 3
Vx1 = data_dt_005.Vrx1; Vy1 = data_dt_005.Vry1;
Vx2 = data_dt_005.Vrx2; Vy2 = data_dt_005.Vry2;

W1 = data_dt_005.W1; W2 = data_dt_005.W2;

E3 = m(1)*sum(W1.*(Vx1.^2 + Vy1.^2)) + m(2)*sum(W2.*(Vx2.^2 + Vy2.^2));


% Total Energy
figure
plot(error_list_s1_005(:,1),error_list_s1_005(:,6)+error_list_s2_005(:,6)-E1, ...
    '-x','MarkerIndices',round(linspace(1,Nt_005,5)),......
    'DisplayName',['dt = ',num2str(dt_005)]);
hold on 
plot(error_list_s1_02(:,1),error_list_s1_02(:,6)+error_list_s2_02(:,6)-E2, ...
    '-*','MarkerIndices',round(linspace(1,Nt_04,5)),......
    'DisplayName',['dt = ',num2str(dt_02)]);
hold on 
plot(error_list_s1_04(:,1),error_list_s1_04(:,6)+error_list_s2_04(:,6)-E3, ...
    '-o','MarkerIndices',round(linspace(1,Nt_04,5)),......
    'DisplayName',['dt = ',num2str(dt_04)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$K-K_0$','Interpreter','latex')
title('\textbf{Total Energy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
set(hl,'Position',[0.7 0.65 0.1 0.2])
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_5_Energy_dt_test.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Energy dt test.eps')

% Total Entropy
figure
plot(error_list_s1_005(:,1),error_list_s1_005(:,8)+error_list_s2_005(:,8), ...
    '-x','MarkerIndices',round(linspace(1,Nt_005,5)),...
    'DisplayName',['dt = ',num2str(dt_005)]);
hold on 
plot(error_list_s1_04(:,1),error_list_s1_04(:,8)+error_list_s2_04(:,8), ...
    '-o','MarkerIndices',round(linspace(1,Nt_04,5)),...
    'DisplayName',['dt = ',num2str(dt_04)]);
hold on 
plot(error_list_s1_02(:,1),error_list_s1_02(:,8)+error_list_s2_02(:,8),...
    '-*','MarkerIndices',round(linspace(1,Nt_02,5)),...
    'DisplayName',['dt = ',num2str(dt_02)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^N$','Interpreter','latex')
title('\textbf{Total Entropy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_5_entropy_dt_test.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 entropy dt test.eps')




% Plots of temperature

T0_005 = data_dt_005.T_relax;
T0_005_vec = zeros(length(error_list_s2_005(:,1)),1);
T0_005_vec = T0_005_vec + T0_005;

T0_02 = data_dt_02.T_relax;
T0_02_vec = zeros(length(error_list_s2_02(:,1)),1);
T0_02_vec = T0_02_vec + T0_02;

T0_04 = data_dt_04.T_relax;
T0_04_vec = zeros(length(error_list_s2_04(:,1)),1);
T0_04vec = T0_04_vec + T0_04;

figure 
% plot(error_list_s1_005(:,1),error_list_s1_005(:,7), '-o', ...
%     'MarkerIndices',round(linspace(1,length(error_list_s1_005(:,1)),20)),...
%     'DisplayName',...
%     ['Species 1 Temperature',num2str(dt_005)]) 
% hold on
% %
% plot(error_list_s2_005(:,1),error_list_s2_005(:,7),'-square', ...
%     'MarkerIndices',round(linspace(1,length(error_list_s2_005(:,1)),20))...
%     ,'Displayname', ['Species 2 Temperature',num2str(dt_02)])
% hold on
% plot(error_list_s1_005(:,1),T0_005_vec,'DisplayName','Relaxation temperature')
% hold on
plot(error_list_s1_02(:,1),error_list_s1_02(:,7), '-o', ...
    'MarkerIndices',round(linspace(1,length(error_list_s1_02(:,1)),20)),...
    'DisplayName',...
    ['Species 1 Temperature',num2str(dt_02)]) 
hold on
plot(error_list_s2_02(:,1),error_list_s2_02(:,7),'-square', ...
    'MarkerIndices',round(linspace(1,length(error_list_s2_02(:,1)),20))...
    ,'Displayname', ['Species 2 Temperature',num2str(dt_02)])
hold on
plot(error_list_s1_02(:,1),T0_02_vec,'DisplayName','Relaxation temperature')
hold off
title('\textbf{Temperature Relaxation}','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$E^n$','Interpreter','latex')
hl = legend('Show','Location','southeast');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_5_Temperature.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Temperature_dt_test.eps')


% Plots of velocities error list 4 and 5
u0 = data_dt_02.u0;

ux_relax = data_dt_02.ux_relax;
uy_relax = data_dt_02.uy_relax;

ux0vec = zeros(length(error_list_s1_02(:,1)),1);
uy0vec = zeros(length(error_list_s1_02(:,1)),1);

ux0vec = ux0vec + ux_relax;
uy0vec = uy0vec + uy_relax;

figure 
plot(error_list_s1_02(:,1),error_list_s1_02(:,4),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_02(:,1)),20))...
    ,'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_02(:,1),error_list_s2_02(:,4),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_02(:,1)),20))...
    ,'DisplayName','Species 2 velocity ')
hold on 
plot(error_list_s1_02(:,1),ux0vec,'DisplayName','Relaxation velocity')
hold off
title('\textbf{Velocity Relaxation $v_1$ Direction}','Interpreter','latex')
ylabel('$U$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_5_Velocity_v1.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Velocity v1.eps')

figure 
plot(error_list_s1_02(:,1),error_list_s1_02(:,5),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_02(:,1)),20)),...
    'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_02(:,1),error_list_s2_02(:,5),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_02(:,1)),20)),...
    'DisplayName','Species 2 velocity')
hold on 
plot(error_list_s2_02(:,1),uy0vec,'DisplayName','Relaxation velocity')
hold off
title('\textbf{Velocity Relaxation $v_2$ Direction}','Interpreter','latex')
ylabel('$U$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_5_Velocity_v2.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Velocity v2.eps')












