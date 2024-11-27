clc 
clear all
close all

data_dt_0025 = load('multi_species_particle_2d_symmetric_parallel_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.0025tmax5.mat');
data_dt_005 = load('multi_species_particle_2d_symmetric_parallel_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.005tmax5.mat')';
data_dt_01 = load('multi_species_particle_2d_symmetric_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.01tmax5.mat');
data_dt_02 = load('multi_species_particle_2d_symmetric_parallel_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.02tmax5.mat');
data_dt_04 = load('multi_species_particle_2d_symmetric_parallel_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.04tmax5.mat');
data_dt_08 = load('multi_species_particle_2d_symmetric_parallel_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.079365tmax5.mat');
data_dt_16 = load('multi_species_particle_2d_symmetric_parallel_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.15625tmax5.mat');

Nt_0025 = data_dt_0025.Nt;
Nt_005 = data_dt_005.Nt;
Nt_01 = data_dt_01.Nt;
Nt_02 = data_dt_02.Nt;
Nt_04 = data_dt_04.Nt;
Nt_08 = data_dt_08.Nt;
Nt_16 = data_dt_16.Nt;

error_list_s1_0025 = data_dt_0025.error_list_s1;
error_list_s1_005 = data_dt_005.error_list_s1;
error_list_s1_01 = data_dt_01.error_list_s1;
error_list_s1_02 = data_dt_02.error_list_s1;
error_list_s1_04 = data_dt_04.error_list_s1;
error_list_s1_08 = data_dt_08.error_list_s1;
error_list_s1_16 = data_dt_16.error_list_s1;

error_list_s2_0025 = data_dt_0025.error_list_s2;
error_list_s2_005 = data_dt_005.error_list_s2;
error_list_s2_01 = data_dt_01.error_list_s2;
error_list_s2_02 = data_dt_02.error_list_s2;
error_list_s2_04 = data_dt_04.error_list_s2;
error_list_s2_08 = data_dt_08.error_list_s2;
error_list_s2_16 = data_dt_16.error_list_s2;

dt_0025 = data_dt_0025.dt;
dt_005 = data_dt_005.dt;
dt_01 = data_dt_01.dt;
dt_02 = data_dt_02.dt;
dt_04 = data_dt_04.dt;
dt_08 = data_dt_08.dt;
dt_16 = data_dt_16.dt;

%%%% Plot of the L2 error versus time species 1
figure
plot(error_list_s1_0025(:,1),error_list_s1_0025(:,4),'DisplayName',[' dt = ',num2str(dt_0025)]);
hold on 
plot(error_list_s1_005(:,1),error_list_s1_005(:,4),'DisplayName',[' dt = ',num2str(dt_005)]);
hold on 
plot(error_list_s1_01(:,1),error_list_s1_01(:,4),'DisplayName',[' dt = ',num2str(dt_01)]);
hold on 
plot(error_list_s1_02(:,1),error_list_s1_02(:,4),'DisplayName',[' dt = ',num2str(dt_02)]);
hold on 
plot(error_list_s1_04(:,1),error_list_s1_04(:,4),'DisplayName',[' dt = ',num2str(dt_04)]);
hold on 
plot(error_list_s1_08(:,1),error_list_s1_08(:,4),'DisplayName',[' dt = ',num2str(dt_08)]);
hold on 
plot(error_list_s1_16(:,1),error_list_s1_16(:,4),'DisplayName',[' dt = ',num2str(dt_16)]);
xlabel('Time')
ylabel('L2 norm of error')
title('Species 1')
legend
hold off


%%%% Plot of the L2 error versus time species 2
figure
plot(error_list_s2_0025(:,1),error_list_s2_0025(:,4),'DisplayName',[' dt = ',num2str(dt_0025)]);
hold on 
plot(error_list_s2_005(:,1),error_list_s2_005(:,4),'DisplayName',[' dt = ',num2str(dt_005)]);
hold on 
plot(error_list_s2_01(:,1),error_list_s2_01(:,4),'DisplayName',[' dt = ',num2str(dt_01)]);
hold on 
plot(error_list_s2_02(:,1),error_list_s2_02(:,4),'DisplayName',[' dt = ',num2str(dt_02)]);
hold on 
plot(error_list_s2_04(:,1),error_list_s2_04(:,4),'DisplayName',[' dt = ',num2str(dt_04)]);
hold on 
plot(error_list_s2_08(:,1),error_list_s2_08(:,4),'DisplayName',[' dt = ',num2str(dt_08)]);
hold on 
plot(error_list_s2_16(:,1),error_list_s2_16(:,4),'DisplayName',[' dt = ',num2str(dt_16)]);
xlabel('Time')
ylabel('L2 norm of error')
title('Species 2')
legend
hold off


%Total momentum plots
%x-direction
figure
plot(error_list_s1_0025(:,1),error_list_s1_0025(:,6)...
    +error_list_s2_0025(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_0025)]);
hold on
plot(error_list_s1_005(:,1),error_list_s1_005(:,6)...
    +error_list_s2_005(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_005)]);
hold on
plot(error_list_s1_01(:,1),error_list_s1_01(:,6)...
    +error_list_s2_01(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_01)]);
hold on
plot(error_list_s1_02(:,1),error_list_s1_02(:,6)...
    +error_list_s2_02(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_02)]);
hold on
plot(error_list_s1_04(:,1),error_list_s1_04(:,6)...
    +error_list_s2_04(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_04)]);
hold on
plot(error_list_s1_08(:,1),error_list_s1_08(:,6)...
    +error_list_s2_08(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_08)]);
hold on
plot(error_list_s1_08(:,1),error_list_s1_08(:,6)...
    +error_list_s2_08(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_08)]);
hold on
plot(error_list_s1_16(:,1),error_list_s1_16(:,6)...
    +error_list_s2_16(:,6)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_16)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
% ax = gca;
% ax.YAxis.Exponent = -15;
% ylim([-L,L])
title('\textbf{Total Momentum in $v_1$ Direction}','Interpreter','latex')
hl = legend('show','Location','northwest');
set(hl,'Interpreter','latex')
f = gcf;
exportgraphics(f,'Example 1b Total Momentum v1 direction dt test.eps')


%y-direction
figure
plot(error_list_s1_0025(:,1),error_list_s1_0025(:,7)...
    +error_list_s2_0025(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_0025)]);
hold on
plot(error_list_s1_005(:,1),error_list_s1_005(:,7)...
    +error_list_s2_005(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_005)]);
hold on
plot(error_list_s1_01(:,1),error_list_s1_01(:,7)...
    +error_list_s2_01(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_01)]);
hold on
plot(error_list_s1_02(:,1),error_list_s1_02(:,7)...
    +error_list_s2_02(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_02)]);
hold on
plot(error_list_s1_04(:,1),error_list_s1_04(:,7)...
    +error_list_s2_04(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_04)]);
hold on
plot(error_list_s1_08(:,1),error_list_s1_08(:,7)...
    +error_list_s2_08(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_08)]);
hold on
plot(error_list_s1_16(:,1),error_list_s1_16(:,7)...
    +error_list_s2_16(:,7)...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_16)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
%ax = gca;
%ax.YAxis.Exponent = -15;
%ylim([-M,M])
title('\textbf{Total Momentum in $v_2$ Direction}','Interpreter','latex')
hl = legend('show','Location','northwest');
set(hl,'Interpreter','latex')
f = gcf;
exportgraphics(f,'Example 1b Total Momentum v2 direction dt test.eps')


%Weights, initial velocities, initial energy ect...

m1 = 2; m2 = 1;

Vx1 = data_dt_0025.Vrx1; Vx2 = data_dt_0025.Vrx2;
Vy1 = data_dt_0025.Vry1; Vy2 = data_dt_0025.Vry2;

W1 = data_dt_0025.W1; W2 = data_dt_0025.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 1
%energy
E_1 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_005.Vrx1; Vx2 = data_dt_005.Vrx2;
Vy1 = data_dt_005.Vry1; Vy2 = data_dt_005.Vry2;

W1 = data_dt_005.W1; W2 = data_dt_005.W2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 2
%energy
E_2 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_01.Vrx1; Vx2 = data_dt_01.Vrx2;
Vy1 = data_dt_01.Vry1; Vy2 = data_dt_01.Vry2;

W1 = data_dt_01.W1; W2 = data_dt_01.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 3
%energy
E_3 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_02.Vrx1; Vx2 = data_dt_02.Vrx2;
Vy1 = data_dt_02.Vry1; Vy2 = data_dt_02.Vry2;

W1 = data_dt_02.W1; W2 = data_dt_02.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 4
%energy
E_4 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_04.Vrx1; Vx2 = data_dt_04.Vrx2;
Vy1 = data_dt_04.Vry1; Vy2 = data_dt_04.Vry2;

W1 = data_dt_04.W1; W2 = data_dt_04.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 5
%energy
E_5 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_08.Vrx1; Vx2 = data_dt_08.Vrx2;
Vy1 = data_dt_08.Vry1; Vy2 = data_dt_08.Vry2;

W1 = data_dt_08.W1; W2 = data_dt_08.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 6
%energy
E_6 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of total energy

Vx1 = data_dt_16.Vrx1; Vx2 = data_dt_16.Vrx2;
Vy1 = data_dt_16.Vry1; Vy2 = data_dt_16.Vry2;

W1 = data_dt_16.W1; W2 = data_dt_16.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 7
%energy
E_7 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(error_list_s1_0025(:,1),error_list_s1_0025(:,8) + error_list_s2_0025(:,8)-E_1,...
    '-x','MarkerIndices',round(linspace(1,Nt_0025,5))...
    ,'DisplayName',['$\Delta t =$ ',num2str(dt_0025)]);
hold on 
plot(error_list_s1_005(:,1),error_list_s1_005(:,8) + error_list_s2_005(:,8)-E_2,...
    '-square','MarkerIndices',round(linspace(1,Nt_005,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_005)]);
hold on 
plot(error_list_s1_01(:,1),error_list_s1_01(:,8) + error_list_s2_01(:,8)-E_3,...
    '-^','MarkerIndices',round(linspace(1,Nt_01,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_01)]);
hold on 
plot(error_list_s1_02(:,1),error_list_s1_02(:,8) + error_list_s2_02(:,8)-E_4,...
    '-*','MarkerIndices',round(linspace(1,Nt_02,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_02)]);
hold on 
plot(error_list_s1_04(:,1),error_list_s1_04(:,8) + error_list_s2_04(:,8)-E_5,...
    '-o','MarkerIndices',round(linspace(1,Nt_04,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_04)]);
hold on 
plot(error_list_s1_08(:,1),error_list_s1_08(:,8) + error_list_s2_08(:,8)-E_6,...
    '-diamond','MarkerIndices',round(linspace(1,Nt_08,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_08)]);
hold on 
plot(error_list_s1_16(:,1),error_list_s1_16(:,8) + error_list_s2_16(:,8)-E_7,...
    '-+','MarkerIndices',round(linspace(1,Nt_16,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_16)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$K-K_0$','Interpreter','latex')
title('\textbf{Total Energy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex','Location','northwest')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_1b_energy_dt_test.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 1b energy dt test.eps')


%Plot of total entropy
figure
plot(error_list_s1_0025(:,1),error_list_s1_0025(:,9) + error_list_s2_0025(:,9),...
    '-x','MarkerIndices',round(linspace(1,Nt_0025,5))...
    ,'DisplayName',['$\Delta t =$ ',num2str(dt_0025)]);
hold on 
plot(error_list_s1_005(:,1),error_list_s1_005(:,9) + error_list_s2_005(:,9),...
    '-square','MarkerIndices',round(linspace(1,Nt_005,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_005)]);
hold on 
plot(error_list_s1_01(:,1),error_list_s1_01(:,9) + error_list_s2_01(:,9),...
    '-^','MarkerIndices',round(linspace(1,Nt_01,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_01)]);
hold on 
plot(error_list_s1_02(:,1),error_list_s1_02(:,9) + error_list_s2_02(:,9),...
    '-*','MarkerIndices',round(linspace(1,Nt_02,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_02)]);
hold on 
plot(error_list_s1_04(:,1),error_list_s1_04(:,9) + error_list_s2_04(:,9),...
    '-o','MarkerIndices',round(linspace(1,Nt_04,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_04)]);
hold on 
plot(error_list_s1_08(:,1),error_list_s1_08(:,9) + error_list_s2_08(:,9),...
    '-diamond','MarkerIndices',round(linspace(1,Nt_08,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_08)]);
hold on 
plot(error_list_s1_16(:,1),error_list_s1_16(:,9) + error_list_s2_16(:,9),...
    '-+','MarkerIndices',round(linspace(1,Nt_16,5))...
    ,'DisplayName',[' $\Delta t =$ ',num2str(dt_16)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^{N}$','Interpreter','latex')
title('\textbf{Total Entropy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex','Location','northeast')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_1b_entropy_dt_test.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 1b entropy dt test.eps')





