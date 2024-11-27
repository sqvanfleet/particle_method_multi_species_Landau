clc 
clear all
close all


data_dt_0008 = ...
    load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.00076923tmax5.mat');
data_dt_0015 = ...
    load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.0015385tmax5.mat');
data_dt_001 = ...
    load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.0010256tmax5.mat');


Nt_0008 = data_dt_0008.Nt;
Nt_0015 = data_dt_0015.Nt;
Nt_001 = data_dt_001.Nt;


n_0008 = data_dt_0008.n;
n_0015 = data_dt_0015.n;
n_001 = data_dt_001.n;


dt_0008 = data_dt_0008.dt;
dt_0015 = data_dt_0015.dt;
dt_001 = data_dt_001.dt;


error_list_s1_0008 = data_dt_0008.error_list_s1;
error_list_s2_0008 = data_dt_0008.error_list_s2;
error_list_s1_0015 = data_dt_0015.error_list_s1;
error_list_s2_0015 = data_dt_0015.error_list_s2;
error_list_s1_001 = data_dt_001.error_list_s1;
error_list_s2_001 = data_dt_001.error_list_s2;


%%%% Plot of the L2 error versus time species 1
figure
plot(error_list_s1_0008(:,1),error_list_s1_0008(:,4),'DisplayName',[' dt = ',num2str(dt_0008)]);
hold on 
plot(error_list_s1_0015(:,1),error_list_s1_0015(:,4),'DisplayName',[' dt = ',num2str(dt_0008)]);
hold on 
plot(error_list_s1_001(:,1),error_list_s1_001(:,4),'DisplayName',[' dt = ',num2str(dt_001)]);
hold off
xlabel('Time')
ylabel('L2 norm of error')
title('Species 2')
legend


%%%% Plot of the L2 error versus time species 2
figure
plot(error_list_s2_0008(:,1),error_list_s2_0008(:,4),'DisplayName',[' dt = ',num2str(dt_0008)]);
hold on 
plot(error_list_s2_0015(:,1),error_list_s2_0015(:,4),'DisplayName',[' dt = ',num2str(dt_0008)]);
hold on 
plot(error_list_s2_001(:,1),error_list_s2_001(:,4),'DisplayName',[' dt = ',num2str(dt_001)]);
hold off
xlabel('Time')
ylabel('L2 norm of error')
title('Species 2')
legend

%Total momentum plots
%x-direction
figure
plot(error_list_s1_0008(:,1),error_list_s1_0008(:,6)...
    +error_list_s2_0008(:,6),...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_0008)]);
hold on
plot(error_list_s1_001(:,1),error_list_s1_001(:,6)...
    +error_list_s2_001(:,6),...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_001)]);
hold on
plot(error_list_s1_0015(:,1),error_list_s1_0015(:,6)...
    +error_list_s2_0015(:,6), ...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_0015)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
%ax = gca;
%ax.YAxis.Exponent = -15;
%ylim([-M,M])
title('\textbf{Total Momentum in $v_2$ Direction}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
exportgraphics(f,'Example 1 Total Momentum v1 direction dt test.eps')

%Total momentum plots
%y-direction
figure
plot(error_list_s1_0008(:,1),error_list_s1_0008(:,7)...
    +error_list_s2_0008(:,7), ...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_0008)]);
hold on
plot(error_list_s1_001(:,1),error_list_s1_001(:,7)...
    +error_list_s2_001(:,7), ...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_001)]);
hold on
plot(error_list_s1_0015(:,1),error_list_s1_0015(:,7)...
    +error_list_s2_0015(:,7), ...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_0015)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
%ax = gca;
%ax.YAxis.Exponent = -15;
%ylim([-M,M])
title('\textbf{Total Momentum in $v_2$ Direction}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
exportgraphics(f,'Example_1_Total_Momentum_v2_direction_dt_test.eps')


%Weights, initial velocities, initial energy ect...

m1 = 2; m2 = 1;

Vx1 = data_dt_0008.Vrx1; Vx2 = data_dt_0008.Vrx2;
Vy1 = data_dt_0008.Vry1; Vy2 = data_dt_0008.Vry2;

W1 = data_dt_0008.W1; W2 = data_dt_0008.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 1
%energy
E_1 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_0015.Vrx1; Vx2 = data_dt_0015.Vrx2;
Vy1 = data_dt_0015.Vry1; Vy2 = data_dt_0015.Vry2;

W1 = data_dt_0015.W1; W2 = data_dt_0015.W2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 2
%energy
E_2 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vx1 = data_dt_001.Vrx1; Vx2 = data_dt_001.Vrx2;
Vy1 = data_dt_001.Vry1; Vy2 = data_dt_001.Vry2;

W1 = data_dt_001.W1; W2 = data_dt_001.W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moments for data set 3
%energy
E_3 = m1*sum(W1.*(Vx1.^2 + Vy1.^2)) + m2*sum(W2.*(Vx2.^2 + Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot of total energy
figure
plot(error_list_s1_0008(:,1),error_list_s1_0008(:,8) + error_list_s2_0008(:,8)-E_1,...
    '-x','MarkerIndices',round(linspace(1,Nt_0008,5)),...
    'DisplayName',['$\Delta t =$ ',num2str(dt_0008)]);
hold on 
plot(error_list_s1_001(:,1),error_list_s1_001(:,8) + error_list_s2_001(:,8)-E_3,...
    '-o','MarkerIndices',round(linspace(1,Nt_001,5)),...
   'DisplayName',[' $\Delta t =$ ',num2str(dt_001)]);
hold on 
plot(error_list_s1_0015(:,1),error_list_s1_0015(:,8) + error_list_s2_0015(:,8)-E_2,...
    '-*','MarkerIndices',round(linspace(1,Nt_0015,5)),...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_0015)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$K-K_0$','Interpreter','latex')
title('\textbf{Total Energy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex','Location','southwest')
%ylim([3.9999896442,3.9999896452])
% ax = gca;
% ax.YAxis.Exponent = -4;
% ytickformat('%.5f')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_1_energy_dt_test.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 1 energy dt test.eps')


%Plot of total entropy
figure
plot(error_list_s1_0008(:,1),error_list_s1_0008(:,9) + error_list_s2_0008(:,9),...
    '-x','MarkerIndices',round(linspace(1,Nt_0008,5)),...
    'DisplayName',['$\Delta t =$ ',num2str(dt_0008)]);
hold on 
plot(error_list_s1_001(:,1),error_list_s1_001(:,9) + error_list_s2_001(:,9),...
   '-o','MarkerIndices',round(linspace(1,Nt_001,5)),...
   'DisplayName',[' $\Delta t =$ ',num2str(dt_001)]);
hold on 
plot(error_list_s1_0015(:,1),error_list_s1_0015(:,9) + error_list_s2_0015(:,9),...
    '-*','MarkerIndices',round(linspace(1,Nt_0015,5)),...
    'DisplayName',[' $\Delta t =$ ',num2str(dt_0015)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^N$','Interpreter','latex')
title('\textbf{Total Entropy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex','Location','northeast')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_1_entropy_dt_test.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 1 entropy dt test.eps')











