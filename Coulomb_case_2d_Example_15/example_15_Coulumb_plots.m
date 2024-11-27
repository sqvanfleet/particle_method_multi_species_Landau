clc
clear all
close all


%gamma = 0
%data_n_40 = load("multi_species_particle_2d_Coulomb_n_40_Example_8_gamma_0_dv1_0.2_dv2_0.2_dt_0.01tmax15.mat");

%gamma = -1
%data_n_40 = load("multi_species_particle_2d_Coulomb_n_40_Example_4_gamma_-1_dv1_0.2_dv2_0.2_dt_0.01tmax5.mat");

%gamma = -2
data_n_40_same = load("multi_species_particle_2d_Coulomb_n_60_Example_8_gamma_-3_dv1_0.13333_dv2_0.13333_dt_0.01tmax50.mat");

%gamma = -3
data_n_40_different = load("multi_species_particle_2d_Coulomb_n_60_Example_15_gamma_-3_dv1_0.083333_dv2_0.11826_dt_0.01tmax50.mat");

error_list_s1_same = data_n_40_same.error_list_s1;
error_list_s2_same = data_n_40_same.error_list_s2;

error_list_s1_different = data_n_40_different.error_list_s1;
error_list_s2_different = data_n_40_different.error_list_s2;

% %Species 1 energy 
% figure 
% plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,6));
% title('Energy Species1')
% 
% %Species 2 energy 
% figure 
% plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,6));
% title('Energy Species1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = data_n_40_same.m;
%Initial temperature same domain
vx1 = data_n_40_same.Vrx1; vy1 = data_n_40_same.Vry1;
vx2 = data_n_40_same.Vrx2; vy2 = data_n_40_same.Vry2;

W1 = data_n_40_same.W1; W2 = data_n_40_same.W2;

E_same = m(1)*sum(W1.*(vx1.^2 + vy1.^2)) + m(2)*sum(W2.*(vx2.^2 + vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = data_n_40_different.m;
%Initial temperature same domain
vx1 = data_n_40_different.Vrx1; vy1 = data_n_40_different.Vry1;
vx2 = data_n_40_different.Vrx2; vy2 = data_n_40_different.Vry2;

W1 = data_n_40_different.W1; W2 = data_n_40_different.W2;

E_different = m(1)*sum(W1.*(vx1.^2 + vy1.^2)) + m(2)*sum(W2.*(vx2.^2 + vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%



% Total Energy
figure
plot(error_list_s1_same(:,1),error_list_s1_same(:,6)+error_list_s2_same(:,6)-E_same,...
    '-o','MarkerIndices',round(linspace(1,length(error_list_s1_same(:,1)),20)),...
    'DisplayName','Total Energy Same Domain')
hold on 
plot(error_list_s1_different(:,1),error_list_s1_different(:,6)+error_list_s2_different(:,6)-E_different,...
    '-square','MarkerIndices',round(linspace(1,length(error_list_s1_different(:,1)),20)),...
    'DisplayName','Total Energy Different Domain')
xlabel('$t$','Interpreter','latex')
ylabel('$K-K_0$','Interpreter','latex')
ylim(1e-4*[0,10])
title('\textbf{Total Energy}','Interpreter','latex')
hl = legend('Show','Location','southeast');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_6_Energy.pdf','-dpdf','-r0')
% exportgraphics(f,'Example_6_Energy.eps')


T0_same = data_n_40_same.T_relax;
T0vec_same = zeros(length(error_list_s2_same(:,1)),1);
T0vec_same = T0vec_same + T0_same;

% Plots of temperature same computational domains 
figure 
plot(error_list_s1_same(:,1),error_list_s1_same(:,7), '-o', ...
    'MarkerIndices',round(linspace(1,length(error_list_s1_same(:,1)),20)),...
    'DisplayName',...
    'Species 1 Temperature') 
hold on
plot(error_list_s2_same(:,1),error_list_s2_same(:,7),'-square', ...
    'MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20))...
    ,'Displayname', 'Species 2 Temperature')
hold on
plot(error_list_s1_same(:,1),T0vec_same,'DisplayName','Relaxation temperature')
hold off
title('\textbf{Temperature Relaxation}','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
ylim([0.12,0.32])
hl = legend('Show','Location','southeast');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_6_Temperature_same_domain.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Temperature same domain.eps','Resolution',300)

T0_different = data_n_40_different.T_relax;
T0vec_different = zeros(length(error_list_s2_different(:,1)),1);
T0vec_different = T0vec_different + T0_different;

% Plots of temperature different computational domains
figure 
plot(error_list_s1_different(:,1),error_list_s1_different(:,7), '-o', ...
    'MarkerIndices',round(linspace(1,length(error_list_s1_different(:,1)),20)),...
    'DisplayName',...
    'Species 1 Temperature') 
hold on
plot(error_list_s2_different(:,1),error_list_s2_different(:,7),'-square', ...
    'MarkerIndices',round(linspace(1,length(error_list_s2_different(:,1)),20))...
    ,'Displayname', 'Species 2 Temperature')
hold on
plot(error_list_s1_different(:,1),T0vec_different,'DisplayName','Relaxation temperature')
hold off
title('\textbf{Temperature Relaxation}','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
hl = legend('Show','Location','southeast');
set(hl,'Interpreter','latex')
ylim([0.1,0.35])
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_6_Temperature_different_domain.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Temperature different domain.eps','Resolution',300)

%Velocity relaxation plots
%x direction
u0 = data_n_40_same.u0;

ux_relax = data_n_40_same.ux_relax;
uy_relax = data_n_40_same.uy_relax;

ux0vec = zeros(length(error_list_s1_same(:,1)),1);
uy0vec = zeros(length(error_list_s1_same(:,1)),1);

ux0vec = ux0vec + ux_relax;
uy0vec = uy0vec + uy_relax;


figure 
plot(error_list_s1_same(:,1),error_list_s1_same(:,4),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20))...
    ,'DisplayName','Species 1 velocity same')
hold on
plot(error_list_s2_same(:,1),error_list_s2_same(:,4),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20))...
    ,'DisplayName','Species 2 velocity same')
hold on 
% plot(error_list_s1_different(:,1),error_list_s1_different(:,4),'-o',...
%     'MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20))...
%     ,'DisplayName','Species 1 velocity different')
% hold on
% plot(error_list_s2_different(:,1),error_list_s2_different(:,4),'-square',...
%     'MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20))...
%     ,'DisplayName','Species 2 velocity different')
% hold on 
plot(error_list_s1_same(:,1),ux0vec,'DisplayName','Relaxation velocity')
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
print(f,'Example_6_Velocity_v1_same_domain.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Velocity v1.eps')

%y direction

figure 
plot(error_list_s1_same(:,1),error_list_s1_same(:,5),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s1_same(:,1)),20)),...
    'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_same(:,1),error_list_s2_same(:,5),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20)),...
    'DisplayName','Species 2 velocity')
hold on 
plot(error_list_s1_same(:,1),uy0vec,'DisplayName','Relaxation velocity')
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
print(f,'Example_6_Velocity_v2_same_domain.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Velocity v2 same domain.eps')

%Velocity relaxation plots
%x direction



figure 
plot(error_list_s1_different(:,1),error_list_s1_different(:,4),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_different(:,1)),20))...
    ,'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_different(:,1),error_list_s2_different(:,4),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_different(:,1)),20))...
    ,'DisplayName','Species 2 velocity ')
hold on 
plot(error_list_s1_different(:,1),ux0vec,'DisplayName','Relaxation velocity')
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
print(f,'Example_6_Velocity_v1_different_domain.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Velocity v1 different domain.eps')

%y direction

figure 
plot(error_list_s1_different(:,1),error_list_s1_different(:,5),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s1_different(:,1)),20)),...
    'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_different(:,1),error_list_s2_different(:,5),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_different(:,1)),20)),...
    'DisplayName','Species 2 velocity')
hold on 
plot(error_list_s1_different(:,1),uy0vec,'DisplayName','Relaxation velocity')
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
print(f,'Example_6_Velocity_v2_different_domain.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Velocity v2 different domain.eps')

% Plots of entropy and total entropy 
figure  
plot(error_list_s1_same(:,1),error_list_s1_same(:,8)+error_list_s2_same(:,8),...
    '-o','MarkerIndices',round(linspace(1,length(error_list_s2_same(:,1)),20)),...
    'DisplayName','Total Entropy Same Domain')
hold on 
plot(error_list_s1_different(:,1),error_list_s1_different(:,8)+error_list_s2_different(:,8),...
    '-square','MarkerIndices',round(linspace(1,length(error_list_s2_different(:,1)),20)),...
    'DisplayName','Total Entropy Different Domain')
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
print(f,'Example_6_Entropy.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 6 Entropy.eps')


