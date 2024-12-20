clc
clear all
close all


%gamma = 0
%data_n_40 = load("multi_species_particle_2d_Coulomb_n_40_Example_8_gamma_0_dv1_0.2_dv2_0.2_dt_0.01tmax15.mat");

%gamma = -1
%data_n_40 = load("multi_species_particle_2d_Coulomb_n_40_Example_4_gamma_-1_dv1_0.2_dv2_0.2_dt_0.01tmax5.mat");

%gamma = -2
%data_n_40 = load("multi_species_particle_2d_Coulomb_n_60_Example_9_gamma_-2_dv1_0.13333_dv2_0.13333_dt_0.01tmax50.mat");

%gamma = -3
data_n_40 = load("multi_species_particle_2d_Coulomb_n_60_Example_9_gamma_-3_dv1_0.13333_dv2_0.13333_dt_0.01tmax50.mat");

error_list_s1_n_40 = data_n_40.error_list_s1;
error_list_s2_n_40 = data_n_40.error_list_s2;

% %Species 1 energy 
% figure 
% plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,6));
% title('Energy Species1')
% 
% %Species 2 energy 
% figure 
% plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,6));
% title('Energy Species1')

% Total Energy
figure
plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,6)+error_list_s2_n_40(:,6))
xlabel('$t$','Interpreter','latex')
ylabel('$K$','Interpreter','latex')
ylim([1.125,1.1254])
title('\textbf{Total Energy}','Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 5 Energy.eps')


T0 = data_n_40.T_relax;
T0vec = zeros(length(error_list_s2_n_40(:,1)),1);
T0vec = T0vec + T0;

m = data_n_40.m;
%Total momentum plots
%x-direction
figure
plot(error_list_s1_n_40(:,1),m(1)*error_list_s1_n_40(:,2)...
    .*error_list_s1_n_40(:,4)+m(2)*error_list_s2_n_40(:,2).*...
    error_list_s2_n_40(:,4)-.25)
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
title('\textbf{Total Momentum in $v_1$ Direction}','Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 5 Total Momentum v1 direction.eps')


%y-direction
figure
plot(error_list_s1_n_40(:,1),m(1)*error_list_s1_n_40(:,2)...
    .*error_list_s1_n_40(:,5)+m(2)*error_list_s2_n_40(:,2).*...
    error_list_s2_n_40(:,5)-.25)
% ax = gca;
% ax.YAxis.Exponent = 2;
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
title('\textbf{Total Momentum in $v_2$ Direction}','Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 5 Total Momentum v2 direction.eps')

% Plots of temperature
figure 
plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,7), '-o', ...
    'MarkerIndices',round(linspace(1,length(error_list_s1_n_40(:,1)),20)),...
    'DisplayName',...
    'Species 1 Temperature') 
hold on
plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,7),'-square', ...
    'MarkerIndices',round(linspace(1,length(error_list_s2_n_40(:,1)),20))...
    ,'Displayname', 'Species 2 Temperature')
hold on
plot(error_list_s1_n_40(:,1),T0vec,'DisplayName','Relaxation temperature')
hold off
title('\textbf{Temperature Relaxation}','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$E^n$','Interpreter','latex')
hl = legend('Show','Location','southeast');
set(hl,'Interpreter','latex')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example 5 Temperature.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Temperature.eps')


% Plots of velocities error list 4 and 5
u0 = data_n_40.u0;

ux_relax = data_n_40.ux_relax;
uy_relax = data_n_40.uy_relax;

ux0vec = zeros(length(error_list_s1_n_40(:,1)),1);
uy0vec = zeros(length(error_list_s1_n_40(:,1)),1);

ux0vec = ux0vec + ux_relax;
uy0vec = uy0vec + uy_relax;


figure 
plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,4),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_n_40(:,1)),20))...
    ,'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,4),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_n_40(:,1)),20))...
    ,'DisplayName','Species 2 velocity ')
hold on 
plot(error_list_s1_n_40(:,1),ux0vec,'DisplayName','Relaxation velocity')
hold off
title('\textbf{Velocity Relaxation $v_1$ Direction}','Interpreter','latex')
ylabel('$U$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example 5 Velocity v1.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Velocity v1.eps')



figure 
plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,5),'-o',...
    'MarkerIndices',round(linspace(1,length(error_list_s1_n_40(:,1)),20)),...
    'DisplayName','Species 1 velocity')
hold on
plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,5),'-square',...
    'MarkerIndices',round(linspace(1,length(error_list_s2_n_40(:,1)),20)),...
    'DisplayName','Species 2 velocity')
hold on 
plot(error_list_s1_n_40(:,1),uy0vec,'DisplayName','Relaxation velocity')
hold off
title('\textbf{Velocity Relaxation $v_2$ Direction}','Interpreter','latex')
ylabel('$U$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example 5 Velocity v2.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 5 Velocity v2.eps')

% Plots of entropy and total entropy 
figure 
% plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,8),'DisplayName',...
%     'Species 1 Entropy')
% hold on 
% plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,8),'DisplayName'...
%     ,'Species 2 Entropy')
% hold on  
plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,8)+error_list_s2_n_40(:,8),'DisplayName','Total Entropy')
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^N$','Interpreter','latex')
title('\textbf{Total Entropy}','Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 5 Entropy.eps')


vr1x = data_n_40.vr1x;
vr2x = data_n_40.vr2x;

vr1y = data_n_40.vr1y;
vr2y = data_n_40.vr2y;

dv1 = data_n_40.dv1;
dv2 = data_n_40.dv2;

Nr = data_n_40.Nr;
% Solution at final time
f1 = data_n_40.f1;
f2 = data_n_40.f2;

% Initial Condition
f1_0 = data_n_40.f1_0;
f2_0 = data_n_40.f2_0;


Nt = data_n_40.Nt;

figure
plot(vr1x,f1(:,Nr/2),'DisplayName',['f1(:,Nr/2) at t = ',num2str(error_list_s1_n_40(Nt,1))]) 
legend
hold on 
plot(vr1x,f1_0(:,Nr/2), 'DisplayName','Initial Condition')
title('f1(:,Nr/2)')

figure
plot(vr2x,f2(:,Nr/2),'DisplayName',['f2(:,Nr/2) at t = ',num2str(error_list_s2_n_40(Nt,1))])
legend
hold on 
plot(vr2x,f2_0(:,Nr/2),'DisplayName','Initial Condition')
title('f2(:,Nr/2)')

figure
plot(vr1y,f1(Nr/2,:),'DisplayName',['f1(Nr/2,:) at t = ',num2str(error_list_s1_n_40(Nt,1))]) 
legend
hold on 
plot(vr1y,f1_0(Nr/2,:), 'DisplayName','Initial Condition')
title('f1(Nr/2,:)')

figure
plot(vr2y,f2(Nr/2,:),'DisplayName',['f2(Nr/2,:) at t = ',num2str(error_list_s2_n_40(Nt,1))])
legend
hold on 
plot(vr2y,f2_0(Nr/2,:),'DisplayName','Initial Condition')
title('f2(Nr/2,:)')





