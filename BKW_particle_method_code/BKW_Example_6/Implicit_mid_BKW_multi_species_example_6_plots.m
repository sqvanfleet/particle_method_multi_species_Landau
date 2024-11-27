clc 
clear all
close all


%data_n_20 = load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_20_Example_6_dv1_0.3_dv2_0.4_dt_0.01tmax5.mat');
data_n_40 = load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_40_Example_6_dv1_0.15_dv2_0.2_dt_0.0025tmax5.mat');
data_n_45 = load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_45_Example_6_dv1_0.13333_dv2_0.17778_dt_0.002tmax5.mat');
data_n_50 = load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_50_Example_6_dv1_0.12_dv2_0.16_dt_0.0015385tmax5.mat');
data_n_55 = load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_55_Example_6_dv1_0.10909_dv2_0.14545_dt_0.00125tmax5.mat');
data_n_60 = load('implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_60_Example_6_dv1_0.1_dv2_0.13333_dt_0.001tmax5.mat');

%Nt_20 = data_n_20.Nt;
Nt_40 = data_n_40.Nt;
Nt_45 = data_n_45.Nt;
Nt_50 = data_n_50.Nt;
Nt_55 = data_n_55.Nt;
Nt_60 = data_n_60.Nt;

%n_20 = data_n_20.n;
n_40 = data_n_40.n;
n_45 = data_n_45.n;
n_50 = data_n_50.n;
n_55 = data_n_55.n;
n_60 = data_n_60.n;


%dt_20 = data_n_20.dt;
dt_40 = data_n_40.dt;
dt_45 = data_n_45.dt;
dt_50 = data_n_50.dt;
dt_55 = data_n_55.dt;
dt_60 = data_n_60.dt;

%error_list_s1_20 = data_n_20.error_list_s1;
%error_list_s2_20 = data_n_20.error_list_s2;
error_list_s1_40 = data_n_40.error_list_s1;
error_list_s2_40 = data_n_40.error_list_s2;
error_list_s1_45 = data_n_45.error_list_s1;
error_list_s2_45 = data_n_45.error_list_s2;
error_list_s1_50 = data_n_50.error_list_s1;
error_list_s2_50 = data_n_50.error_list_s2;
error_list_s1_55 = data_n_55.error_list_s1;
error_list_s2_55 = data_n_55.error_list_s2;
error_list_s1_60 = data_n_60.error_list_s1;
error_list_s2_60 = data_n_60.error_list_s2;


%%%% Plot of the L2 error versus time species 1
figure()
% plot(error_list_s1_20(:,1),error_list_s1_20(:,4),'DisplayName',['n = ',num2str(n_20),' dt = ',num2str(dt_20)]);
% hold on 
plot(error_list_s1_40(:,1),error_list_s1_40(:,4),'DisplayName',['n = ',num2str(n_40),' dt = ',num2str(dt_40)]);
hold on 
plot(error_list_s1_45(:,1),error_list_s1_45(:,4),'DisplayName',['n = ',num2str(n_45),' dt = ',num2str(dt_45)]);
hold on 
plot(error_list_s1_50(:,1),error_list_s1_50(:,4),'DisplayName',['n = ',num2str(n_50),' dt = ',num2str(dt_50)]);
hold on 
plot(error_list_s1_55(:,1),error_list_s1_55(:,4),'DisplayName',['n = ',num2str(n_55),' dt = ',num2str(dt_55)]);
hold on 
plot(error_list_s1_60(:,1),error_list_s1_60(:,4),'DisplayName',['n = ',num2str(n_55),' dt = ',num2str(dt_60)]);
hold off
xlabel('Time')
ylabel('L2 norm of error')
title('Species 1')
legend



%%%% Plot of the L2 error versus time species 2
figure
% plot(error_list_s2_20(:,1),error_list_s2_20(:,4),'DisplayName',['n = ',num2str(n_20),' dt = ',num2str(dt_20)]);
% hold on 
plot(error_list_s2_40(:,1),error_list_s2_40(:,4),'DisplayName',['n = ',num2str(n_40),' dt = ',num2str(dt_40)]);
hold on 
plot(error_list_s2_45(:,1),error_list_s2_45(:,4),'DisplayName',['n = ',num2str(n_45),' dt = ',num2str(dt_45)]);
hold on 
plot(error_list_s2_50(:,1),error_list_s2_50(:,4),'DisplayName',['n = ',num2str(n_50),' dt = ',num2str(dt_50)]);
hold on 
plot(error_list_s2_55(:,1),error_list_s2_55(:,4),'DisplayName',['n = ',num2str(n_55),' dt = ',num2str(dt_55)]);
hold on 
plot(error_list_s2_60(:,1),error_list_s2_60(:,4),'DisplayName',['n = ',num2str(n_55),' dt = ',num2str(dt_60)]);
xlabel('Time')
ylabel('L2 norm of error')
title('Species 2')
legend
hold off




% %%%% Log Log plot of the L1, L2, and Linfty error to display order of
% %%%% convergence
P = 5; % Number of data sets
Linf_error_s1 = zeros(P,1);
L1_error_s1 = zeros(P,1);
L2_error_s1 = zeros(P,1);
Linf_error_s2 = zeros(P,1);
L1_error_s2 = zeros(P,1);
L2_error_s2 = zeros(P,1);
h_list_s1 = zeros(P,1);
h_list_s2 = zeros(P,1);



% % The first column is h
h_list_s1(1) = data_n_40.dv1;
h_list_s1(2) = data_n_45.dv1;
h_list_s1(3) = data_n_50.dv1;
h_list_s1(4) = data_n_55.dv1;
h_list_s1(5) = data_n_60.dv1;

% 
% % The first column is h
h_list_s2(1) = data_n_40.dv2;
h_list_s2(2) = data_n_45.dv2;
h_list_s2(3) = data_n_50.dv2;
h_list_s2(4) = data_n_55.dv2;
h_list_s2(5) = data_n_60.dv2;


% %the second column of error_list is the Linf error
Linf_error_s1(1) = error_list_s1_40(Nt_40,2);
Linf_error_s1(2) = error_list_s1_45(Nt_45,2);
Linf_error_s1(3) = error_list_s1_50(Nt_50,2);
Linf_error_s1(4) = error_list_s1_55(Nt_55,2);
Linf_error_s1(5) = error_list_s1_60(Nt_60,2);


% %the third column of error_list is the L1 error
L1_error_s1(1) = error_list_s1_40(Nt_40,3);
L1_error_s1(2) = error_list_s1_45(Nt_45,3);
L1_error_s1(3) = error_list_s1_50(Nt_50,3);
L1_error_s1(4) = error_list_s1_55(Nt_55,3);
L1_error_s1(5) = error_list_s1_60(Nt_60,3);

% %the forth column of error_list is the L2 error
L2_error_s1(1) = error_list_s1_40(Nt_40,4);
L2_error_s1(2) = error_list_s1_45(Nt_45,4);
L2_error_s1(3) = error_list_s1_50(Nt_50,4);
L2_error_s1(4) = error_list_s1_55(Nt_55,4);
L2_error_s1(5) = error_list_s1_60(Nt_60,4);

% %the second column of error_list is the Linf error
Linf_error_s2(1) = error_list_s2_40(Nt_40,2);
Linf_error_s2(2) = error_list_s2_45(Nt_45,2);
Linf_error_s2(3) = error_list_s2_50(Nt_50,2);
Linf_error_s2(4) = error_list_s2_55(Nt_55,2);
Linf_error_s2(5) = error_list_s2_60(Nt_60,2);


% %the third column of error_list is the L1 error
L1_error_s2(1) = error_list_s2_40(Nt_40,3);
L1_error_s2(2) = error_list_s2_45(Nt_45,3);
L1_error_s2(3) = error_list_s2_50(Nt_50,3);
L1_error_s2(4) = error_list_s2_55(Nt_55,3);
L1_error_s2(5) = error_list_s2_60(Nt_60,3);


% %the forth column of error_list is the L2 error
% 
L2_error_s2(1) = error_list_s2_40(Nt_40,4);
L2_error_s2(2) = error_list_s2_45(Nt_45,4);
L2_error_s2(3) = error_list_s2_50(Nt_50,4);
L2_error_s2(4) = error_list_s2_55(Nt_55,4);
L2_error_s2(5) = error_list_s2_60(Nt_60,4);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose which data set
lower_limit = 1;
upper_limit = P;
ind = lower_limit:upper_limit;




x_s1 = ones(length(h_list_s1(ind)),2);
x_s1(:,2) = log(h_list_s1(ind));
% 
x_s2 = ones(length(h_list_s2(ind)),2);
x_s2(:,2) = log(h_list_s2(ind));
% 
% 
Linf_roc_s1 = x_s1\log(Linf_error_s1(ind)); %slope using least squares fitting
L1_roc_s1 = x_s1\log(L1_error_s1(ind)); %slope using least squares fitting
L2_roc_s1 = x_s1\log(L2_error_s1(ind)); %slope using a least squares fitting
% 
Linf_roc_s2 = x_s2\log(Linf_error_s2(ind)); %slope using least squares fitting
L1_roc_s2 = x_s2\log(L1_error_s2(ind)); %slope using least squares fitting
L2_roc_s2 = x_s2\log(L2_error_s2(ind)); %slope using a least squares fitting
% 
% 
% %Rate of convergence plot species 1
figure
% plot(-log(h_list_s1),-log(L2_error_s1),'.','markersize',8,'DisplayName','L2error')
% hold on 
% plot(-log(h_list_s1),-log(L1_error_s1),'.','markersize',8,'DisplayName','L1error')
% hold on
% plot(-log(h_list_s1),-log(Linf_error_s1),'.','markersize',8,'DisplayName','Linferror')
% hold on
plot(-log(h_list_s1(ind)),-L2_roc_s1(2)*log(h_list_s1(ind))-L2_roc_s1(1),'-o','color','r','DisplayName',['$L^2$ error slope $=$ ',num2str(L2_roc_s1(2))])
hold on
plot(-log(h_list_s1(ind)),-L1_roc_s1(2)*log(h_list_s1(ind))-L1_roc_s1(1),'-square','color','b','DisplayName',['$L^1$ error slope $=$ ',num2str(L1_roc_s1(2))])
hold on
plot(-log(h_list_s1(ind)),-Linf_roc_s1(2)*log(h_list_s1(ind))-Linf_roc_s1(1),'-^','color','g','DisplayName',['$L^{\infty}$ error slope $=$ ',num2str(Linf_roc_s1(2))])
hold on
plot(-log(h_list_s1(ind)),-2*log(h_list_s1(ind))-.7,'color','k','LineWidth',2,'DisplayName','Reference Slope = 2')
hold off
xlabel('$\log{h_1}$','Interpreter','latex')
ylabel('$\log(\mbox{error})$','Interpreter','latex')
title('$\textbf{Species 1}$','Interpreter','latex')
hl = legend('show','Location','northwest');
set(hl,'Interpreter','latex')
xlim([1.85,2.35])
ylim([2.5,5.0])
f = gcf;
set(f,'Units','Inches');
fontsize(f,20,"points")
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_1_ROC_species_1.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 1 ROC species 1.eps')

% 
% %Rate of convergence plot species 2
figure
% plot(-log(h_list_s2),-log(L2_error_s2),'.','markersize',8,'DisplayName','L2error')
% hold on 
% plot(-log(h_list_s2),-log(L1_error_s2),'.','markersize',8,'DisplayName','L1error')
% hold on
% plot(-log(h_list_s2),-log(Linf_error_s2),'.','markersize',8,'DisplayName','Linferror')
% hold on
plot(-log(h_list_s2),-L2_roc_s2(2)*log(h_list_s2)-L2_roc_s2(1),'-o',...
    'color','r','DisplayName',['$L^2$ error slope $=$ ',num2str(L2_roc_s2(2))])
hold on
plot(-log(h_list_s2),-L1_roc_s2(2)*log(h_list_s2)-L1_roc_s2(1),'-square',...
    'color','b','DisplayName',['$L^2$ error slope $=$ ',num2str(L1_roc_s2(2))])
hold on
plot(-log(h_list_s2),-Linf_roc_s2(2)*log(h_list_s2)-Linf_roc_s2(1),'-^',...
    'color','g','DisplayName',['$L^2$ error slope $=$ ',num2str(Linf_roc_s2(2))])
hold on
plot(-log(h_list_s2),-2*log(h_list_s2),'LineWidth',2,'color','k',...
    'DisplayName','Reference Slope = 2')
hold off
xlim([1.55,2.05])
ylim([2.6,5.0])
xlabel('$\log{h_2}$','Interpreter','latex')
ylabel('$\log(\mbox{error})$','Interpreter','latex')
title('$\textbf{Species 2}$','Interpreter','latex')
hl = legend('show','Location','northwest');
set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,"points")
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_1_ROC_species_2.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 1 ROC species 2.eps')



%Plot of total energy
figure
plot(error_list_s1_40(:,1),error_list_s1_40(:,8) + error_list_s2_40(:,8),...
    '-o','MarkerIndices',round(linspace(1,Nt_40,5))...
    ,'DisplayName',['$n = $ ',num2str(n_40),' $\Delta t =$ ',num2str(dt_40)]);
hold on 
plot(error_list_s1_45(:,1),error_list_s1_45(:,8) + error_list_s2_45(:,8),...
    '-square','MarkerIndices',round(linspace(1,Nt_45,5))...
    ,'DisplayName',['$n = $ ',num2str(n_45),' $\Delta t =$ ',num2str(dt_45)]);
hold on 
plot(error_list_s1_50(:,1),error_list_s1_50(:,8) + error_list_s2_50(:,8),...
    '-^','MarkerIndices',round(linspace(1,Nt_50,5))...
    ,'DisplayName',['$n = $ ',num2str(n_50),' $\Delta t =$ ',num2str(dt_50)]);
hold on 
plot(error_list_s1_55(:,1),error_list_s1_55(:,8) + error_list_s2_55(:,8),...
    '-*','MarkerIndices',round(linspace(1,Nt_55,5))...
    ,'DisplayName',['$n = $ ',num2str(n_55),' $\Delta t =$ ',num2str(dt_55)]);
hold on 
plot(error_list_s1_60(:,1),error_list_s1_60(:,8) + error_list_s2_60(:,8),...
    '-x','MarkerIndices',round(linspace(1,Nt_60,5))...
    ,'DisplayName',['$n = $ ',num2str(n_60),' $\Delta t =$ ',num2str(dt_60)]);
hold off
% ylim([2.])
% ylim([2.9,3.1])
xlabel('$t$','Interpreter','latex')
ylabel('$K$','Interpreter','latex')
%ytickformat('%.7f')
%ylim([3.9999894,3.99999015])
title('\textbf{Total Energy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
% set(fenergy,'Units','pixels')
% set(fenergy,'Position',[440 378 560 420])
exportgraphics(f,'Example 1 energy.eps')



%Plot of total entropy
figure
plot(error_list_s1_40(:,1), error_list_s1_40(:,9) + error_list_s2_40(:,9),...
    '-o','MarkerIndices',round(linspace(1,Nt_40,5)),...
    'DisplayName',['$n = $',' ',num2str(n_40)]);
hold on 
plot(error_list_s1_45(:,1), error_list_s1_45(:,9) + error_list_s2_45(:,9),...
    '-square','MarkerIndices',round(linspace(1,Nt_45,5)),...
    'DisplayName',['$n = $',' ',num2str(n_45)]);
hold on 
plot(error_list_s1_50(:,1), error_list_s1_50(:,9) + error_list_s2_50(:,9),...
    '-^','MarkerIndices',round(linspace(1,Nt_50,5)),...
    'DisplayName',['$n = $',' ',num2str(n_50)]);
hold on 
plot(error_list_s1_55(:,1), error_list_s1_55(:,9) + error_list_s2_55(:,9),...
    '-*','MarkerIndices',round(linspace(1,Nt_55,5)),...
    'DisplayName',['$n = $',' ',num2str(n_55)]);
hold on 
plot(error_list_s1_60(:,1), error_list_s1_60(:,9) + error_list_s2_60(:,9),...
    '-x','MarkerIndices',round(linspace(1,Nt_60,5)),...
    'DisplayName',['$n = $ ',num2str(n_60)]);
hold off
% ylim([2.9,3.1])
xlabel('$t$','Interpreter','latex')
ylabel('$E^{N}$','Interpreter','latex')
title('\textbf{Total Entropy}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 1 entropy.eps')

L = 10*max(abs(error_list_s1_60(:,6)+error_list_s2_60(:,6)));
M = 10*max(abs(error_list_s1_60(:,7)+error_list_s2_60(:,6)));

%Total momentum plots
%x-direction
figure
% plot(error_list_s1_40(:,1),error_list_s1_40(:,6)...
%     +error_list_s2_40(:,6),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_40),' $\Delta t =$ ',num2str(dt_40)]);
% hold on
% plot(error_list_s1_45(:,1),error_list_s1_45(:,6)...
%     +error_list_s2_45(:,6),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_45),' $\Delta t =$ ',num2str(dt_45)]);
% hold on
% plot(error_list_s1_50(:,1),error_list_s1_50(:,6)...
%     +error_list_s2_50(:,6),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_50),' $\Delta t =$ ',num2str(dt_50)]);
% hold on
% plot(error_list_s1_55(:,1),error_list_s1_55(:,6)...
%     +error_list_s2_55(:,6),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_55),' $\Delta t =$ ',num2str(dt_55)]);
% hold on
plot(error_list_s1_60(:,1),error_list_s1_60(:,6)...
    +error_list_s2_60(:,6),'-x','MarkerIndices',round(linspace(1,Nt_60,5))...
    ,'DisplayName',['$n = $ ',num2str(n_60),' $\Delta t =$ ',num2str(dt_60)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
% ax = gca;
% ax.YAxis.Exponent = -15;
% ylim([-L,L])
title('\textbf{Total Momentum in $v_1$ Direction}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 1 Total Momentum v1 direction.eps')

%Total momentum plots
%y-direction
figure
% plot(error_list_s1_40(:,1),error_list_s1_40(:,7)...
%     +error_list_s2_40(:,7),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_40),' $\Delta t =$ ',num2str(dt_40)]);
% hold on
% plot(error_list_s1_45(:,1),error_list_s1_45(:,7)...
%     +error_list_s2_45(:,7),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_45),' $\Delta t =$ ',num2str(dt_45)]);
% hold on
% plot(error_list_s1_50(:,1),error_list_s1_50(:,7)...
%     +error_list_s2_50(:,7),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_50),' $\Delta t =$ ',num2str(dt_50)]);
% hold on
% plot(error_list_s1_55(:,1),error_list_s1_55(:,7)...
%     +error_list_s2_55(:,7),'-*','MarkerIndices',round(linspace(1,Nt_60,5))...
%     ,'DisplayName',['$n = $ ',num2str(n_55),' $\Delta t =$ ',num2str(dt_55)]);
% hold on
plot(error_list_s1_60(:,1),error_list_s1_60(:,7)...
    +error_list_s2_60(:,7),'-x','MarkerIndices',round(linspace(1,Nt_60,5))...
    ,'DisplayName',['$n = $ ',num2str(n_60),' $\Delta t =$ ',num2str(dt_60)]);
hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
% ax = gca;
% ax.YAxis.Exponent = -15;
% ylim([-M,M])
title('\textbf{Total Momentum in $v_2$ Direction}','Interpreter','latex')
hl = legend('show');
set(hl,'Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 1 Total Momentum v2 direction.eps')


%Plot of numerical and exact solutions at final time

% vr1_60 = data_n_60.vr1; vr2_60 = data_n_60.vr2;
% vr1_80 = data_n_80.vr1; vr2_80 = data_n_80.vr2;
% vr1_100 = data_n_100.vr1; vr2_100 = data_n_100.vr2;
% vr1_120 = data_n_120.vr1; vr2_120 = data_n_120.vr2;
% % vr1_150 = data_n_150.vr1; vr2_150 = data_n_150.vr2;
% 
% Nr_60 = data_n_60.Nr; Nr_80 = data_n_80.Nr; 
% Nr_100 = data_n_100.Nr;
% Nr_120 = data_n_120.Nr; 
% % Nr_150 = data_n_150.Nr; 
% 
% n_60 = data_n_60.n; n_80 = data_n_80.n; 
% n_100 = data_n_100.n;
% n_120 = data_n_120.n; 
% % n_150 = data_n_150.n; 
% 
% f1_60 = data_n_60.f1; f2_60 = data_n_60.f2;
% f1_80 = data_n_80.f1; f2_80 = data_n_80.f2;
% f1_100 = data_n_100.f1; f2_100 = data_n_100.f2;
% f1_120 = data_n_120.f1; f2_120 = data_n_120.f2;
% % f1_150 = data_n_150.f1; f2_150 = data_n_150.f2;
% 
% f1_exact_120 = data_n_120.f1_exact; f2_exact_120 = data_n_120.f2_exact;
% 
% figure
% plot(vr1_60,f1_60(:,Nr_60/2),'DisplayName',['Particle Solution n = ',num2str(n_60)])
% hold on
% plot(vr1_80,f1_80(:,Nr_80/2),'DisplayName',['Particle Solution n = ',num2str(n_80)])
% hold on
% plot(vr1_100,f1_100(:,Nr_100/2),'DisplayName',['Particle Solution n = ',num2str(n_100)])
% hold on
% plot(vr1_120,f1_120(:,Nr_120/2),'DisplayName',['Particle Solution n = ',num2str(n_120)])
% % hold on
% % plot(vr1_150,f1_150(:,Nr_150/2),'DisplayName',['Particle Solution n = ',num2str(n_150)])
% hold on 
% plot(vr1_120,f1_exact_120(:,Nr_120/2),'DisplayName',['Exact Solution n = ',num2str(n_120)])
% title('Species 1 at final time')
% legend
% hold off
% 
% figure
% plot(vr2_60,f2_60(:,Nr_60/2),'DisplayName',['Particle Solution n = ',num2str(n_60)])
% hold on
% plot(vr2_80,f2_80(:,Nr_80/2),'DisplayName',['Particle Solution n = ',num2str(n_80)])
% hold on
% plot(vr2_100,f2_100(:,Nr_100/2),'DisplayName',['Particle Solution n = ',num2str(n_100)])
% hold on
% plot(vr2_120,f2_120(:,Nr_120/2),'DisplayName',['Particle Solution n = ',num2str(n_120)])
% % hold on
% % plot(vr2_150,f2_150(:,Nr_150/2),'DisplayName',['Particle Solution n = ',num2str(n_150)])
% hold on 
% plot(vr2_120,f2_exact_120(:,Nr_120/2),'DisplayName',['Exact Solution n = ',num2str(n_120)])
% title('Species 2 at final time')
% legend
% hold off
% 
% 
% 
