clc 
clear all
close all


myVars = {'error_list_s1','error_list_s2','dv1','dv2','Nt','total_moments',...
    'f1','f2','f1_exact','f2_exact','vr1','vr2','Nr','n','dt','m'};


data_dt_001 = load('multi_species_particle_2d_symmetric_n_60_Example_11_dv1_0.013333_dv2_0.13333_dt_0.001tmax5.mat',myVars{:});
% data_dt_0025 = load('multi_species_particle_2d_symmetric_n_40_Example_11_dv1_0.02_dv2_0.2_dt_0.0025tmax5.mat',myVars{:});
% data_dt_0005 = load('multi_species_particle_2d_symmetric_n_40_Example_11_dv1_0.02_dv2_0.2_dt_0.0005tmax5.mat',myVars{:});
data_dt_00025 = load('multi_species_particle_2d_symmetric_n_60_Example_11_dv1_0.013333_dv2_0.13333_dt_0.00025tmax5.mat',myVars{:});
data_dt_0001667 = load('multi_species_particle_2d_symmetric_n_100_Example_11_dv1_0.008_dv2_0.08_dt_0.00016667tmax5.mat',myVars{:});
data_dt_000125 = load('multi_species_particle_2d_symmetric_n_100_Example_11_dv1_0.008_dv2_0.08_dt_0.000125tmax5.mat',myVars{:});
data_dt_0001 = load('multi_species_particle_2d_symmetric_parallel_n_100_Example_11_dv1_0.008_dv2_0.08_dt_0.0001tmax5.mat');

% Nt_001 = data_dt_001.Nt;
% Nt_0025 = data_dt_0025.Nt;
% Nt_0005 = data_dt_0005.Nt;
% Nt_00025 = data_dt_00025.Nt;
Nt_0001667 = data_dt_0001667.Nt;
Nt_000125 = data_dt_000125.Nt;
Nt_0001 = data_dt_0001.Nt;


error_list_s1_001 = data_dt_001.error_list_s1;
% error_list_s1_0025 = data_dt_0025.error_list_s1;
% error_list_s1_0005 = data_dt_0005.error_list_s1;
error_list_s1_00025 = data_dt_00025.error_list_s1;
error_list_s1_0001667 = data_dt_0001667.error_list_s1;
error_list_s1_000125 = data_dt_000125.error_list_s1;
error_list_s1_0001 = data_dt_0001.error_list_s1;


error_list_s2_001 = data_dt_001.error_list_s2;
% error_list_s2_0025 = data_dt_0025.error_list_s2;
% error_list_s2_0005 = data_dt_0005.error_list_s2;
error_list_s2_00025 = data_dt_00025.error_list_s2;
error_list_s2_0001667 = data_dt_0001667.error_list_s2;
error_list_s2_000125 = data_dt_000125.error_list_s2;
error_list_s2_0001 = data_dt_0001.error_list_s2;


% column 1 of error list is a list of times
% column 4 of error list is a list of L2 error


%%%% Plot of the L2 error versus time species 1

figure
%plot(error_list_s1_001(:,1),error_list_s1_001(:,4),'DisplayName','dt = 0.001');
% hold on
% plot(error_list_s1_0025(:,1),error_list_s1_0025(:,4),'DisplayName','dt = 0.0025');
% hold on
% plot(error_list_s1_0005(:,1),error_list_s1_0005(:,4),'DisplayName','dt = 0.0005');
%hold on
plot(error_list_s1_0001667(:,1),error_list_s1_0001667(:,4),'DisplayName','dt = 0.0001667');
hold on 
plot(error_list_s1_000125(:,1),error_list_s1_000125(:,4),'DisplayName','dt = 0.000125');
xlabel('$t$','Interpreter','latex')
ylabel('$L^2$ norm of error','Interpreter','latex')
title('\textbf{Species 1}','Interpreter','latex')
legend
hold off



%%%% Plot of the L2 error versus time species 2

figure
%plot(error_list_s2_001(:,1),error_list_s2_001(:,4),'DisplayName','dt = 0.001');
% hold on
% plot(error_list_s2_0025(:,1),error_list_s2_0025(:,4),'DisplayName','dt = 0.0025');
% hold on
% plot(error_list_s2_0005(:,1),error_list_s2_0005(:,4),'DisplayName','dt = 0.0005');
%hold on
plot(error_list_s2_0001667(:,1),error_list_s2_0001667(:,4),'DisplayName','dt = 0.0001667');
hold on 
plot(error_list_s2_000125(:,1),error_list_s2_000125(:,4),'DisplayName','dt = 0.000125');
xlabel('Time')
ylabel('L2 Error')
title('Species 2')
%legend
hold off


%%%% Plot of the L1 error versus time species 1

figure
%plot(error_list_s1_001(:,1),error_list_s1_001(:,3),'DisplayName','dt = 0.001');
% hold on
% plot(error_list_s1_0025(:,1),error_list_s1_0025(:,3),'DisplayName','dt = 0.0025');
% hold on
% plot(error_list_s1_0005(:,1),error_list_s1_0005(:,3),'DisplayName','dt = 0.0005');
%hold on
plot(error_list_s1_0001667(:,1),error_list_s1_0001667(:,3),'DisplayName','dt = 0.0001667');
hold on 
plot(error_list_s1_000125(:,1),error_list_s1_000125(:,3),'DisplayName','dt = 0.000125');
xlabel('Time')
ylabel('L1 Error')
title('L1 error plot Species 1')
legend
hold off



%%%% Plot of the L1 error versus time species 2

figure
%plot(error_list_s2_001(:,1),error_list_s2_001(:,3),'DisplayName','dt = 0.001');
% hold on
% plot(error_list_s2_0025(:,1),error_list_s2_0025(:,3),'DisplayName','dt = 0.0025');
% hold on
% plot(error_list_s2_0005(:,1),error_list_s2_0005(:,3),'DisplayName','dt = 0.0005');
%hold on
plot(error_list_s2_0001667(:,1),error_list_s2_0001667(:,3),'DisplayName','dt = 0.0001667');
hold on 
plot(error_list_s2_000125(:,1),error_list_s2_000125(:,3),'DisplayName','dt = 0.000125');
xlabel('Time')
ylabel('L1 Error')
title('L1 error plot Species 1')
legend
hold off



% %%%% Log Log plot of the L1, L2, and Linfty error to display order of
% %%%% convergence
% P = 4; % Number of data sets
% Linf_error_s1 = zeros(P,1);
% L1_error_s1 = zeros(P,1);
% L2_error_s1 = zeros(P,1);
% Linf_error_s2 = zeros(P,1);
% L1_error_s2 = zeros(P,1);
% L2_error_s2 = zeros(P,1);
% h_list_s1 = zeros(P,1);
% h_list_s2 = zeros(P,1);
% 
% 
% 
% % The first column is h
% h_list_s1(1) = data_n_60.dv1;
% h_list_s1(2) = data_n_80.dv1;
% h_list_s1(3) = data_n_100.dv1;
% h_list_s1(4) = data_n_120.dv1;
% % h_list_s1(5) = data_n_150.dv1;
% 
% % The first column is h
% h_list_s2(1) = data_n_60.dv2;
% h_list_s2(2) = data_n_80.dv2;
% h_list_s2(3) = data_n_100.dv2;
% h_list_s2(4) = data_n_120.dv2;
% % h_list_s2(5) = data_n_150.dv2;
% 
% 
% %the second column of error_list is the Linf error
% % Linf_error_s1(1) = error_list_s1_40(Nt_40,2);
% Linf_error_s1(1) = error_list_s1_60(Nt_60,2);
% Linf_error_s1(2) = error_list_s1_80(Nt_80,2);
% Linf_error_s1(3) = error_list_s1_100(Nt_100,2);
% Linf_error_s1(4) = error_list_s1_120(Nt_120,2);
% % Linf_error_s1(5) = error_list_s1_150(Nt_150,2);
% 
% %the third column of error_list is the L1 error
% % L1_error_s1(1) = error_list_s1_40(Nt_40,3);
% L1_error_s1(1) = error_list_s1_60(Nt_60,3);
% L1_error_s1(2) = error_list_s1_80(Nt_80,3);
% L1_error_s1(3) = error_list_s1_100(Nt_100,3);
% L1_error_s1(4) = error_list_s1_120(Nt_120,3);
% % L1_error_s1(5) = error_list_s1_150(Nt_150,3);
% 
% %the forth column of error_list is the L2 error
% 
% % L2_error_s1(1) = error_list_s1_40(Nt_40,4);
% L2_error_s1(1) = error_list_s1_60(Nt_60,4);
% L2_error_s1(2) = error_list_s1_80(Nt_80,4);
% L2_error_s1(3) = error_list_s1_100(Nt_100,4);
% L2_error_s1(4) = error_list_s1_120(Nt_120,4);
% % L2_error_s1(5) = error_list_s1_150(Nt_150,4);
% 
% 
% %the second column of error_list is the Linf error
% % Linf_error_s2(1) = error_list_s2_40(Nt_40,2);
% Linf_error_s2(1) = error_list_s2_60(Nt_60,2);
% Linf_error_s2(2) = error_list_s2_80(Nt_80,2);
% Linf_error_s2(3) = error_list_s2_100(Nt_100,2);
% Linf_error_s2(4) = error_list_s2_120(Nt_120,2);
% % Linf_error_s2(5) = error_list_s2_150(Nt_150,2);
% 
% 
% %the third column of error_list is the L1 error
% % L1_error_s2(1) = error_list_s2_40(Nt_40,3);
% L1_error_s2(1) = error_list_s2_60(Nt_60,3);
% L1_error_s2(2) = error_list_s2_80(Nt_80,3);
% L1_error_s2(3) = error_list_s2_100(Nt_100,3);
% L1_error_s2(4) = error_list_s2_120(Nt_120,3);
% % L1_error_s2(5) = error_list_s2_150(Nt_150,3);
% 
% 
% %the forth column of error_list is the L2 error
% 
% % L2_error_s2(1) = error_list_s2_40(Nt_40,4);
% L2_error_s2(1) = error_list_s2_60(Nt_60,4);
% L2_error_s2(2) = error_list_s2_80(Nt_80,4);
% L2_error_s2(3) = error_list_s2_100(Nt_100,4);
% L2_error_s2(4) = error_list_s2_120(Nt_120,4);
% % L2_error_s2(5) = error_list_s2_150(Nt_150,4);
% 
% x_s1 = ones(length(h_list_s1),2);
% x_s1(:,2) = log(h_list_s1);
% 
% x_s2 = ones(length(h_list_s2),2);
% x_s2(:,2) = log(h_list_s2);
% 
% 
% Linf_roc_s1 = x_s1\log(Linf_error_s1); %slope using least squares fitting
% L1_roc_s1 = x_s1\log(L1_error_s1); %slope using least squares fitting
% L2_roc_s1 = x_s1\log(L2_error_s1); %slope using a least squares fitting
% 
% Linf_roc_s2 = x_s2\log(Linf_error_s2); %slope using least squares fitting
% L1_roc_s2 = x_s2\log(L1_error_s2); %slope using least squares fitting
% L2_roc_s2 = x_s2\log(L2_error_s2); %slope using a least squares fitting
% 
% 
% %Rate of convergence plot species 1
% figure
% plot(-log(h_list_s1),-log(L2_error_s1),'.','markersize',8,'DisplayName','L2error')
% hold on 
% plot(-log(h_list_s1),-log(L1_error_s1),'.','markersize',8,'DisplayName','L1error')
% hold on
% plot(-log(h_list_s1),-log(Linf_error_s1),'.','markersize',8,'DisplayName','Linferror')
% hold on
% plot(-log(h_list_s1),-L2_roc_s1(2)*log(h_list_s1)-L2_roc_s1(1),'DisplayName',['L2 roc = ',num2str(L2_roc_s1(2))])
% hold on
% plot(-log(h_list_s1),-L1_roc_s1(2)*log(h_list_s1)-L1_roc_s1(1),'DisplayName',['L1 roc = ',num2str(L1_roc_s1(2))])
% hold on
% plot(-log(h_list_s1),-Linf_roc_s1(2)*log(h_list_s1)-Linf_roc_s1(1),'DisplayName',['Linf roc = ',num2str(Linf_roc_s1(2))])
% hold on
% plot(-log(h_list_s1),-2*log(h_list_s1),'DisplayName','Reference Slope = 2')
% xlabel('log(h)')
% ylabel('log(error)')
% title('Rate of convergence Species 1 m1 > m2')
% legend
% 
% %Rate of convergence plot species 2
% figure
% plot(-log(h_list_s2),-log(L2_error_s2),'.','markersize',8,'DisplayName','L2error')
% hold on 
% plot(-log(h_list_s2),-log(L1_error_s2),'.','markersize',8,'DisplayName','L1error')
% hold on
% plot(-log(h_list_s2),-log(Linf_error_s2),'.','markersize',8,'DisplayName','Linferror')
% hold on
% plot(-log(h_list_s2),-L2_roc_s2(2)*log(h_list_s2)-L2_roc_s2(1),'DisplayName',['L2 roc = ',num2str(L2_roc_s2(2))])
% hold on
% plot(-log(h_list_s2),-L1_roc_s2(2)*log(h_list_s2)-L1_roc_s2(1),'DisplayName',['L1 roc = ',num2str(L1_roc_s2(2))])
% hold on
% plot(-log(h_list_s2),-Linf_roc_s2(2)*log(h_list_s2)-Linf_roc_s2(1),'DisplayName',['Linf roc = ',num2str(Linf_roc_s2(2))])
% hold on
% plot(-log(h_list_s2),-2*log(h_list_s2),'DisplayName','Reference Slope = 2')
% xlabel('log(h)')
% ylabel('log(error)')
% title('Rate of convergence Species 2 M1 > M2')
% legend
total_moments_000125 = data_dt_000125.total_moments;

dt_00025 = data_dt_00025.dt;
dt_000125 = data_dt_000125.dt;
dt_0001 = data_dt_0001.dt;
%n_001 = data_dt_001.n; 
n_0001667 = data_dt_0001667.n; 
n_000125 = data_dt_000125.n; 
n_00025 = data_dt_00025.n;
n_0001 = data_dt_0001.n;


m = data_dt_000125.m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the initial energy for the time evolution plot

Vx1 = data_dt_0001.Vrx1; Vy1 = data_dt_0001.Vry1;
Vx2 = data_dt_0001.Vrx2; Vy2 = data_dt_0001.Vry2;

W1 = data_dt_0001.W1; W2 = data_dt_0001.W2;

E1 = m(1)*sum(W1.*(Vx1.^2 + Vy1.^2)) + m(2)*sum(W2.*(Vx2.^2+Vy2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %Plot of total energy
%figure('Renderer','painters','OuterPosition',[10,10,700,700])
figure
% options.units = 'centimeters';
% options.FontMode = 'Fixed'; %This line and next are the important bit
% options.FixedFontSize = 12;
% options.Width = 25;
% options.Height = 25;
% options.format = 'eps'; %or whatever options you'd like
% options.FontName = 'arial'; 
% options.Renderer = 'painters';
% plot(error_list_s1_000125(:,1),error_list_s1_000125(:,8)+...
%     error_list_s2_000125(:,8),'DisplayName',['$n = $',num2str(n_000125),...
%      ', $\Delta t = $',num2str(dt_000125)]);
% hold on 
% plot(error_list_s1_000125(:,1),m(1)*error_list_s1_000125(:,8)+...
%     m(2)*error_list_s2_000125(:,8),'DisplayName',['$n = $',num2str(n_000125),...
%      ', $\Delta t = $',num2str(dt_000125)]);
% hold on
plot(error_list_s1_0001(:,1),error_list_s1_0001(:,8)+...
    error_list_s2_0001(:,8)-E1,'DisplayName',[...
     '$\Delta t = $',num2str(dt_0001)]);
% hold off
xlabel('$t$','Interpreter','latex')
ylabel('$K-K_0$','Interpreter','latex')
% ytickformat('%,.2f')
title('\textbf{Total Energy}','Interpreter','latex')
% hl = legend('Show');
% set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_3_Total_Energy.pdf','-dpdf','-r0')
% hgexport(f,'Example 3 Total Energy.eps')
% exportgraphics(f,'Example 3 Total Energy.pdf')





% Plot of total entropy

%Plot of total entropy and species entropy
%figure('Renderer','painters','OuterPosition',[10,10,700,700])
figure
plot(error_list_s1_000125(:,1),error_list_s1_000125(:,9)+...
    error_list_s2_000125(:,9),'DisplayName','Total Entropy');
% hold on
% plot(error_list_s1_000125(:,1),error_list_s1_000125(:,9),'DisplayName','Species 1 Entropy');
% hold on 
% plot(error_list_s2_000125(:,1),error_list_s2_000125(:,9),'DisplayName','Species 2 Entropy 2');
%hold off
xlabel('$t$','Interpreter','latex')
ylabel('$E^N$','Interpreter','latex')
title('\textbf{Total Entropy}','Interpreter','latex')
% hl = legend('show');
% set(hl,'Interpreter','latex')
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_3_Total_Entropy.pdf','-dpdf','-r0')
% exportgraphics(f,'Example 3 Total Entropy.pdf')
% hgexport(f,'Example 3 Total Entropy.eps')
%print(f,'Example 3 Total Entropy.eps')


%Total momentum plots
%x-direction
figure
plot(error_list_s1_0001(:,1),error_list_s1_0001(:,6)...
    +error_list_s2_0001(:,6))
hold on 
% plot(error_list_s1_00025(:,1),error_list_s1_00025(:,6)...
%     +error_list_s2_00025(:,6))
% hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
title('\textbf{Total Momentum in $v_1$ Direction}','Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 2 Total Momentum v1 direction.eps')

%y-direction
figure
plot(error_list_s1_0001(:,1),error_list_s1_0001(:,7)...
    +error_list_s2_0001(:,7))
hold on 
% plot(error_list_s1_000125(:,1),error_list_s1_000125(:,7)...
%     +error_list_s2_000125(:,7))
% hold off
xlabel('$t$','Interpreter','latex')
ylabel('$P$','Interpreter','latex')
title('\textbf{Total Momentum in $v_2$ Direction}','Interpreter','latex')
f = gcf;
% set(f,'Units','pixels')
% set(f,'Position',[440 378 560 420])
exportgraphics(f,'Example 2 Total Momentum v2 direction.eps')


%Plot of numerical and exact solutions at final time

%vr1_dt_001 = data_dt_001.vr1; vr2_dt_001 = data_dt_001.vr2;

vr1_dt_0001667 = data_dt_0001667.vr1; vr2_dt_0001667 = data_dt_0001667.vr2;
vr1_dt_000125 = data_dt_000125.vr1; vr2_dt_000125 = data_dt_000125.vr2;
vr1_dt_0001 = data_dt_0001.vr1; vr2_dt_0001 = data_dt_0001.vr2;

%Nr_001 = data_dt_001.Nr; 
Nr_0001667 = data_dt_0001667.Nr; 
Nr_000125 = data_dt_000125.Nr; 
Nr_0001 = data_dt_0001.Nr;



%f1_001 = data_dt_001.f1; f2_001 = data_dt_001.f2;
f1_0001667 = data_dt_0001667.f1; f2_0001667 = data_dt_0001667.f2;
f1_000125 = data_dt_000125.f1; f2_000125 = data_dt_000125.f2;
f1_0001 = data_dt_0001.f1; f2_0001 = data_dt_0001.f2;

f1_exact_0001667 = data_dt_0001667.f1_exact; f2_exact_0001667 = data_dt_0001667.f2_exact;
f1_exact_000125 = data_dt_000125.f1_exact; f2_exact_000125 = data_dt_000125.f2_exact;
f1_exact_0001 = data_dt_0001.f1_exact; f2_exact_0001 = data_dt_0001.f2_exact;
figure
% options.units = 'centimeters';
% options.FontMode = 'Fixed'; %This line and next are the important bit
% options.FixedFontSize = 12;
% options.Width = 25;
% options.Height = 20;
% options.format = 'eps'; %or whatever options you'd like
% options.FontName = 'arial'; 
% options.Renderer = 'painters';
% plot(vr1_dt_001,f1_001(:,Nr_001/2),'DisplayName','Particle Solution dt = 0.001')
% hold on
plot(vr1_dt_0001,f1_000125(:,Nr_0001/2),'-o',...
    'MarkerIndices',round(linspace(1,length(vr1_dt_0001),30)),...
     'DisplayName','Particle Solution')
hold on
plot(vr1_dt_0001,f1_exact_000125(:,Nr_0001/2),'DisplayName','Exact Solution')
hold off
title('$f_1(:,n/2)$','Interpreter','latex')
hl = legend('Show');
set(hl,'Interpreter','latex')
ylim([0.0,18])
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_3_species_1.pdf','-dpdf','-r0')
%hgexport(f,'Example 3 species 1 pdf.eps',options)
% exportgraphics(f,'Example 3 species 1 pdf.eps')


figure
% options.units = 'centimeters';
% options.FontMode = 'Fixed'; %This line and next are the important bit
% options.FixedFontSize = 12;
% options.Width = 25;
% options.Height = 20;
% options.format = 'eps'; %or whatever options you'd like
% options.FontName = 'arial'; 
% options.Renderer = 'painters';
%plot(vr2_dt_001,f2_001(:,Nr_001/2),'DisplayName','Particle Solution dt = 0.001')
%hold on
plot(vr2_dt_0001,f2_0001(:,Nr_000125/2),'-o',...
    'MarkerIndices',round(linspace(1,length(vr2_dt_000125),30)),...
    'DisplayName','Particle Solution')
hold on
plot(vr2_dt_0001,f2_exact_0001(:,Nr_0001/2),'DisplayName','Exact Solution')
hold off
title('$f_2(:,n/2)$','Interpreter','latex')
hl = legend('Show');
set(hl,'Interpreter','latex')
ylim([0.0,0.18])
f = gcf;
fontsize(f,20,'points')
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Example_3_species_2.pdf','-dpdf','-r0')
%exportgraphics(f,'Example 3 species 2 pdf.eps')