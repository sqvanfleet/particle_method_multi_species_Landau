clc 
clear all
close all


myVars = {'error_list_s1','error_list_s2','dv1','dv2','Nt','total_moments','f1','f2','f1_exact','f2_exact','vr1','vr2','Nr','n','m'};



data_n_60 = load('multi_species_particle_2d_symmetric_parallel_n_20_Example_6_dv1_0.3_dv2_0.4_dt_0.001tmax0.5.mat');
data_n_80 = load('multi_species_particle_2d_symmetric_n_80_Example_6_dv1_0.075_dv2_0.1_dt_0.01tmax5.mat',myVars{:});
data_n_100 = load('multi_species_particle_2d_symmetric_n_100_Example_6_dv1_0.06_dv2_0.08_dt_0.01tmax5.mat',myVars{:});
data_n_120 = load('multi_species_particle_2d_symmetric_n_120_Example_6_dv1_0.05_dv2_0.066667_dt_0.01tmax5.mat',myVars{:});
% data_n_150 = load('multi_species_particle_2d_symmetric_n_150_Example_6_dv1_0.04_dv2_0.053333_dt_0.01tmax5.mat',myVars{:});


Nt_60 = data_n_60.Nt;
Nt_80 = data_n_80.Nt;
Nt_100 = data_n_100.Nt;
Nt_120 = data_n_120.Nt;
% Nt_150 = data_n_150.Nt;



error_list_s1_60 = data_n_60.error_list_s1;
error_list_s1_80 = data_n_80.error_list_s1;
error_list_s1_100 = data_n_100.error_list_s1;
error_list_s1_120 = data_n_120.error_list_s1;
% error_list_s1_150 = data_n_150.error_list_s1;


error_list_s2_60 = data_n_60.error_list_s2;
error_list_s2_80 = data_n_80.error_list_s2;
error_list_s2_100 = data_n_100.error_list_s2;
error_list_s2_120 = data_n_120.error_list_s2;
% error_list_s2_150 = data_n_150.error_list_s2;


%total_moments_60 = data_n_60.total_moments;
%total_moments_80 = data_n_80.total_moments;
%total_moments_100 = data_n_100.total_moments;
%total_moments_120 = data_n_120.total_moments;
% total_moments_150 = data_n_150.total_moments;

% column 1 of error list is a list of times
% column 4 of error list is a list of L2 error


%%%% Plot of the L2 error versus time species 1

figure
plot(error_list_s1_60(:,1),error_list_s1_60(:,4),'DisplayName','n = 60');
hold on
plot(error_list_s1_80(:,1),error_list_s1_80(:,4),'DisplayName','n = 80');
hold on 
plot(error_list_s1_100(:,1),error_list_s1_100(:,4),'DisplayName','n = 100');
hold on 
plot(error_list_s1_120(:,1),error_list_s1_120(:,4),'DisplayName','n = 120');
% hold on 
% plot(error_list_s1_150(:,1),error_list_s1_150(:,4),'DisplayName','n = 150');
xlabel('Time')
ylabel('L2 Error')
title('L2 error plot Species 1')
legend
hold off



%%%% Plot of the L2 error versus time species 2

figure
plot(error_list_s2_60(:,1),error_list_s2_60(:,4),'DisplayName','n = 60');
hold on
plot(error_list_s2_80(:,1),error_list_s2_80(:,4),'DisplayName','n = 80');
hold on 
plot(error_list_s2_100(:,1),error_list_s2_100(:,4),'DisplayName','n = 100');
hold on 
plot(error_list_s2_120(:,1),error_list_s2_120(:,4),'DisplayName','n = 120');
% hold on
% plot(error_list_s2_150(:,1),error_list_s2_150(:,4),'DisplayName','n = 150');
xlabel('Time')
xlabel('Time')
ylabel('L2 Error')
title('L2 error plot Species 2')
legend
hold off

%%%% Plot of the L1 error versus time species 1

figure
plot(error_list_s1_60(:,1),error_list_s1_60(:,3),'DisplayName','n = 60');
hold on
plot(error_list_s1_80(:,1),error_list_s1_80(:,3),'DisplayName','n = 80');
hold on 
plot(error_list_s1_100(:,1),error_list_s1_100(:,3),'DisplayName','n = 100');
hold on 
plot(error_list_s1_120(:,1),error_list_s1_120(:,3),'DisplayName','n = 120');
% hold on 
% plot(error_list_s1_150(:,1),error_list_s1_150(:,3),'DisplayName','n = 150');
xlabel('Time')
ylabel('L1 Error')
title('L1 error plot Species 1')
legend
hold off



%%%% Plot of the L1 error versus time species 2

figure
plot(error_list_s2_60(:,1),error_list_s2_60(:,3),'DisplayName','n = 60');
hold on
plot(error_list_s2_80(:,1),error_list_s2_80(:,3),'DisplayName','n = 80');
hold on 
plot(error_list_s2_100(:,1),error_list_s2_100(:,3),'DisplayName','n = 100');
hold on 
plot(error_list_s2_120(:,1),error_list_s2_120(:,3),'DisplayName','n = 120');
% hold on
% plot(error_list_s2_150(:,1),error_list_s2_150(:,3),'DisplayName','n = 150');
xlabel('Time')
xlabel('Time')
ylabel('L1 Error')
title('L1 error plot Species 2')
legend
hold off




%%%% Log Log plot of the L1, L2, and Linfty error to display order of
%%%% convergence
P = 4; % Number of data sets
Linf_error_s1 = zeros(P,1);
L1_error_s1 = zeros(P,1);
L2_error_s1 = zeros(P,1);
Linf_error_s2 = zeros(P,1);
L1_error_s2 = zeros(P,1);
L2_error_s2 = zeros(P,1);
h_list_s1 = zeros(P,1);
h_list_s2 = zeros(P,1);



% The first column is h
h_list_s1(1) = data_n_60.dv1;
h_list_s1(2) = data_n_80.dv1;
h_list_s1(3) = data_n_100.dv1;
h_list_s1(4) = data_n_120.dv1;
% h_list_s1(5) = data_n_150.dv1;

% The first column is h
h_list_s2(1) = data_n_60.dv2;
h_list_s2(2) = data_n_80.dv2;
h_list_s2(3) = data_n_100.dv2;
h_list_s2(4) = data_n_120.dv2;
% h_list_s2(5) = data_n_150.dv2;


%the second column of error_list is the Linf error
% Linf_error_s1(1) = error_list_s1_40(Nt_40,2);
Linf_error_s1(1) = error_list_s1_60(Nt_60,2);
Linf_error_s1(2) = error_list_s1_80(Nt_80,2);
Linf_error_s1(3) = error_list_s1_100(Nt_100,2);
Linf_error_s1(4) = error_list_s1_120(Nt_120,2);
% Linf_error_s1(5) = error_list_s1_150(Nt_150,2);

%the third column of error_list is the L1 error
% L1_error_s1(1) = error_list_s1_40(Nt_40,3);
L1_error_s1(1) = error_list_s1_60(Nt_60,3);
L1_error_s1(2) = error_list_s1_80(Nt_80,3);
L1_error_s1(3) = error_list_s1_100(Nt_100,3);
L1_error_s1(4) = error_list_s1_120(Nt_120,3);
% L1_error_s1(5) = error_list_s1_150(Nt_150,3);

%the forth column of error_list is the L2 error

% L2_error_s1(1) = error_list_s1_40(Nt_40,4);
L2_error_s1(1) = error_list_s1_60(Nt_60,4);
L2_error_s1(2) = error_list_s1_80(Nt_80,4);
L2_error_s1(3) = error_list_s1_100(Nt_100,4);
L2_error_s1(4) = error_list_s1_120(Nt_120,4);
% L2_error_s1(5) = error_list_s1_150(Nt_150,4);


%the second column of error_list is the Linf error
% Linf_error_s2(1) = error_list_s2_40(Nt_40,2);
Linf_error_s2(1) = error_list_s2_60(Nt_60,2);
Linf_error_s2(2) = error_list_s2_80(Nt_80,2);
Linf_error_s2(3) = error_list_s2_100(Nt_100,2);
Linf_error_s2(4) = error_list_s2_120(Nt_120,2);
% Linf_error_s2(5) = error_list_s2_150(Nt_150,2);


%the third column of error_list is the L1 error
% L1_error_s2(1) = error_list_s2_40(Nt_40,3);
L1_error_s2(1) = error_list_s2_60(Nt_60,3);
L1_error_s2(2) = error_list_s2_80(Nt_80,3);
L1_error_s2(3) = error_list_s2_100(Nt_100,3);
L1_error_s2(4) = error_list_s2_120(Nt_120,3);
% L1_error_s2(5) = error_list_s2_150(Nt_150,3);


%the forth column of error_list is the L2 error

% L2_error_s2(1) = error_list_s2_40(Nt_40,4);
L2_error_s2(1) = error_list_s2_60(Nt_60,4);
L2_error_s2(2) = error_list_s2_80(Nt_80,4);
L2_error_s2(3) = error_list_s2_100(Nt_100,4);
L2_error_s2(4) = error_list_s2_120(Nt_120,4);
% L2_error_s2(5) = error_list_s2_150(Nt_150,4);

x_s1 = ones(length(h_list_s1),2);
x_s1(:,2) = log(h_list_s1);

x_s2 = ones(length(h_list_s2),2);
x_s2(:,2) = log(h_list_s2);


Linf_roc_s1 = x_s1\log(Linf_error_s1); %slope using least squares fitting
L1_roc_s1 = x_s1\log(L1_error_s1); %slope using least squares fitting
L2_roc_s1 = x_s1\log(L2_error_s1); %slope using a least squares fitting

Linf_roc_s2 = x_s2\log(Linf_error_s2); %slope using least squares fitting
L1_roc_s2 = x_s2\log(L1_error_s2); %slope using least squares fitting
L2_roc_s2 = x_s2\log(L2_error_s2); %slope using a least squares fitting


%Rate of convergence plot species 1
figure
plot(-log(h_list_s1),-log(L2_error_s1),'.','markersize',8,'DisplayName','L2error')
hold on 
plot(-log(h_list_s1),-log(L1_error_s1),'.','markersize',8,'DisplayName','L1error')
hold on
plot(-log(h_list_s1),-log(Linf_error_s1),'.','markersize',8,'DisplayName','Linferror')
hold on
plot(-log(h_list_s1),-L2_roc_s1(2)*log(h_list_s1)-L2_roc_s1(1),'DisplayName',['L2 roc = ',num2str(L2_roc_s1(2))])
hold on
plot(-log(h_list_s1),-L1_roc_s1(2)*log(h_list_s1)-L1_roc_s1(1),'DisplayName',['L1 roc = ',num2str(L1_roc_s1(2))])
hold on
plot(-log(h_list_s1),-Linf_roc_s1(2)*log(h_list_s1)-Linf_roc_s1(1),'DisplayName',['Linf roc = ',num2str(Linf_roc_s1(2))])
hold on
plot(-log(h_list_s1),-2*log(h_list_s1),'DisplayName','Reference Slope = 2')
xlabel('log(h)')
ylabel('log(error)')
title('Rate of convergence Species 1 m1 > m2')
legend

%Rate of convergence plot species 2
figure
plot(-log(h_list_s2),-log(L2_error_s2),'.','markersize',8,'DisplayName','L2error')
hold on 
plot(-log(h_list_s2),-log(L1_error_s2),'.','markersize',8,'DisplayName','L1error')
hold on
plot(-log(h_list_s2),-log(Linf_error_s2),'.','markersize',8,'DisplayName','Linferror')
hold on
plot(-log(h_list_s2),-L2_roc_s2(2)*log(h_list_s2)-L2_roc_s2(1),'DisplayName',['L2 roc = ',num2str(L2_roc_s2(2))])
hold on
plot(-log(h_list_s2),-L1_roc_s2(2)*log(h_list_s2)-L1_roc_s2(1),'DisplayName',['L1 roc = ',num2str(L1_roc_s2(2))])
hold on
plot(-log(h_list_s2),-Linf_roc_s2(2)*log(h_list_s2)-Linf_roc_s2(1),'DisplayName',['Linf roc = ',num2str(Linf_roc_s2(2))])
hold on
plot(-log(h_list_s2),-2*log(h_list_s2),'DisplayName','Reference Slope = 2')
xlabel('log(h)')
ylabel('log(error)')
title('Rate of convergence Species 2 M1 > M2')
legend

m = [2,1];

%Plot of total energy
figure
plot(error_list_s1_120(:,1),m(1)*error_list_s1_120(:,8)+m(2)*error_list_s2_120(:,8),'DisplayName','n = 120');
% ylim([2.9,3.1])
xlabel('Time')
ylabel('Energy')
title('Total Energy Plot')
legend
hold off

%Test plot of total energy
% figure
% plot(error_list_s1_120(:,1),error_list_s2_120(:,8)+m(1)*error_list_s1_120(:,8))
% xlabel('Time')
% ylabel('Energy')
% title('Total Energy Plot 2')
% legend
% hold off
% 
%Plot of total entropy and species entropy
figure
%plot(error_list_s1_120(:,1),total_moments_120(:,5),'DisplayName','Total Entropy');
%hold on
plot(error_list_s1_120(:,1),error_list_s1_120(:,9),'DisplayName','Species 1 Entrpy');
hold on 
plot(error_list_s2_120(:,1),error_list_s2_120(:,9),'DisplayName','Species 2 Entrpy 2');
% ylim([2.9,3.1])
xlabel('Time')
ylabel('Entropy')
title('Entropy plot n = 120')
legend
hold off



%Plot of numerical and exact solutions at final time

vr1_60 = data_n_60.vr1; vr2_60 = data_n_60.vr2;
vr1_80 = data_n_80.vr1; vr2_80 = data_n_80.vr2;
vr1_100 = data_n_100.vr1; vr2_100 = data_n_100.vr2;
vr1_120 = data_n_120.vr1; vr2_120 = data_n_120.vr2;
% vr1_150 = data_n_150.vr1; vr2_150 = data_n_150.vr2;

Nr_60 = data_n_60.Nr; Nr_80 = data_n_80.Nr; 
Nr_100 = data_n_100.Nr;
Nr_120 = data_n_120.Nr; 
% Nr_150 = data_n_150.Nr; 

n_60 = data_n_60.n; n_80 = data_n_80.n; 
n_100 = data_n_100.n;
n_120 = data_n_120.n; 
% n_150 = data_n_150.n; 

f1_60 = data_n_60.f1; f2_60 = data_n_60.f2;
f1_80 = data_n_80.f1; f2_80 = data_n_80.f2;
f1_100 = data_n_100.f1; f2_100 = data_n_100.f2;
f1_120 = data_n_120.f1; f2_120 = data_n_120.f2;
% f1_150 = data_n_150.f1; f2_150 = data_n_150.f2;

f1_exact_120 = data_n_120.f1_exact; f2_exact_120 = data_n_120.f2_exact;

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


