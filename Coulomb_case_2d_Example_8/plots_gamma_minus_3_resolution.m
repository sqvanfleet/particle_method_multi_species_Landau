%gamma = -3 plots

clc
clear all
close all

data_n_60 = load("multi_species_particle_2d_Coulomb_n_60_Example_8_gamma_-3_dv1_0.13333_dv2_0.13333_dt_0.005tmax30.mat");
data_n_80 = load("multi_species_particle_2d_Coulomb_n_80_Example_8_gamma_-3_dv1_0.1_dv2_0.1_dt_0.005tmax30.mat");

error_list_s1_n_60 = data_n_60.error_list_s1;
error_list_s2_n_60 = data_n_60.error_list_s2;

error_list_s1_n_80 = data_n_80.error_list_s1;
error_list_s2_n_80 = data_n_80.error_list_s2;



n_60 = data_n_60.n;
n_80 = data_n_80.n;


dt_005 = data_n_60.dt;
dt_01 = data_n_80.dt;

% % Species 1 energy 
% figure 
% plot(error_list_s1_n_40(:,1),error_list_s1_n_40(:,6));
% title('Energy Species 1')
% 
% % Species 2 energy 
% figure 
% plot(error_list_s2_n_40(:,1),error_list_s2_n_40(:,6));
% title('Energy Species 2')

% Total Energy
figure
plot(error_list_s1_n_60(:,1),error_list_s1_n_60(:,6)+error_list_s2_n_60(:,6),'DisplayName',['Total Energy n = ',num2str(n_60)])
hold on
plot(error_list_s1_n_80(:,1),error_list_s1_n_80(:,6)+error_list_s2_n_80(:,6),'DisplayName',['Total Energy n = ',num2str(n_80)])
legend
title('Total Energy')
hold off

T0 = data_n_60.T_relax;
T0vec_n_60 = zeros(length(error_list_s2_n_60(:,1)),1);
T0vec_n_60 = T0vec_n_60 + T0;

T0 = data_n_80.T_relax;
T0vec_n_80 = zeros(length(error_list_s2_n_80(:,1)),1);
T0vec_n_80 = T0vec_n_80 + T0;



% Plots of temperature n = 60
figure 
plot(error_list_s1_n_60(:,1),error_list_s1_n_60(:,7),'DisplayName',['Species 1 Temp n = ',num2str(n_60)]) 
hold on
plot(error_list_s2_n_60(:,1),error_list_s2_n_60(:,7),'Displayname',['Species 2 Temp n = ',num2str(n_60)])
hold on
plot(error_list_s1_n_60(:,1),T0vec_n_60,'DisplayName','Relaxation Temperature')
xlabel('Time')
ylabel('Temperature')
title('Temperature Relaxation')
legend('Location','southeast')
hold off

% Plots of temperature n = 80
figure 
plot(error_list_s1_n_80(:,1),error_list_s1_n_80(:,7),'DisplayName',['Species 1 Temp n = ',num2str(n_80)]) 
hold on
plot(error_list_s2_n_80(:,1),error_list_s2_n_80(:,7),'Displayname',['Species 2 Temp n = ',num2str(n_80)])
hold on
plot(error_list_s1_n_80(:,1),T0vec_n_80,'DisplayName','Relaxation Temperature')
xlabel('Time')
ylabel('Temperature')
title('Temperature Relaxation')
legend('Location','southeast')
hold off



% Plots of entropy and total entropy 

figure 
plot(error_list_s1_n_60(:,1),error_list_s1_n_60(:,8),'DisplayName',['Species 1 Entropy',num2str(n_60)])
hold on 
plot(error_list_s2_n_60(:,1),error_list_s2_n_60(:,8),'DisplayName',['Species 2 Entropy',num2str(n_60)])
hold on  
plot(error_list_s1_n_60(:,1),error_list_s1_n_60(:,8)+error_list_s2_n_60(:,8),'DisplayName',['Total Entropy',num2str(n_60)])
xlabel('Time')
ylabel('Entropy')
title('Time Evolution of Entropy')
legend
hold off

figure 
plot(error_list_s1_n_80(:,1),error_list_s1_n_80(:,8),'DisplayName',['Species 1 Entropy',num2str(n_80)])
hold on 
plot(error_list_s2_n_80(:,1),error_list_s2_n_80(:,8),'DisplayName',['Species 2 Entropy',num2str(n_80)])
hold on  
plot(error_list_s1_n_80(:,1),error_list_s1_n_80(:,8)+error_list_s2_n_80(:,8),'DisplayName',['Total Entropy',num2str(n_80)])
xlabel('Time')
ylabel('Entropy')
title('Time Evolution of Entropy')
legend
hold off


%Solution at time Final time 

vr1_60 = data_n_60.vr1;
vr2_60 = data_n_60.vr2;

vr1_80 = data_n_80.vr1;
vr2_80 = data_n_80.vr2;

dv1_60 = data_n_60.dv1;
dv2_60 = data_n_60.dv2;

Nr_60 = data_n_60.Nr;
f1_60 = data_n_60.f1;
f2_60 = data_n_60.f2;

f1_0_60 = data_n_60.f1_0;
f2_0_60 = data_n_60.f2_0;

Nt_60 = data_n_60.Nt;

dv1_80 = data_n_80.dv1;
dv2_80 = data_n_80.dv2;

Nr_80 = data_n_80.Nr;
f1_80 = data_n_80.f1;
f2_80 = data_n_80.f2;

f1_0_80 = data_n_80.f1_0;
f2_0_80 = data_n_80.f2_0;

Nt_80 = data_n_80.Nt;

%species 1 x
figure
plot(vr1_60,f1_60(:,Nr_60/2),'DisplayName',['f1(:,Nr/2) at t = ',num2str(error_list_s1_n_60(Nt_60,1))]) 
hold on 
plot(vr1_60,f1_0_60(:,Nr_60/2), 'DisplayName','Initial Condition')
legend
title(['f1(:,Nr/2) n = ',num2str(n_60)])

figure
plot(vr1_80,f1_80(:,Nr_80/2),'DisplayName',['f1(:,Nr/2) at t = ',num2str(error_list_s1_n_80(Nt_80,1))]) 
hold on 
plot(vr1_80,f1_0_80(:,Nr_80/2), 'DisplayName','Initial Condition')
legend
title(['f1(:,Nr/2) n = ',num2str(n_80)])

%species 1 y
figure
plot(vr1_60,f1_60(Nr_60/2,:),'DisplayName',['f1(Nr/2,:) at t = ',num2str(error_list_s1_n_60(Nt_60,1))]) 
hold on 
plot(vr1_60,f1_0_60(Nr_60/2,:), 'DisplayName','Initial Condition')
legend
title(['f1(Nr/2,:) n = ',num2str(n_60)])

figure
plot(vr1_80,f1_80(Nr_80/2,:),'DisplayName',['f1(Nr/2,:) at t = ',num2str(error_list_s1_n_80(Nt_80,1))]) 
hold on 
plot(vr1_80,f1_0_80(Nr_80/2,:), 'DisplayName','Initial Condition')
legend
title(['f1(Nr/2,:) n = ',num2str(n_80)])

%species 2 x
figure
plot(vr2_60,f2_60(:,Nr_60/2),'DisplayName',['f2(:,Nr/2) at t = ',num2str(error_list_s2_n_60(Nt_60,1))]) 
hold on 
plot(vr2_60,f2_0_60(:,Nr_60/2), 'DisplayName','Initial Condition')
legend
title(['f2(:,Nr/2) n = ',num2str(n_60)])

figure
plot(vr2_80,f2_80(:,Nr_80/2),'DisplayName',['f2(:,Nr/2) at t = ',num2str(error_list_s2_n_80(Nt_80,1))]) 
hold on 
plot(vr2_80,f2_0_80(:,Nr_80/2), 'DisplayName','Initial Condition')
legend
title(['f2(:,Nr/2) n = ',num2str(n_80)])


%species 2 y
figure
plot(vr2_60,f2_60(Nr_60/2,:),'DisplayName',['f2(Nr/2,:) at t = ',num2str(error_list_s2_n_60(Nt_60,1))]) 
hold on 
plot(vr2_60,f2_0_60(Nr_60/2,:), 'DisplayName','Initial Condition')
legend
title(['f2(Nr/2,:) n = ',num2str(n_60)])

figure
plot(vr2_80,f2_80(Nr_80/2,:),'DisplayName',['f2(Nr/2,:) at t = ',num2str(error_list_s2_n_80(Nt_80,1))]) 
hold on 
plot(vr2_80,f2_0_80(Nr_80/2,:), 'DisplayName','Initial Condition')
legend
title(['f2(Nr/2,:) n = ',num2str(n_80)])
