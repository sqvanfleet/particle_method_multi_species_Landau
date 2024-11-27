clc
clear all
close all

% 2D Maxwell molecule 
gamma = 0;

% Example 1
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_1_mass_ratio_2();

% Example 2 
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_2_mass_ratio_2();

% Example 3
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_3();

% Example 4
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_4();

% Example 5
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_5();

% Example 6
[Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_6();

% Example 7
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_7();

% Example 8
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_8();

% Example 9
%[Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_9();

% Example 10
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_10();

% Example 11
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_11();

% Example 12
%[Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_12();

% Example 13
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_13();

% Example 14
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_14();

% Example 15
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_15();

% Example 16
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_16();

% Example 17
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_17();

n = 55;
% total number of particles
Np = n^2; 

% initial mesh species 1
dv1 = 2*Vmax1/n;
v1 = (-Vmax1+dv1/2):dv1:(Vmax1-dv1/2); %List of cell centers species 1
[vx1,vy1] = ndgrid(v1); % species 1 initial velocities
Vx1 = vx1(:); Vy1 = vy1(:); %species 1 initial velocities

dv2 = 2*Vmax2/n;
v2 = (-Vmax2+dv2/2):dv2:(Vmax2-dv2/2); %cell centers species 2
[vx2,vy2] = ndgrid(v2); % species 2 initial velocities
Vx2 = vx2(:); Vy2 = vy2(:); %species 2 initial velocities

% initial weights
t0 = 0; 
[f1_0,f2_0] = multi_species_exact_2d(t0,vx1,vy1,vx2,vy2,beta,m,n1);
w1 = dv1^2*f1_0; w2 = dv2^2*f2_0;
W1 = w1(:); % Weights for species one
W2 = w2(:); % Weights for species two

% reconstruction mesh species 1
Nr = n;
dvr1 = 2*Vmax1/Nr;
vr1 = (-Vmax1+dvr1/2):dvr1:(Vmax1-dvr1/2); 
[vrx1,vry1] = ndgrid(vr1);
Vrx1 = vrx1(:);
Vry1 = vry1(:);

% reconstruction mesh species 2
Nr = n;
dvr2 = 2*Vmax2/Nr;
vr2 = (-Vmax2+dvr2/2):dvr2:(Vmax2-dvr2/2); 
[vrx2,vry2] = ndgrid(vr2);
Vrx2 = vrx2(:);
Vry2 = vry2(:);


% reconstructed solution
f1 = zeros(Nr,Nr); f2 = zeros(Nr,Nr);
% choosing epsilon
epsilon1 = 4*(0.4*(dv1)^0.99)^2;
epsilon2 = 4*(0.4*(dv2)^0.99)^2;
parfor i = 1:Nr
    for j = 1:Nr
        f1(i,j) = sum(W1.*psi_2d(vrx1(i,j)-Vx1,vry1(i,j)-Vy1,epsilon1));
        f2(i,j) = sum(W2.*psi_2d(vrx2(i,j)-Vx2,vry2(i,j)-Vy2,epsilon2));
    end
end

% initial comparison
[f1_exact,f2_exact] = multi_species_exact_2d(t0,vrx1,vry1,vrx2,vry2,beta,m,n1);
% figure
% plot(vr1,f1_exact(:,Nr/2),'DisplayName','Exact Solution') 
% hold on
% plot(vr1,f1(:,Nr/2),'DisplayName','Particle Method')
% legend
% title('Species 1')
% 
% figure
% plot(vr2,f2_exact(:,Nr/2),'DisplayName','Exact Solution')
% hold on 
% plot(vr2,f2(:,Nr/2),'DisplayName','Particle Method')
% legend
% title('Species 2')

tmax = 5;
dt = 0.01/8;
%Nt = 1;
Nt = round((tmax-t0)/dt);
error_list_s1 = zeros(Nt,9);
error_list_s2 = zeros(Nt,9);


max_iter = 300;
tol = 1e-8;
Fixed_Point_list = zeros(Nt,2);

tic %total time

for nt = 1:Nt
    time = t0+dt*nt;
    %Initial Guess forward Euler step 
    [Ux1,Uy1,Ux2,Uy2] = right_hand_side_implicit_mid_multi_species_parallel(W1,W2,...
        Vx1,Vx2,Vy1,Vy2,Vx1,Vx2,Vy1,Vy2,Vrx1,Vrx2,Vry1,Vry2,...
        dv1,dv2,B1,m,Np,epsilon1,epsilon2);

    Vx1_new = Vx1 + dt*Ux1;
    Vy1_new = Vy1 + dt*Uy1;
    Vx2_new = Vx2 + dt*Ux2;
    Vy2_new = Vy2 + dt*Uy2;

    %Fixed point iteration to solve non-linear equation from implicit
    %midpoint

    for i = 1:max_iter
        Vx1_old = Vx1_new; Vx2_old = Vx2_new;
        Vy1_old = Vy1_new; Vy2_old = Vy2_new;

        [Ux1,Uy1,Ux2,Uy2] = ...
            right_hand_side_implicit_mid_multi_species_parallel(W1,W2,Vx1,Vx2,...
            Vy1,Vy2,Vx1_new,Vx2_new,Vy1_new,Vy2_new,Vrx1,Vrx2,...
            Vry1,Vry2,dv1,dv2,B1,m,Np,epsilon1,epsilon2);
        
        Vx1_new = Vx1 + dt*Ux1;
        Vx2_new = Vx2 + dt*Ux2;
        Vy1_new = Vy1 + dt*Uy1;
        Vy2_new = Vy2 + dt*Uy2;

        absres = sqrt(sum((Vx1_new - Vx1_old).^2 +(Vx2_new - Vx2_old).^2 + ...
            (Vy1_new - Vy1_old).^2 + (Vy2_new - Vy2_old).^2));

	    relres = absres/(sqrt(sum(Vx1_new.^2 + Vy1_new.^2 + Vx2_new.^2 + Vy2_new.^2)));

        disp(['res = ',num2str(relres)]);

        if relres < tol
            disp(['fixed point iteration took ',num2str(i),' iterations at time ', num2str(time)])
            break
        end
        if i == max_iter
            Fixed_Point_list(nt,1) = time;
            Fixed_Point_list(nt,2) = relres;
            disp(['maximum number of iterations reached at time ',num2str(time)])
        end


    end
    
    Vx1 = Vx1_new; Vx2 = Vx2_new;
    Vy1 = Vy1_new; Vy2 = Vy2_new;

    
    %reconstruction
    parfor i = 1:Nr
        for j = 1:Nr
            f1(i,j) = sum(W1.*psi_2d(vrx1(i,j)-Vx1,vry1(i,j)-Vy1,epsilon1));
            f2(i,j) = sum(W2.*psi_2d(vrx2(i,j)-Vx2,vry2(i,j)-Vy2,epsilon2));
        end
    end


    % plot

    [f1_exact,f2_exact] = multi_species_exact_2d(time,vrx1,vry1,vrx2,vry2,beta,m,n1);
%     figure(1)
%     plot(vr1,f1_exact(:,Nr/2),'DisplayName','Exact Solution') 
%     hold on
%     plot(vr1,f1(:,Nr/2),'DisplayName','Particle Method')
%     title('Species 1')
%     legend
%     hold off
%     drawnow
%     
%     figure(2)
%     plot(vr2,f2_exact(:,Nr/2),'DisplayName','Exact Solution')
%     hold on 
%     plot(vr2,f2(:,Nr/2),'DisplayName','Particle Method')
%     title('Species 2')
%     legend
%     hold off
%     drawnow

    disp("current time: ");
    disp(time);

    %save the moments and errors at each time step
    error_list_s1(nt,1) = time;
    error_list_s2(nt,1) = time;

    f1_error = f1-f1_exact;
    f2_error = f2-f2_exact;

    %relative error for species 1
    Linf_error = max(max(abs(f1_error)))/max(max(abs(f1_exact)));
    L1_error = sum(sum(abs(f1_error)))/sum(sum(abs(f1_exact)));
    L2_error = sqrt(sum(sum(f1_error.^2))/sum(sum(f1_exact.^2))); 
    error_list_s1(nt,2) = Linf_error;
    error_list_s1(nt,3) = L1_error;
    error_list_s1(nt,4) = L2_error;

    % moments for species 1

    rho = sum(W1); %mass
    m1 = m(1)*sum(W1.*Vx1); %momentum in Vx component   
    m2 = m(1)*sum(W1.*Vy1); %momentum in Vy component
    E = m(1)*sum(W1.*(Vx1.^2+Vy1.^2)); %energy
    inside = zeros(Nr^2,1);
    parfor i = 1:Nr^2
        inside(i) = sum(W1.*psi_2d(Vrx1(i) - Vx1,Vry1(i) - Vy1,epsilon1));
    end

    eta = dv1^2*sum(inside.*log(inside)); %entropy
    error_list_s1(nt,5) = rho;
    error_list_s1(nt,6) = m1;
    error_list_s1(nt,7) = m2;
    error_list_s1(nt,8) = E;
    error_list_s1(nt,9) = eta;

    %relative error for species 2

    Linf_error = max(max(abs(f2_error)))/max(max(abs(f2_exact)));
    L1_error = sum(sum(abs(f2_error)))/sum(sum(abs(f2_exact)));
    L2_error = sqrt(sum(sum(f2_error.^2))/sum(sum(f2_exact.^2))); 
    error_list_s2(nt,2) = Linf_error;
    error_list_s2(nt,3) = L1_error;
    error_list_s2(nt,4) = L2_error;

    % moments for species 2

    rho = m(2)*sum(W2); %mass
    m1 = m(2)*sum(W2.*Vx2); %momentum in Vx component   
    m2 = m(2)*sum(W2.*Vy2); %momentum in Vy component
    E = m(2)*sum(W2.*(Vx2.^2+Vy2.^2)); %total energy
    inside = zeros(Nr^2,1);
    parfor i = 1:Nr^2
        inside(i) = sum(W2.*psi_2d(Vrx2(i) - Vx2,Vry2(i) - Vy2,epsilon2));
    end

    eta = dv2^2*sum(inside.*log(inside)); %entropy
    error_list_s2(nt,5) = rho;
    error_list_s2(nt,6) = m1;
    error_list_s2(nt,7) = m2;
    error_list_s2(nt,8) = E;
    error_list_s2(nt,9) = eta;

end

r = Fixed_Point_list(:,2)>tol;
if max(r) == 1
    disp('Fixed point iteration achieved the maximum number of iterations at time ')
    disp(Fixed_Point_list(r,:))
end

if max(r) == 0
    disp('Fixed point iteration converged at every time step')
end

Total_time = toc;

disp(['Total time = ' num2str(Total_time)])

filename = ['implicit_midpoint_BKW_multi_species_particle_2d_symmetric_n_',num2str(n),'_Example_',num2str(example),'_dv1_',num2str(dv1),'_dv2_',num2str(dv2),'_dt_',num2str(dt),'tmax',num2str(tmax),'.mat'];
save(filename)


