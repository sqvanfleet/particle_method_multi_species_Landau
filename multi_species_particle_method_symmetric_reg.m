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
% [Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_6();

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
[Vmax1,Vmax2,B1,beta,m,n1,example] = multi_species_example_17();

n = 20;
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
for i = 1:Nr
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
dt = 0.001;
Nt = round((tmax-t0)/dt);
error_list_s1 = zeros(Nt,9);
error_list_s2 = zeros(Nt,9);


for nt = 1:Nt
    
    inside1 = zeros(Nr^2,1);
    inside2 = zeros(Nr^2,1);

    for p = 1:Np
        inside1(p) = sum(W1.*psi_2d(Vrx1(p)-Vx1,Vry1(p)-Vy1,epsilon1));
        inside2(p) = sum(W2.*psi_2d(Vrx2(p)-Vx2,Vry2(p)-Vy2,epsilon2));
    end
    inside1 = log(inside1);
    inside2 = log(inside2);

    term1_x = zeros(Np,1);
    term1_y = zeros(Np,1);
    term2_x = zeros(Np,1);
    term2_y = zeros(Np,1);
    for p = 1:Np
        [A,B] = gpsi_2d(Vx1(p)-Vrx1,Vy1(p)-Vry1,epsilon1);
        term1_x(p) = dv1^2*sum(A.*inside1);
        term1_y(p) = dv1^2*sum(B.*inside1);


        [A,B] = gpsi_2d(Vx2(p)-Vrx2,Vy2(p)-Vry2,epsilon2);
        term2_x(p) = dv2^2*sum(A.*inside2);
        term2_y(p) = dv2^2*sum(B.*inside2);
    end

    Ux1 = zeros(Np,1);
    Uy1 = zeros(Np,1);
    Ux2 = zeros(Np,1);
    Uy2 = zeros(Np,1);
    for p = 1:Np
        %A matrix i = 1 and j = 1
        len = sqrt((Vx1(p)-Vx1).^2 +(Vy1(p)-Vy1).^2);
        Const = B1(1,1)/m(1);
        A11 = Const*(len.^2 - (Vx1(p)-Vx1).^2);
        A12 = -Const*(Vx1(p)-Vx1).*(Vy1(p)-Vy1);
        A21 = A12;
        A22 = Const*(len.^2 - (Vy1(p) - Vy1).^2);

        termA = term1_x(p)/m(1) - term1_x/m(1);
        termB = term1_y(p)/m(1) - term1_y/m(1);

        Ux1(p) = -sum(W1.*(A11.*termA + A12.*termB));
        Uy1(p) = -sum(W1.*(A21.*termA + A22.*termB));

        %A matrix i=1 j=2
        len = sqrt((Vx1(p)-Vx2).^2 + (Vy1(p) - Vy2).^2);
        Const = B1(1,2)/m(1);
        A11 = Const*(len.^2 - (Vx1(p) - Vx2).^2);
        A12 = -Const*(Vx1(p) - Vx2).*(Vy1(p) - Vy2);
        A21 = A12;
        A22 = Const*(len.^2 - (Vy1(p) - Vy2).^2);

        termA = term1_x(p)/m(1) - term2_x/m(2);
        termB = term1_y(p)/m(1) - term2_y/m(2);

        Ux1(p) = Ux1(p) - sum(W2.*(A11.*termA + A12.*termB));
        Uy1(p) = Uy1(p) - sum(W2.*(A21.*termA + A22.*termB));

        %A matrix i=2 j=1
        len = sqrt((Vx2(p) - Vx1).^2 + (Vy2(p) - Vy1).^2);
        Const = B1(2,1)/m(2);
        A11 = Const*(len.^2 - (Vx2(p) - Vx1).^2);
        A12 = -Const*(Vx2(p) - Vx1).*(Vy2(p) - Vy1);
        A21 = A12;
        A22 = Const*(len.^2 - (Vy2(p) - Vy1).^2);

        termA = term2_x(p)/m(2) - term1_x/m(1);
        termB = term2_y(p)/m(2) - term1_y/m(1);

        Ux2(p) = -sum(W1.*(A11.*termA + A12.*termB));
        Uy2(p) = -sum(W1.*(A21.*termA + A22.*termB));

        %A matrix i = 2 j = 2
        len = sqrt((Vx2(p) - Vx2).^2 + (Vy2(p) - Vy2).^2);
        Const = B1(2,2)/m(2);
        A11 = Const*(len.^2 - (Vx2(p) - Vx2).^2);
        A12 = -Const*(Vx2(p) - Vx2).*(Vy2(p) - Vy2);
        A21 = A12;
        A22 = Const*(len.^2 - (Vy2(p) - Vy2).^2);

        termA = term2_x(p)/m(2) - term2_x/m(2);
        termB = term2_y(p)/m(2) - term2_y/m(2);

        Ux2(p) = Ux2(p) - sum(W2.*(A11.*termA + A12.*termB));
        Uy2(p) = Uy2(p) - sum(W2.*(A21.*termA + A22.*termB));

    end

    Vx1 = Vx1 + dt*Ux1;
    Vy1 = Vy1 + dt*Uy1;

    Vx2 = Vx2 + dt*Ux2;
    Vy2 = Vy2 + dt*Uy2;
    
    %reconstruction
    for i = 1:Nr
        for j = 1:Nr
            f1(i,j) = sum(W1.*psi_2d(vrx1(i,j)-Vx1,vry1(i,j)-Vy1,epsilon1));
            f2(i,j) = sum(W2.*psi_2d(vrx2(i,j)-Vx2,vry2(i,j)-Vy2,epsilon2));
        end
    end


    % plot
    time = t0+dt*nt;
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
    for i = 1:Nr^2
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
    for i = 1:Nr^2
        inside(i) = sum(W2.*psi_2d(Vrx2(i) - Vx2,Vry2(i) - Vy2,epsilon2));
    end

    eta = dv2^2*sum(inside.*log(inside)); %entropy
    error_list_s2(nt,5) = rho;
    error_list_s2(nt,6) = m1;
    error_list_s2(nt,7) = m2;
    error_list_s2(nt,8) = E;
    error_list_s2(nt,9) = eta;

end

% figure
% plot(error_list_s1(:,1),error_list_s1(:,4),'DisplayName',['n = ',num2str(n)]);
% xlabel('Time')
% ylabel('L2 Error')
% title('L2 error plot Species 1')
% legend
% hold off
% 
% %%%% Plot of the L2 error versus time species 2
% 
% figure
% plot(error_list_s2(:,1),error_list_s2(:,4),'DisplayName',['n = ',num2str(n)]);
% xlabel('Time')
% ylabel('L2 Error')
% title('L2 error plot Species 1')
% legend
% hold off

% Compute the total moments
% total_moments = zeros(Nt,5);
% total_moments(:,1) = error_list_s1(:,5) + error_list_s2(:,5); %total mass
% total_moments(:,2) = error_list_s1(:,6) + error_list_s2(:,6); %total momentum Vx component
% total_moments(:,3) = error_list_s1(:,7) + error_list_s2(:,7); %total momentum Vy component
% total_moments(:,4) = error_list_s1(:,8) + error_list_s2(:,8); %total energy
% total_moments(:,5) = error_list_s1(:,9) + error_list_s2(:,9); %total entropy

filename = ['multi_species_particle_2d_symmetric_n_',num2str(n),'_Example_',num2str(example),'_dv1_',num2str(dv1),'_dv2_',num2str(dv2),'_dt_',num2str(dt),'tmax',num2str(tmax),'.mat'];
save(filename)


