}{clc
clear all
close all

% Coulumb case
gamma = -3;

Vmax1 = 4;
Vmax2 = 4;

%Example 9 parameters
m = [1,1];
n1 = [1,1];
u1 = [.5,.25];
u2 = [-.25,0];
T = [.25,.125];
B1 = [1/32,1/32;1/32,1/32];

% number of particles in one direction
n = 40;
% total number of particles
Np = n^2; 

% initial particle velocites species 1
dv1 = 2*Vmax1/n;
%List of cell centers species 1 x direction
v1x = (-Vmax1+u1(1)+dv1/2):dv1:(Vmax1+u1(1)-dv1/2);
%List of cell centers species 1 y direction
v1y = (-Vmax1+u1(2)+dv1/2):dv1:(Vmax1+u1(2)-dv1/2);
[vx1,vy1] = ndgrid(v1x,v1y); % species 1 initial velocities
Vx1 = vx1(:); Vy1 = vy1(:); %species 1 initial velocities

% initial particle velocities species 2
dv2 = 2*Vmax2/n;
%list of cell centers species 2 x direction
v2x = (-Vmax2+u2(1)+dv2/2):dv2:(Vmax2+u2(1)-dv2/2); 
%list of cell centers species 2 y direction
v2y = (-Vmax2+u2(2)+dv2/2):dv2:(Vmax2+u2(2)-dv2/2); 
[vx2,vy2] = ndgrid(v2x,v2y); % species 2 initial velocities
Vx2 = vx2(:); Vy2 = vy2(:); %species 2 initial velocities

% Initial Conditions
[f1_0,f2_0,example] = Example_9_coulumb_case_2d(vx1,vy1,vx2,vy2,m,n1,u1,u2,T);
% initial weights
t0 = 0; 
w1 = dv1^2*f1_0; w2 = dv2^2*f2_0;
W1 = w1(:); % Weights for species one
W2 = w2(:); % Weights for species two

%Check initial moments
n0 = [sum(W1),sum(W2)];
u0 = [sum(W1.*Vx1),sum(W1.*Vy1);...
        sum(W2.*Vx2),sum(W2.*Vy2)];

rho = [m(1)*n0(1),m(2)*n0(2)];
total_rho = sum(rho);

u0(1,:) = (1/n0(1))*u0(1,:);
u0(2,:) = (1/n0(2))*u0(2,:);

T0 = [sum(W1.*((Vx1-u0(1,1)).^2+(Vy1-u0(1,2)).^2)),...
    sum(W2.*((Vx2-u0(2,1)).^2+(Vy2-u0(2,2)).^2))];

T0 = (m./(2*n0)).*T0;

total_n = sum(n0);

ux_relax = (1/total_rho)*(rho(1)*u0(1,1)+rho(2)*u0(2,1));
uy_relax = (1/total_rho)*(rho(1)*u0(1,2)+rho(2)*u0(2,2));

E0 = [m(1)*sum(W1.*(Vx1.^2+Vy1.^2)),m(2)*sum(W2.*(Vx2.^2+Vy2.^2))];
total_E = sum(E0);

T_relax = (1/(2*total_n))*(total_E-total_rho*(ux_relax^2+uy_relax^2));

% reconstruction mesh species 1
Nr = n;
dvr1 = 2*Vmax1/Nr;
% Cell center x direction
vr1x = (-Vmax1+u1(1)+dvr1/2):dvr1:(Vmax1+u1(1)-dvr1/2); 
% Cell center y direction
vr1y = (-Vmax1+u1(2)+dvr1/2):dvr1:(Vmax1+u1(2)-dvr1/2);
[vrx1,vry1] = ndgrid(vr1x,vr1y);
Vrx1 = vrx1(:); %cell centers x direction species 1
Vry1 = vry1(:); %cell centers y direction species 2

% reconstruction mesh species 2
Nr = n;
dvr2 = 2*Vmax2/Nr;
% Cell center x direction
vr2x = (-Vmax2+u2(1)+dvr2/2):dvr2:(Vmax2+u2(1)-dvr2/2); 
% Cell center y direction
vr2y = (-Vmax2+u2(2)+dvr2/2):dvr2:(Vmax2+u2(2)-dvr2/2);
[vrx2,vry2] = ndgrid(vr2x,vr2y);
Vrx2 = vrx2(:); %cell centers x direction species 1
Vry2 = vry2(:); %cell centers y direction species 2

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

%initial comparison in the x direction
% figure
% plot(vr1x,f1_0(:,Nr/2),'DisplayName','Initial Condition') 
% hold on
% plot(vr1x,f1(:,Nr/2),'DisplayName','Particle Method')
% legend
% title('f1 initial (:,Nr/2)')
% 
% figure
% plot(vr2x,f2_0(:,Nr/2),'DisplayName','Initial Condition')
% hold on 
% plot(vr2x,f2(:,Nr/2),'DisplayName','Particle Method')
% legend
% title('f2 initial (:,Nr/2)')
% 
% % initial comparison in the y direction
% figure
% plot(vr1y,f1_0(Nr/2,:),'DisplayName','Initial Condition') 
% hold on
% plot(vr1y,f1(Nr/2,:),'DisplayName','Particle Method')
% legend
% title('f1 initial (Nr/2,:)')
% 
% figure
% plot(vr2y,f2_0(Nr/2,:),'DisplayName','Initial Condition')
% hold on 
% plot(vr2y,f2(Nr/2,:),'DisplayName','Particle Method')
% legend
% title('f2 initial (Nr/2,:)')



tmax = 40;
dt = 0.01;
%Nt = 5;
Nt = round((tmax-t0)/dt);
error_list_s1 = zeros(Nt,8);
error_list_s2 = zeros(Nt,8);


inside1 = zeros(Nr^2,1);
inside2 = zeros(Nr^2,1);

Ux1 = zeros(Np,1);
Uy1 = zeros(Np,1);
Ux2 = zeros(Np,1);
Uy2 = zeros(Np,1);

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

        x_comp = (Vx1(p)-Vx1);
        y_comp = (Vy1(p)-Vy1);

        for q = 1:Np
            if min(abs(x_comp(q)) < (100*eps))
                x_comp(q) = 0;
            end
        end

        for q = 1:Np
            if min(abs(y_comp(q)) < (100*eps))
                y_comp(q) = 0;
            end
        end


        len = sqrt(x_comp.^2 + y_comp.^2);
        Const = B1(1,1)/m(1);
        A11 = Const*(len+eps).^gamma.*(len.^2 - x_comp.^2);
        A12 = -Const*(len+eps).^gamma.*x_comp.*y_comp;
        A21 = A12;
        A22 = Const*(len+eps).^gamma.*(len.^2 - y_comp.^2);

        termA = term1_x(p)/m(1) - term1_x/m(1);
        termB = term1_y(p)/m(1) - term1_y/m(1);

        Ux1(p) = -sum(W1.*(A11.*termA + A12.*termB));
        Uy1(p) = -sum(W1.*(A21.*termA + A22.*termB));

    end

    for p = 1:Np

        %A matrix i=1 j=2

        x_comp = (Vx1(p)-Vx2);
        y_comp = (Vy1(p)-Vy2);

        for q = 1:Np
            if min(abs(x_comp(q)) < (100*eps))
                x_comp(q) = 0;
            end
        end

        for q = 1:Np
            if min(abs(y_comp(q)) < (100*eps))
                y_comp(q) = 0;
            end
        end


        len = sqrt(x_comp.^2 + y_comp.^2);
        Const = B1(1,2)/m(1);
        A11 = Const*(len+eps).^gamma.*(len.^2 - x_comp.^2);
        A12 = -Const*(len+eps).^gamma.*x_comp.*y_comp;
        A21 = A12;
        A22 = Const*(len+eps).^gamma.*(len.^2 - y_comp.^2);

        termA = term1_x(p)/m(1) - term2_x/m(2);
        termB = term1_y(p)/m(1) - term2_y/m(2);

        Ux1(p) = Ux1(p) - sum(W2.*(A11.*termA + A12.*termB));
        Uy1(p) = Uy1(p) - sum(W2.*(A21.*termA + A22.*termB));

    end

    for p = 1:Np

        %A matrix i=2 j=1
        x_comp = (Vx2(p)-Vx1);
        y_comp = (Vy2(p)-Vy1);

        for q = 1:Np
            if min(abs(x_comp(q)) < (100*eps))
                x_comp(q) = 0;
            end
        end

        for q = 1:Np
            if min(abs(y_comp(q)) < (100*eps))
                y_comp(q) = 0;
            end
        end


        len = sqrt(x_comp.^2 + y_comp.^2);
        Const = B1(2,1)/m(2);
        A11 = Const*(len+eps).^gamma.*(len.^2 - x_comp.^2);
        A12 = -Const*(len+eps).^gamma.*x_comp.*y_comp;
        A21 = A12;
        A22 = Const*(len+eps).^gamma.*(len.^2 - y_comp.^2);

        termA = term2_x(p)/m(2) - term1_x/m(1);
        termB = term2_y(p)/m(2) - term1_y/m(1);

        Ux2(p) = -sum(W1.*(A11.*termA + A12.*termB));
        Uy2(p) = -sum(W1.*(A21.*termA + A22.*termB));

    end


    for p = 1:Np

        %A matrix i = 2 j = 2
        
        x_comp = (Vx2(p)-Vx2);
        y_comp = (Vy2(p)-Vy2);

        for q = 1:Np
            if min(abs(x_comp(q)) < (100*eps))
                x_comp(q) = 0;
            end
        end

        for q = 1:Np
            if min(abs(y_comp(q)) < (100*eps))
                y_comp(q) = 0;
            end
        end


        len = sqrt(x_comp.^2 + y_comp.^2);
        Const = B1(2,2)/m(2);
        A11 = Const*(len+eps).^gamma.*(len.^2 - x_comp.^2);
        A12 = -Const*(len+eps).^gamma.*x_comp.*y_comp;
        A21 = A12;
        A22 = Const*(len+eps).^gamma.*(len.^2 - y_comp.^2);


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

%     figure(5)
%     plot(vr1x,f1(:,Nr/2),'DisplayName',['f1(:,Nr/2) at time t = ',num2str(time)])
%     title('f1(:,Nr/2')
%     legend
%     drawnow
%     
%     figure(6)
%     plot(vr2x,f2(:,Nr/2),'DisplayName',['f2(:,Nr/2) at time t = ',num2str(time)])
%     title('f2(:,Nr/2)')
%     legend
%     drawnow
% 
%     figure(7)
%     plot(vr1y,f1(Nr/2,:),'DisplayName',['f1(Nr/2,:) at time t = ',num2str(time)])
%     title('f1(Nr/2,:)')
%     legend
%     drawnow
% 
%     figure(8)
%     plot(vr2y,f2(Nr/2,:),'DisplayName',['f2(Nr/2,:) at time t = ',num2str(time)])
%     title('f2(Nr/2,:)')
%     legend
%     drawnow
   

    disp("current time: ");
    disp(time);

    %save the moments and errors at each time step
    error_list_s1(nt,1) = time;
    error_list_s2(nt,1) = time;



    % moments for species 1
    
    nn = sum(W1);%density
    rho = m(1)*nn; %mass
    m1 = (1/nn)*sum(W1.*Vx1); %velocity in Vx component   
    m2 = (1/nn)*sum(W1.*Vy1); %velocity in Vy component
    E = m(1)*sum(W1.*(Vx1.^2+Vy1.^2)); %energy
    Temperature = (m(1)/(2*nn))*sum(W1.*((Vx1-m1).^2 + (Vy1-m2).^2)); %Temperature
    inside = zeros(Nr^2,1);
    for i = 1:Nr^2
        inside(i) = sum(W1.*psi_2d(Vrx1(i) - Vx1,Vry1(i) - Vy1,epsilon1));
    end
    eta = dv1^2*sum(inside.*log(inside)); %entropy

    error_list_s1(nt,2) = nn;
    error_list_s1(nt,3) = rho;
    error_list_s1(nt,4) = m1;
    error_list_s1(nt,5) = m2;
    error_list_s1(nt,6) = E;
    error_list_s1(nt,7) = Temperature;
    error_list_s1(nt,8) = eta;



    % moments for species 2

    nn = sum(W2);
    rho = m(2)*nn; %mass
    m1 = (1/nn)*sum(W2.*Vx2); %momentum in Vx component   
    m2 = (1/nn)*sum(W2.*Vy2); %momentum in Vy component
    E = m(2)*sum(W2.*(Vx2.^2+Vy2.^2)); %energy
    Temperature = (m(2)/(2*nn))*sum(W2.*((Vx2-m1).^2 + (Vy2-m2).^2));
    inside = zeros(Nr^2,1);
    for i = 1:Nr^2
        inside(i) = sum(W2.*psi_2d(Vrx2(i) - Vx2,Vry2(i) - Vy2,epsilon2));
    end
    eta = dv2^2*sum(inside.*log(inside)); %entropy

    
    
    error_list_s2(nt,2) = nn;
    error_list_s2(nt,3) = rho;
    error_list_s2(nt,4) = m1;
    error_list_s2(nt,5) = m2;
    error_list_s2(nt,6) = E;
    error_list_s2(nt,7) = Temperature;
    error_list_s2(nt,8) = eta;


end


% Compute the total moments
% total_moments = zeros(Nt,7);
% total_moments(:,1) = error_list_s1(:,2) + error_list_s2(:,2); %total n 
% total_moments(:,2) = error_list_s1(:,3) + error_list_s2(:,3); %total rho 
% total_moments(:,3) = (1./total_moments(:,2)).*...
%     (error_list_s1(:,3).*error_list_s1(:,4) + error_list_s2(:,3).*error_list_s2(:,4)); %total u x-direction
% total_moments(:,4) = (1./total_moments(:,2)).*...
%     (error_list_s1(:,3).*error_list_s1(:,5) + error_list_s2(:,3).*error_list_s2(:,5)); %total u y-direction
% total_moments(:,5) = error_list_s1(:,6) + error_list_s2(:,6); %total energy
% %total_moments(:,6) = (1./(2*total_moments(:,1))).*...
% %    (m(1)*sum(W1.*((Vx1-total_moments(:,3)).^2+(Vy1-total_moments(:,4)).^2))+...
% %     m(2)*sum(W2.*((Vx2-total_moments(:,3)).^2+(Vy2-total_moments(:,4)).^2))); %Temperature
% total_moments(:,7) = error_list_s1(:,8) + error_list_s2(:,8); %Total Entropy 


filename = ['multi_species_particle_2d_Coulomb_n_',num2str(n),'_Example_',num2str(example),'_gamma_',num2str(gamma),'_dv1_',num2str(dv1),'_dv2_',num2str(dv2),'_dt_',num2str(dt),'tmax',num2str(tmax),'.mat'];
save(filename)




