function [Ux1,Uy1,Ux2,Uy2] = right_hand_side_implicit_mid_multi_species_parallel(W1,...
    W2,Vx1,Vx2,Vy1,Vy2,Vx1_new,Vx2_new,Vy1_new,Vy2_new,Vrx1,Vrx2,...
    Vry1,Vry2,dv1,dv2,B1,m,Np,epsilon1,epsilon2)


%midpoint calculation
Vmid_x1 = 0.5*(Vx1_new + Vx1); Vmid_x2 = 0.5*(Vx2_new + Vx2);
Vmid_y1 = 0.5*(Vy1_new + Vy1); Vmid_y2 = 0.5*(Vy2_new + Vy2);

inside1 = zeros(Np,1);
inside2 = zeros(Np,1);

parfor p = 1:Np
    inside1(p) = sum(W1.*psi_2d(Vrx1(p)-Vmid_x1,Vry1(p)-Vmid_y1,epsilon1));
    inside2(p) = sum(W2.*psi_2d(Vrx2(p)-Vmid_x2,Vry2(p)-Vmid_y2,epsilon2));
end
inside1 = log(inside1);
inside2 = log(inside2);

term1_x = zeros(Np,1);
term1_y = zeros(Np,1);
term2_x = zeros(Np,1);
term2_y = zeros(Np,1);
parfor p = 1:Np
    [A,B] = gpsi_2d(Vmid_x1(p)-Vrx1,Vmid_y1(p)-Vry1,epsilon1);
    term1_x(p) = dv1^2*sum(A.*inside1);
    term1_y(p) = dv1^2*sum(B.*inside1);


    [A,B] = gpsi_2d(Vmid_x2(p)-Vrx2,Vmid_y2(p)-Vry2,epsilon2);
    term2_x(p) = dv2^2*sum(A.*inside2);
    term2_y(p) = dv2^2*sum(B.*inside2);
end

Ux1 = zeros(Np,1);
Uy1 = zeros(Np,1);
Ux2 = zeros(Np,1);
Uy2 = zeros(Np,1);
parfor p = 1:Np
    %A matrix i = 1 and j = 1
    len = sqrt((Vmid_x1(p)-Vmid_x1).^2 +(Vmid_y1(p)-Vmid_y1).^2);
    Const = B1(1,1)/m(1);
    A11 = Const*(len.^2 - (Vmid_x1(p)-Vmid_x1).^2);
    A12 = -Const*(Vmid_x1(p)-Vmid_x1).*(Vmid_y1(p)-Vmid_y1);
    A21 = A12;
    A22 = Const*(len.^2 - (Vmid_y1(p) - Vmid_y1).^2);

    termA = term1_x(p)/m(1) - term1_x/m(1);
    termB = term1_y(p)/m(1) - term1_y/m(1);

    Ux1(p) = -sum(W1.*(A11.*termA + A12.*termB));
    Uy1(p) = -sum(W1.*(A21.*termA + A22.*termB));

    %A matrix i=1 j=2
    len = sqrt((Vmid_x1(p)-Vmid_x2).^2 + (Vmid_y1(p) - Vmid_y2).^2);
    Const = B1(1,2)/m(1);
    A11 = Const*(len.^2 - (Vmid_x1(p) - Vmid_x2).^2);
    A12 = -Const*(Vmid_x1(p) - Vmid_x2).*(Vmid_y1(p) - Vmid_y2);
    A21 = A12;
    A22 = Const*(len.^2 - (Vmid_y1(p) - Vmid_y2).^2);

    termA = term1_x(p)/m(1) - term2_x/m(2);
    termB = term1_y(p)/m(1) - term2_y/m(2);

    Ux1(p) = Ux1(p) - sum(W2.*(A11.*termA + A12.*termB));
    Uy1(p) = Uy1(p) - sum(W2.*(A21.*termA + A22.*termB));

    %A matrix i=2 j=1
    len = sqrt((Vmid_x2(p) - Vmid_x1).^2 + (Vmid_y2(p) - Vmid_y1).^2);
    Const = B1(2,1)/m(2);
    A11 = Const*(len.^2 - (Vmid_x2(p) - Vmid_x1).^2);
    A12 = -Const*(Vmid_x2(p) - Vmid_x1).*(Vmid_y2(p) - Vmid_y1);
    A21 = A12;
    A22 = Const*(len.^2 - (Vmid_y2(p) - Vmid_y1).^2);

    termA = term2_x(p)/m(2) - term1_x/m(1);
    termB = term2_y(p)/m(2) - term1_y/m(1);

    Ux2(p) = -sum(W1.*(A11.*termA + A12.*termB));
    Uy2(p) = -sum(W1.*(A21.*termA + A22.*termB));

    %A matrix i = 2 j = 2
    len = sqrt((Vmid_x2(p) - Vmid_x2).^2 + (Vmid_y2(p) - Vmid_y2).^2);
    Const = B1(2,2)/m(2);
    A11 = Const*(len.^2 - (Vmid_x2(p) - Vmid_x2).^2);
    A12 = -Const*(Vmid_x2(p) - Vmid_x2).*(Vmid_y2(p) - Vmid_y2);
    A21 = A12;
    A22 = Const*(len.^2 - (Vmid_y2(p) - Vmid_y2).^2);

    termA = term2_x(p)/m(2) - term2_x/m(2);
    termB = term2_y(p)/m(2) - term2_y/m(2);

    Ux2(p) = Ux2(p) - sum(W2.*(A11.*termA + A12.*termB));
    Uy2(p) = Uy2(p) - sum(W2.*(A21.*termA + A22.*termB));

end

