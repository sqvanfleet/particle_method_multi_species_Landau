function [f1,f2] = multi_species_exact_2d(t,vx1,vy1,vx2,vy2,beta,m,n)


C = 1/2; %Integration constant and beta from BKW solution

d = 2; %dimension
K = 1 - C*exp(-2*beta*(d-1)*t);
Q = (1-K)/(2*K);



f1(:,:) = n(1)*(m(1)/(2*pi*K))^(d/2)*exp(-m(1)*(vx1.^2+vy1.^2)/(2*K)).*...
    (1-d*Q+m(1)*Q*(vx1.^2+vy1.^2)/K);

f2(:,:) = n(2)*(m(2)/(2*pi*K))^(d/2)*exp(-m(2)*(vx2.^2+vy2.^2)/(2*K)).*...
    (1-d*Q+m(2)*Q*(vx2.^2+vy2.^2)/K);






