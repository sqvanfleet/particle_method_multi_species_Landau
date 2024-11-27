function [f1,f2,example] = Example_9_coulumb_case_2d(vx1,vy1,vx2,vy2,m,n1,u1,u2,T)

example = 9;


f1 = n1(1)*(m(1)/(2*pi*T(1)))*...
    exp(-m(1)*((vx1-u1(1)).^2+(vy1-u1(2)).^2)/(2*T(1)));


f2 = n1(2)*(m(2)/(2*pi*T(2)))*...
    exp(-m(2)*((vx2-u2(1)).^2+(vy2-u2(2)).^2)/(2*T(2)));