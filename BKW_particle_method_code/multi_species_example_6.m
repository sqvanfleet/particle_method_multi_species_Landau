function [Vmax1,Vmax2,B,beta,m,n,example] = multi_species_example_6()

% mass ratio 2.  With beta = 1/16 we have a similar relaxation rate as the
% single species examples.  The difference between example 5 and example 6
% is that we will choose a different Vmax1.  In otherwords we will use
% different epsilon valuse for each species.  Here we use a smaller Vmax1
% for the heavier species.

example = 6;
m = [2,1];
Vmax1 = 3;
Vmax2 = 4;
B = [1/8,1/16;1/16,1/32];
n = [1,1];
beta = 1/16;