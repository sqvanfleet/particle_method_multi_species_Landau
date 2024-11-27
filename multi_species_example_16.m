function [Vmax1,Vmax2,B,beta,m,n,example] = multi_species_example_16()

% mass ratio 10.  With beta = 1/16 we have a similar relaxation rate as the
% single species examples.

example = 16;
m = [20,1];
Vmax1 = 4; 
Vmax2 = 4;
B = [1/2,49/40;49/40,1/800];
n = [1,1];
beta = 1/16;