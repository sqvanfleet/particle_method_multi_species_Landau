function [Vmax1,Vmax2,B,beta,m,n,example] = multi_species_example_11()

% mass ratio 100.  With beta = 1/16 we have a similar relaxation rate as the
% single species examples.

example = 11;
m = [100,1];
Vmax1 = .4; 
Vmax2 = 4;
B = [1/2,1249/200;1249/200,1/20000];
n = [1,1];
beta = 1/16;