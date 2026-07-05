function Y = CA_Weight(X, C)


k = 0.1;

A = 2.5;


Y =  A./(1+exp(-k*(X-C)));





end