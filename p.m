function pp=p(x)
%filename p.m
%crossbridge force as a function of crossbridge displacement
%x is vector of crossbridge displacements
%pp=p(x) is a vector of corresponding forces
global p1 GAMMA;
pp=p1*(exp(GAMMA*x)-1);
