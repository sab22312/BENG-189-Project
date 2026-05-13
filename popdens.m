function [u,xc]=popdens(x,a,nbins)
%filename popdens.m
%find population density function of attached crossbridges
%input variables:
%  x=array of crossbridge displacements
%  a=array of crossbridge states (0=detached,1=attached)
%  nbins=number of bins 
%output variables:
%  xc=array of (equally spaced) bin centers
%  u=array of crossbridge population densities
n0=length(x);                  %total number of crossbridges
Nabridges=length(x(find(a)));  %number of attached crossbridges
%nabridges=number of attached crossbridges in each bin; 
%xc=centers of the bins:
[nabridges,xc]=hist(x(find(a)),nbins);
%normalize nabridges to find crossbridge population density
xstep=xc(2)-xc(1);
u=(nabridges/Nabridges)/xstep;
