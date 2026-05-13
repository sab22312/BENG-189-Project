%filename: crossbridge.m
%units:
% length: nm = nanometers  (nano=10^-9)
% time:   s  = seconds
% velocity:  nm/s     
% force:  pN = picoNewtons (pico=10^-12; Newton=Kg m/s^2)
clear all
clf
global p1 GAMMA;           %constants in p(x)
global Tstart V;           %constants in v(t)
n0=10000                   %number of crossbridges
alpha=14    % /s  probability per unit time for attachment
beta=126    % /s  probability per unit time for detachment
klokmax=3000               %total number of time steps
dt=0.01/(alpha+beta)       % s  duration of time step
A=5         % nm  displacement of newly attached crossbridge
%constants used in the fuction p(x): 
p1=4        % pN     
GAMMA=0.322 % /nm
%test the function p(x) to make sure p(0)=0
%to within 10*eps, where eps=machine precision:
if(abs(p(0)/p(A))>10*eps)
  error('p(0) is not zero')
end
%constants used in the function v(t):
Tstart=10/(alpha+beta) % s     time at which motion begins
V=500                  % nm/s  shortening velocity for t>Tstart
%all crossbridges are initially detached:
a=zeros(1,n0);
x=zeros(1,n0);
%number of bins for population density function:
nbins=round(sqrt(n0*alpha/(alpha+beta)));
%
for klok=1:klokmax          %loop over time steps
  t=klok*dt;                %current time
  dx= -v(t)*dt;             %change in half-sarcomere length
%Note that the following ``for'' loop is commented out.
%It explains the logic of the program but would execute slowly.
%The code actually used is a vectorized version which appears
%immediately below the ``for'' loop.
% for i=1:n0                %loop over crossbridges
%   if(a(i))                %if crossbridge is attached
%     x(i)=x(i)+dx;         %crossbridge displacement
%     if(rand<beta*dt)      %if crossbridge detaches
%       a(i)=0;             %record its state as detached
%       x(i)=0;             %and its displacement as zero
%     end
%   else                    %crossbridge is detached
%     if(rand<alpha*dt)     %if crossbridge attaches
%       a(i)=1;             %record its state as attached
%       x(i)=A;             %and its displacement as A
%     end
%   end
% end
%The following vector statements are equivalent 
%to the above commented-out ``for'' loop:
  x(find(a))=x(find(a))+dx;   %displace only attached crossbridges
  prob=(beta*dt)*a+(alpha*dt)*(1-a);%probability of changing state
  change=rand(1,n0)<prob;  %decide which crossbridges change state
  a=xor(change,a);         %change the state of those crossbridges
  x(find(a&change))=A;        %for newly attached bridges, set x=A
  x(find(~a))=0;                %for all detached bridges, set x=0
%store results for future plotting:
  tsave(klok)=t;          %current time
  Usave(klok)=sum(a)/n0;  %fraction attached
  Psave(klok)=sum(p(x));  %force
  P=Psave(klok)           %write force on screen
end
figure(1)
subplot(2,1,1),plot(tsave,Usave)
subplot(2,1,2),plot(tsave,Psave)
figure(2)
[u,xc]=popdens(x,a,nbins);  %evaluate population density
plot(xc,u,'*')              %plot population density
