function vv=v(t)
%filename v.m
%prescribed velocity of shortening at time t
%isometric contraction for t<Tstart
%shortening velocity V for t>Tstart
global Tstart V;
if (t<Tstart)
  vv=0;
else
  vv=V;
end