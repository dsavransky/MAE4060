function dv = circular_phasing_animation(gam,jtarg)
%circular_phasing_animation demonstrates rendezvous on a circular orbit
%
%   circular_phasing_animation animates rendezvous on a circular orbit and
%   returns the total delta v.
%
%   circular_phasing_animation(gam) animates rendezvous between an
%   intercepter and target separated by gam degrees (positive for target
%   leading). Defaults to 30 degrees.
%
%   circular_phasing_animation(...,jtarg) performs the rendezvous in jtarg
%   orbits of the target spacecraft. Must be >= 1.  Defaults to 1.

% Copyright (c) 2016 Dmitry Savransky (ds264@cornell.edu)

if ~exist('gam','var') || isempty(gam)
    gam = 30;
end

if ~exist('jtarg','var') || isempty(jtarg)
    jtarg = 1;
end

assert(jtarg >= 1, 'jtarg must be >= 1.')

gam = gam*pi/180;

mu = 1; %gravitational parameter
a = 3;  %semi-major axis
n = sqrt(mu/a^3); %target orbit

if gam > pi, gam = gam-2*pi; end

th = linspace(0,2*pi,1000);
dt = (th(2)-th(1))/n; %time step on target orbit
pint = 3/length(th); %want one orbit to take 3 seconds to animate

Tphase = (2*pi*jtarg - gam)/n;
aphase = (mu*(Tphase/(2*pi))^2)^(1/3);
nphase = 2*pi/Tphase;
ephase = abs(a/aphase - 1);
dv = 2*sqrt(2*mu/a - mu/aphase) - sqrt(mu/a);

%if phasing orbit is larger, you start at periapsis, if it's smaller you
%start at apoapsis.
if aphase < a
    Mphase = mod(nphase*(0:dt:Tphase)+pi,2*pi);
    forbind = length(th)/2;
    dvsign = 1;
else
    Mphase = mod(nphase*(0:dt:Tphase),2*pi);
    forbind = length(th);
    dvsign = -1;
end

nuphase = invKepler(Mphase,ephase);
rphase = (aphase*(1 - ephase^2))./(1 + ephase*cos(nuphase));

figure(1)
clf
plot(a*cos(th),a*sin(th),'b--')
hold on
targ = plot(a*cos(th(1)+gam),a*sin(th(1)+gam),'b.','MarkerSize',60);
intr = plot(a*cos(th(1)),a*sin(th(1)),'r.','MarkerSize',40);
introrbit = plot(rphase.*cos(nuphase),rphase.*sin(nuphase),'r--');
axis equal off
set(gca,'XLimMode','manual','YLimMode','manual')
set(introrbit,'Visible','off')

pause(1)

%one orbit of following
for j = 2:forbind
    set(intr,'XData',a*cos(th(j)),'YData',a*sin(th(j)))
    set(targ,'XData',a*cos(th(j)+gam),'YData',a*sin(th(j)+gam))
    pause(pint)
end

set(introrbit,'Visible','on')
dv1 = quiver(a*cos(th(j)),a*sin(th(j)),0,dvsign,'LineWidth',2);
pause(0.5)
delete(dv1)

st = j;
for k = 1:length(nuphase)
    ind = mod(k+st,length(th))+1;
    set(intr,'XData',rphase(k).*cos(nuphase(k)),'YData',rphase(k).*sin(nuphase(k)))
    set(targ,'XData',a*cos(th(ind)+gam),'YData',a*sin(th(ind)+gam))
    pause(pint)
end

dv2 = quiver(a*cos(th(j)),a*sin(th(j)),0,-dvsign,'LineWidth',2);
pause(0.5)
delete(dv2)

set(introrbit,'Visible','off')
st = ind;
for k = 1:length(th)/2
    ind = mod(k+st,length(th))+1;
    set(intr,'XData',a*cos(th(ind)+gam),'YData',a*sin(th(ind)+gam))
    set(targ,'XData',a*cos(th(ind)+gam),'YData',a*sin(th(ind)+gam))
    pause(pint)
end