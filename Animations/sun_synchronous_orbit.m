function sun_synchronous_orbit(varargin)
% Animation for exploring a sun synchronous orbit
% 
%      sun_synchronous_orbit, by itself, animates a sun synchronous orbit
%      over 120 solar days
%
%      solar_vs_sidereal_animation('rotation','on'|'off')
%      toggles nodal regression on (default) or off.
%
% Notes:
%      The perspective of the animation is a viewer situated above the
%      Earth and moving with the Earth on its orbit, but not rotating with
%      the Earth about its axis.  The Earth is oriented with the pole axis
%      along the viewer's vertical axis, so that the frame is effectively
%      equatorial, with the ecliptic oriented at approximately 23.4 degrees
%      from the equator.
%
% Examples:
%      sun_synchronous_orbit()
%      sun_synchronous_orbit('rotation','off')

% Copyright (c) 2020 Dmitry Savransky (ds264@cornell.edu)

%input parsing
p = inputParser;
addParameter(p,'rotation', 'on', @(x) strcmpi(x,'on') || strcmpi(x,'off'));
addParameter(p,'figure', 1,  @(x) isnumeric(x) && floor(x) == x && x>0)
parse(p,varargin{:});

R_E = 6378.1366*1000; %Earth radius (m)
mAU = 149597870700; %AU in m

a = 1.00000261;  %semi-major axis (AU)
a = a*mAU/R_E; %Earth radii

e = 0.01671123;  %eccentricity
I = -0.00001531 * pi/180.; %inclination (rad)
L = 100.46457166 * pi/180; %mean longitude (rad)
varpi = 102.93768193 * pi/180; %long. of periapsis (rad)
Omega = 0; %long. of asending node (rad)
omega = varpi - Omega; %argument of periapsis (rad)

musun = 1.32712440018e20; %m^3/s^2
sun_EM = 328900.56; %mass ratio: sun/(Earth+Moon)
mu = (1 + 1/sun_EM)*musun/(R_E^3); %Earth radii^3/s^2;
n = sqrt(mu/a^3); %s^(-1)
R_S = 6.957e8; %Sun radius (m)

tropyear =  365.242190402*86400; %tropical year in seconds

t = linspace(0,tropyear,500000);

M = n*t;
[E,nu] = invKepler(M,e);
E = E.';

A = [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega));...
    a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega));...
    a.*sin(I).*sin(omega)];

B = [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + ...
    sin(Omega).*cos(I).*cos(omega));...
    a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + ...
    cos(Omega).*cos(I).*cos(omega));...
    a.*sqrt(1-e.^2).*sin(I).*cos(omega)];

r = A*(cos(E) - e) + B*(sin(E)); %orbit path (r_E/Sun) in ecliptic coords

vareps = 84381.412/3600*pi/180; %as->rad
rotMat = [1 0 0; 0 cos(vareps) -sin(vareps); 0 sin(vareps) cos(vareps)]; %eclip->equat
rsun = -rotMat*r; %r_Sun/E in equatorial coords
%% sat orbit

w_e = 7.2921150e-5; %Earth rotation (rad/s)
mue = 3.986004418e5/(R_E/1000)^3; %R_E^3/s^2

a = 1 + 1000*1000/R_E; %1000 km orbit (R_E)
e = 0.1;  %eccentricity
J2 = 1082.6e-6;
n = sqrt(mue/a^3); %s^(-1)
Omega0 = -30*pi/180;
omega = 0;

Omegadot = 2*pi/tropyear;
I = acos( -2/3*Omegadot*a^2*(1 - e^2)/J2/n );

if strcmp(p.Results.rotation,'on')
    Omega = Omega0 + Omegadot*t;
else
    Omega = Omega0;
end

M = mod(n*t,2*pi);
E = invKepler(M,e);
E = E.';

A = [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega));...
    a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega));...
    ones(size(Omega))*a.*sin(I).*sin(omega)];

B = [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + ...
    sin(Omega).*cos(I).*cos(omega));...
    a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + ...
    cos(Omega).*cos(I).*cos(omega));...
    ones(size(Omega))*a.*sqrt(1-e.^2).*sin(I).*cos(omega)];

r = A.*repmat((cos(E) - e),3,1) + B.*repmat(sin(E),3,1); %orbit

%%
earth = imread('landOcean.jpg'); %read Earth map (should be in MATLAB demos dir)

f1 = figure(1);
clf(f1)

[x,y,z] = sphere(100);
props.AmbientStrength = 0;  %no ambient light
props.DiffuseStrength = 1;  %full diffuse illumination
props.SpecularColorReflectance = 0.5;
props.SpecularExponent = 1;
props.SpecularStrength = 0.45;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.Cdata = earth;
Earth = surface(x,y,flip(z),props);
hold on
g = hgtransform;
%hold off
Earth.Parent = g;
view(50,10)
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
axis equal off

% Add sun
sun = light('Position',rsun(:,1),'Style','local'); %light is effective point source at sun location

orbpath = plot3(r(1,1),r(2,1),r(3,1),'b');
orbloc = plot3(r(1,1),r(2,1),r(3,1),'b.','MarkerSize',20);

%%
for j = 2:500:length(t)
    set(g,'Matrix',makehgtform('zrotate',w_e*t(j)))
    pause(1/30)
    set(sun,'Position',rsun(:,j))
    sind = max([1,j-1000]);
    set(orbpath,'XData',r(1,sind:j),'YData',r(2,sind:j),'ZData',r(3,sind:j))
    set(orbloc,'XData',r(1,j),'YData',r(2,j),'ZData',r(3,j))
end
