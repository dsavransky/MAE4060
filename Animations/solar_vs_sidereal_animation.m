function solar_vs_sidereal_animation(varargin)
% Animation for exploring solar vs sidereal time
% 
%      solar_vs_sidereal_animation, by itself, animates the rotation of the
%      Earth over 1 sidereal day.
%
%      solar_vs_sidereal_animation(n) animates the rotation of the Earth
%      over n days.
%
%      solar_vs_sidereal_animation(...,'daytype','solar'|'sidereal')
%      toggles between solar and sidereal days (default sidereal).
%
%      solar_vs_sidereal_animation(...,'animlength',n) sets the duration of
%      the full animation to n seconds (default 5).
%
%      solar_vs_sidereal_animation(...,'figure',f) creates the animation in
%      figure f (default 1).
%
%      solar_vs_sidereal_animation(...,'oneyear',true) animates for exactly
%      one year (defined by the mean (J2000) orbital characteristics of
%      Earth's orbit) 
%
% Notes:
%      The perspective of the animation is a viewer situated above the
%      Earth and moving with the Earth on its orbit, but not rotating with
%      the Earth about its axis.  The Earth is oriented with the pole axis
%      along the viewer's vertical axis, so that the frame is effectively
%      equatorial, with the ecliptic oriented at approximately 23.4 degrees
%      from the equator. The red axes represent the original orientation of
%      the Earth, and pass through the equator, prime meridian, and north
%      pole. The blue axes rotate with the Earth throughout the animation.
%      A rotation of one sidereal day will leave the two sets of axes
%      aligned, whereas one solar day will result in approximately a one
%      degree misalignment between the axes.
%
% Examples:
%      solar_vs_sidereal_animation(1) % 1 sidereal day
%      solar_vs_sidereal_animation(1,'daytype','solar') % 1 solar day
%      solar_vs_sidereal_animation(30,'daytype','solar') % 1 solar month
%      solar_vs_sidereal_animation(365) % 1 year (approx)
%      solar_vs_sidereal_animation('oneyear',true,'figure',2) % 1 year (exact)

% Copyright (c) 2018 Dmitry Savransky (ds264@cornell.edu)

%input parsing
p = inputParser;
addOptional(p,'ndays', 1, @(x) isnumeric(x) && x>=0);
addParameter(p,'daytype', 'sidereal', @(x) strcmpi(x,'sidereal') || strcmpi(x,'solar'));
addParameter(p,'animlength', 5, @(x) isnumeric(x) && x>=0);
addParameter(p,'oneyear', false, @(x) islogical(x));
addParameter(p,'figure', 1,  @(x) isnumeric(x) && floor(x) == x && x>0)
parse(p,varargin{:});
nframes = 30*p.Results.animlength; %30 fps 

%constants
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

w_e = 7.2921150e-5; %Earth rotation (rad/s) 

musun = 1.32712440018e20; %m^3/s^2
sun_EM = 328900.56; %mass ratio: sun/(Earth+Moon)
mu = (1 + 1/sun_EM)*musun/(R_E^3); %Earth rad^3/s^2;
n = sqrt(mu/a^3); %s^(-1)
T_p = 2*pi/n; %orbital period in s

if p.Results.oneyear
    ndays = T_p/86400;
    daytype = 'solar';
else
    ndays = p.Results.ndays;
    daytype = p.Results.daytype;
end

%define time array
if strcmpi(daytype, 'sidereal')
    t = linspace(0,ndays*2*pi/w_e,nframes); %sidereal days
else
    t = linspace(0,ndays*86400,nframes); %solar day
end

M = n*t;
E = invKepler(M,e);
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

earth = imread('landOcean.jpg'); %read Earth map (should be in MATLAB demos dir)

f1 = figure(p.Results.figure);
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
oax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.5,0,0],[0,1.5,0],[0,0,1.5],0,'Linewidth',2,'color','r');
eax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.5,0,0],[0,1.5,0],[0,0,1.5],0,'Linewidth',2,'Parent',g,'color','b');
hold off
Earth.Parent = g;
view(50,10)
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
axis equal off

% Add sun
sun = light('Position',rsun(:,1),'Style','local'); %light is effective point source at sun location

if ndays > 0
    for j = 2:length(t)
        set(g,'Matrix',makehgtform('zrotate',w_e*t(j)))
        pause(1/30)
        set(sun,'Position',rsun(:,j))
        axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
    end
end