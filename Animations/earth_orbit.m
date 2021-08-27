function earth_orbit(h,e,O,I,w,az,el)
%Generate plot of an Earth orbit
%
%INPUT
%   h (float): Orbit altitude in km.  For eccentric orbits, this will be
%              the altitude of the semi-major axis. 
%   e (float): Eccentricity
%   O (float): Longitude of the ascending node in degrees
%   I (float): Inclination in degrees
%   w (float): Argument of periapsis in degrees
%   az(float): View azimuth in degrees
%   el(float): View elevation in degrees
%
%OUTPUT
%  None
%

% Copyright (c) 2021 Dmitry Savransky (ds264@cornell.edu)


% h = 600;
% e = 0.1;
% O = 30*pi/180;
% I = 40*pi/180;
% w = 60*pi/180;

O = O*pi/180;
I = I*pi/180;
w = w*pi/180;

R_E = 6378.1366; %Earth radius (km)
mue = 3.986004418e5; %km^3/s^2  gravitational parameter of the Earth
rotMats = DCMs();
earth = imread('landOcean.jpg'); %read Earth map (should be in MATLAB demos dir)

a = h+R_E;

M = linspace(0,2*pi,101);
[~,nu] = invKepler(M,e);

r = a*(1-e^2)./(1+ e*cos(nu));

rvec = [r.*cos(nu), r.*sin(nu), zeros(length(nu),1)].';

iCp = (rotMats{3}(w)*rotMats{1}(I)*rotMats{3}(O)).';
rvec = iCp*rvec;

zp = mean(rvec,2); %offset such that orbit is centered

%figure(1)
figure('units','normalized','outerposition',[0 0 1 1])
clf

[x,y,z] = sphere(100);
props.AmbientStrength = 1;  %full ambient illumination
props.DiffuseStrength = 1;  %full diffuse illumination
props.SpecularColorReflectance = 0.5;
props.SpecularExponent = 1;
props.SpecularStrength = 0.45;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.Cdata = earth;
Earth = surface(R_E*x - zp(1),R_E*y - zp(2), R_E*flip(z) - zp(3),props);
lightangle(gca,-45,30)

hold on
plot3(rvec(1,:)-zp(1), rvec(2,:)-zp(2), rvec(3,:)-zp(3),'g','LineWidth',3)

%plot apses
plot3(rvec(1,1)-zp(1), rvec(2,1)-zp(2), rvec(3,1)-zp(3),...
    'b.','MarkerSize',20)
plot3(rvec(1,51)-zp(1), rvec(2,51)-zp(2), rvec(3,51)-zp(3),...
    'r.','MarkerSize',20)
hold off

axis equal
view(az,el)
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Visible','off')
%axis([-1,1,-1,1,-1,1]*max(abs(axis())))
