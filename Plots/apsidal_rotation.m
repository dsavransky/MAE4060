% Generate a visualization of apsidal rotation.  Based on Vallado (2013)

% Copyright (c) 2019 Dmitry Savransky (ds264@cornell.edu)

earth = imread('landOcean.jpg');
[x,y,z] = sphere(100);
%%

I = pi/2;
omega0 = -pi/2;
Omega = 0;
e = 0.1;  %eccentricity
J2 = 1082.6e-6;

R_E = 6378.1366; %Earth radius (km)

a = R_E+1000;

w_e = 7.2921150e-5; %Earth rotation (rad/s) 
mue = 3.986004418e5; %km^3/s^2

n = sqrt(mue/a^3); %s^(-1)
T_p = 2*pi/n; %orbital period in s

t = linspace(0,10*T_p,10000);

omegadot = 3/2*J2*n*(R_E/a/(1-e^2))^2*(2 - 5/2*sin(I)^2);
omega = omega0 + 100*omegadot*t; %going to exaggerate rotation by factor of 100

%%
M = mod(n*t,2*pi);
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

r = A.*repmat((cos(E) - e),3,1) + B.*repmat(sin(E),3,1); %orbit

%%
f1 = figure(1);
clf(f1)

props.AmbientStrength = 1;  %full ambient light
props.DiffuseStrength = 1;  %full diffuse illumination
props.SpecularColorReflectance = 0.5;
props.SpecularExponent = 1;
props.SpecularStrength = 0.45;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.Cdata = earth;
Earth = surface(x*R_E,y*R_E,flip(z)*R_E,props);
hold on
plot3(r(1,:),r(2,:),r(3,:),'b','LineWidth',2)
plot3(r(1,1),r(2,1),r(3,1),'m.','MarkerSize',20)
plot3(r(1,end),r(2,end),r(3,end),'r.','MarkerSize',20)


sun = light('Position',[-1e5,-1e5,0],'Style','local'); %light is effective point source at sun location

view(0,0)
axis equal off
hold off

