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

%t = linspace(0,1*86400,2); %solar day

M = 50*pi/180;
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

earth = imread('landOcean.jpg'); %read Earth map (should be in MATLAB demos dir)%%

%%
f1 = figure(1);
clf(f1)

[x,y,z] = sphere(100);
props.AmbientStrength = 0.75;  %some ambient light
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

sun = light('Position',rsun(:,1),'Style','local'); %light is effective point source at sun location



% equatax = plot3([0,1.5],[0,0],[0,0],[0,0],[0,1.05],[0,0],[0,0],[0,0],[0,1.5], 'color','b');
% rotate(equatax,[0,0,1],-40,[0,0,0])
% 
% eclipax = plot3([0,0],[0,1.5],[0,0],[0,0],[0,0],[0,1.5], 'Linewidth',1,'color','r');
% rotate(eclipax,[1,0,0],23.44,[0,0,0])
% rotate(eclipax,[0,0,1],-40,[0,0,0])


axis equal
view(60,10)

th = linspace(0,2*pi,100);
circ = [cos(th);sin(th);zeros(1,length(th))];

equator = plot3(circ(1,:),circ(2,:),circ(3,:),'b','Linewidth',2);
ecliptic = plot3(circ(1,:),circ(2,:),circ(3,:),'r','Linewidth',2);
rotate(ecliptic,[1,0,0],23.44,[0,0,0])
rotate(ecliptic,[0,0,1],-40,[0,0,0])

%% plot prime meridian
pm = plot3(circ(1,:),circ(2,:),circ(3,:),'w','Linewidth',2);
rotate(pm,[1,0,0],90,[0,0,0])

g = hgtransform;

Earth.Parent = g;
equator.Parent = g;
ecliptic.Parent = g;
pm.Parent = g;

set(g,'Matrix',makehgtform('xrotate',-23.44*pi/180));
view(80,25)
%% plot planes

xmn = -2;
xmx = 4.5;
ymn = -2;
ymx = 2;
equatpl = fill3([xmn xmn xmx xmx xmn]*1.5,[ymn ymx ymx ymn ymn]*1.5,[0 0 0 0 0],'b');
set(equatpl,'EdgeColor','none','FaceAlpha',0.3,'FaceLighting','none')
rotate(equatpl,[0,0,1],-40,[0,0,0])

eclippl = fill3([xmn xmn xmx xmx xmn]*1.5,[ymn ymx ymx ymn ymn]*1.5,[0 0 0 0 0],'r');
set(eclippl,'EdgeColor','none','FaceAlpha',0.3,'FaceLighting','none')
rotate(eclippl,[1,0,0],23.44,[0,0,0])
rotate(eclippl,[0,0,1],-40,[0,0,0])

