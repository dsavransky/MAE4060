% Animation for exploring the ecliptic and equatorial planes.  The
% animation shows the nominal orbit of the Earth (using J2000 orbital
% elements) along with the axes of the IERS ecliptic frame (yellow), the
% axes of the orbit's perifocal frame (red) and the axes of the equatorial
% frame (blue).  The animation pauses at various points of interest in the
% orbit (continue by hitting any key in the figure window).

% Copyright (c) 2018 Dmitry Savransky (ds264@cornell.edu)

nframes = 30*10; %30 fps for 10 s

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
T_p = 2*pi/n; %orbital period in s
R_S = 6.957e8; %Sun radius (m)
R_S = R_S/R_E;

t = linspace(0,T_p,nframes);

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
%%
f1 = figure(1);
clf(f1)

emult = 2500;
smult = 20;

earth = imread('landOcean.jpg');
[x,y,z] = sphere(100);
props.AmbientStrength = 0.12;  %no ambient light
props.DiffuseStrength = 1;  %full diffuse illumination
props.SpecularColorReflectance = 0.5;
props.SpecularExponent = 1;
props.SpecularStrength = 0.45;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.Cdata = earth;
Earth = surface(emult*x,emult*y,emult*flip(z),props);
hold on

props2.AmbientStrength = 1;  %no ambient light
props2.DiffuseStrength = 0;  %full diffuse illumination
props2.FaceColor= 'texture';
props2.EdgeColor = 'none';
props2.FaceLighting = 'flat';
props2.FaceColor = 'y';
Sun = surface(smult*R_S*x,smult*R_S*y,smult*R_S*z,props2);
l = light('Position',[0,0,0],'Style','local');

plot3(r(1,:),r(2,:),r(3,:),'b','Linewidth',2)

eclipax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.5*a,0,0],[0,1.5*a,0],[0,0,1.5*a],0,'Linewidth',2,'color','y');
plot3([0,-a],[0,0],[0,0],'y--','Linewidth',2)
plot3([0,0],[0,-a],[0,0],'y--','Linewidth',2)

h = hgtransform;
perifocalax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.25*a,0,0],[0,1*a,0],[0,0,1.5*a],0,'Linewidth',2,'color','r','Parent',h);
plot3([0,-a],[0,0],[0,0],'r--','Linewidth',2,'Parent',h)
set(h,'Matrix',makehgtform('zrotate',omega));

g = hgtransform;
Earth.Parent = g;
equatax = quiver3([0,0,0],[0,0,0],[0,0,0],emult*[2,0,0],emult*[0,2,0],emult*[0,0,2],0,'Linewidth',2,'color','b','Parent',g);
set(g,'Matrix',makehgtform('translate',r(:,1),'xrotate',-vareps));

view(110,20)
axis equal
ax = gca();
ax.Visible = 'off';
set(gcf,'Color','k');
axis manual

axlim = axis();
txt = text(0,1,0,'Perihelion ~ Jan 3','Units','normalized','Color','w','FontSize',36);

btn = uicontrol('Style','togglebutton','String','Continue','FontSize',24,'Position',[0,0,150,40]);

%%
waitfor(btn,'value');

%find points of interest
poi = {'Vernal Equinox ~ Mar 20','Summer Solistice ~Jun 21','Aphelion ~ Jul 3','Autumnal Equinox ~ Sep 22','Winter Solstice ~ Dec 21', 'Perihelion ~ Jan 3'};
[~,vei] = min(abs(nu + (omega - pi/2) - pi/2));
[~,ssi] = min(abs(nu + (omega - pi/2) - pi));
[~,api] = min(abs(nu - pi));
[~,aei] = min(abs(nu + (omega - pi/2) - 3*pi/2));
[~,wsi] = min(abs(nu + (omega - pi/2) - 2*pi));
ioi = [vei,ssi,api,aei,wsi,length(t)];

j = 1;
for k = 1:length(ioi)
    set(btn,'Visible','off')
    set(txt,'String','')
    for j = j+1:ioi(k)
        %set(Earth,'XData',emult*x+r(1,j),'YData',emult*y+r(2,j),'ZData',emult*z+r(3,j))
        set(g,'Matrix',makehgtform('translate',r(:,j),'xrotate',-vareps));
        pause(1/30)
    end
    set(txt,'String',poi{k});
    if k ~= length(ioi)
        set(btn,'Visible','on');
        waitfor(btn,'value');
    end
end