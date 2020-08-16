earth = imread('landOcean.jpg');
[x,y,z] = sphere(100);
rotMats = DCMs();
%%

I = asin(sqrt(4/5));
omega = 90*pi/180;
Omega = 80*pi/180;
e = 0.74;  %eccentricity
J2 = 1082.6e-6;

R_E = 6378.1366; %Earth radius (km)

a = 26600; %km

w_e = 7.2921150e-5; %Earth rotation (rad/s) 
mue = 3.986004418e5; %km^3/s^2

n = sqrt(mue/a^3); %s^(-1)
T_p = 2*pi/n; %orbital period in s

t = linspace(0,2*T_p+60,1000);

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

%%

rrot = zeros(size(r));
for j = 1:length(r)
    rrot(:,j) = rotMats{3}(w_e*t(j))*r(:,j);
end

lons = atan2(rrot(2,:),rrot(1,:));
lats = atan2(sqrt(sum(rrot(1:2,:).^2,1)), rrot(3,:))-pi/2;
%%
f2 = figure(2);
clf
w = worldmap('World');
cdat = load('coast');
plotm(cdat.lat, cdat.long)
hold on
plotm(lats*180/pi,lons*180/pi,'r.','LineWidth',3)
hold off
