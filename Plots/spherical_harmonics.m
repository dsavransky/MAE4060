%%
earth = imread('landOcean.jpg');
[x,y,z] = sphere(100);

rotMats = DCMs();
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
Earth = surface(x,y,flip(z),props);
hold on
oax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.5,0,0],[0,1.5,0],[0,0,1.25],0,'Linewidth',2,'color','r');

sun = light('Position',[1e5,0,0],'Style','local'); %light is effective point source at sun location

view(130,15)
axis equal off

%% sectorals
s = [];
counter = 1;
for j = 0:20:360
    ph = linspace(0,pi,200);
    th = linspace((j-5)*pi/180,(j+5)*pi/180,200);
    [ph,th] = meshgrid(ph,th);
    
    X = cos(th).*sin(ph);
    Y = sin(th).*sin(ph);
    Z = cos(ph);
    
    s(counter) = surface(X,Y,Z,'EdgeColor','none','FaceColor','m','FaceAlpha',0.5);
    counter = counter+1;
end

%% zonals

s = [];
counter = 1;
for j = -170:20:170
    ph = linspace((j-5)*pi/180,(j+5)*pi/180,200);
    th = linspace(0,2*pi,200);
    [ph,th] = meshgrid(ph,th);
    
    X = cos(th).*sin(ph);
    Y = sin(th).*sin(ph);
    Z = cos(ph);
    
    s(counter) = surface(X,Y,Z,'EdgeColor','none','FaceColor','c','FaceAlpha',0.5);
    counter = counter+1;
end