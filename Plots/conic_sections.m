th = linspace(0,2*pi,100);
r= 1;
h = 3;
z = linspace(0,h,100);

[TH,Z] = meshgrid(th,z);
X = (h - Z)/h*r.*cos(TH);
Y = (h - Z)/h*r.*sin(TH);

figure(1)
clf()
props.AmbientStrength = 0.5;  %no ambient light
props.DiffuseStrength = 1;  %full diffuse illumination
props.SpecularColorReflectance = 0.5;
props.SpecularExponent = 1;
props.SpecularStrength = 0.45;
props.FaceColor= 'c';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.FaceAlpha = 0.5;
surface(X,Y,Z,props)
view(110,15)
hold on
plot([-r,r],[0,0],'k')
plot(cos(th),sin(th),'k')
l = light('Position',[5,1,0],'Style','local');
axis equal

p1 = fill3([-1,-1,1,1,-1]*0.25,[-1,1,1,-1,-1]*0.25,[0,0,0,0,0]+h-0.5,'b','FaceAlpha',0.5);

p2 = fill3([-1,-1,1,1,-1]*0.75,[-1,1,1,-1,-1]*0.75,[0,0,0,0,0],'b','FaceAlpha',0.5);
rotate(p2,[1,0,0],25,[0,0,0])
p2.ZData = p2.ZData+h/2;

p3 = fill3([-1,-1,1,1,-1],[-1,1,1,-1,-1]*1.75,[0,0,0,0,0],'b','FaceAlpha',0.5);
rotate(p3,[1,0,0],atan(h)*180/pi,[0,0,0])

p4 = fill3([-1,-1,1,1,-1],[-1,1,1,-1,-1],[0,0,0,0,0],'b','FaceAlpha',0.5);
rotate(p4,[1,0,0],-45,[0,0,0])
p4.YData = p4.YData-0.5;
zlim([0,h])