f1 = figure(1);
clf(f1)

[x,y,z] = sphere(100);
props.AmbientStrength = 0.5;  %some ambient light
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
eax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.5,0,0],[0,1.5,0],[0,0,1.5],0,'Linewidth',2,'color','b','MaxHeadSize',0.1);

%axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
axis equal off
view(-10,15)

sun = light('Position',[1e5,-1e5,0],'Style','local'); 

%%
lon = -80.575 *pi/180;
lat = 28.533*pi/180;
th = lon;
ph = pi/2 - lat;

ls = [cos(th).*sin(ph);sin(th).*sin(ph);cos(ph)]*1.01;
plot3(ls(1),ls(2),ls(3),'k.','MarkerSize',25)

%%
ths = linspace(0,2*pi,100);
circ = [cos(ths);sin(ths);zeros(1,length(ths))];
equat = plot3(circ(1,:),circ(2,:),circ(3,:),'w','Linewidth',1);
pm = plot3(circ(1,:),circ(2,:),circ(3,:),'w','Linewidth',1);
rotate(pm,[1,0,0],90,[0,0,0])
ll = plot3(circ(1,:),circ(2,:),circ(3,:),'w','Linewidth',1);
rotate(ll,[1,0,0],90,[0,0,0])
rotate(ll,[0,0,1],lon*180/pi,[0,0,0])

%%
I = lat;
beta = asin(cos(I)/cos(lat));
lonu = acos(cos(beta)/sin(I));

orb = plot3(circ(1,:),circ(2,:),circ(3,:),'r','Linewidth',2);
rotate(orb,[1,0,0],lat*180/pi,[0,0,0])
rotate(orb,[0,0,1],(lon-lonu)*180/pi,[0,0,0])
