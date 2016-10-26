rotMatE3 = @(ang) [cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
rotMatE2 = @(ang) [cos(ang) 0 -sin(ang);0 1 0;sin(ang) 0 cos(ang)];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];
rotMats = {rotMatE1,rotMatE2,rotMatE3};

I = [1 2 3]; %principal moments of inertia
T = 1;

%w = [2,3,4]; %angular velocity
w = [1,0.5,0.25];
n = w/norm(w); %axis of rotation

[t,ws] = torque_free_motion([0,100],w(:),I);

h = diag(I)*w.'; %angular momentum
th = atan2(h(2),h(1));
ph = acos(h(3)/norm(h));

r1 = rotMatE1(-ph);

close(1)
figure(1)
clf
g = hgtransform;

[xe,ye,ze]=ellipsoid(0,0,0,sqrt(1/I(1)),sqrt(1/I(2)),...
    sqrt(1/I(3)),100);
s1 = surface(xe,ye,ze,'FaceAlpha',0.2,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b','Parent',g);

hold on
%set up lighting
l(1) = light('Position',[0 -1.5 0]);
%l(2) = light('Position',[1.5 0.5 -0.5]);
lighting phong
shading interp
% cs = colormap('gray');
% cs = cs(30:50,:);
%colormap(cs)
view(3)
axis equal

pc = [0.3,0.3,0.3];
xax = quiver3(0,0,0, sqrt(1/I(1))*1.5,0,0,'Linewidth',2,'Color',pc,'Parent',g);
yax = quiver3(0,0,0, 0, sqrt(1/I(2))*1.5,0,'Linewidth',2,'Color',pc,'Parent',g);
zax = quiver3(0,0,0,0,0,sqrt(1/I(3))*1.5,'Linewidth',2,'Color',pc,'Parent',g);

nax = quiver3(0,0,0,n(1),n(2),n(3),0,'Color','r','LineWidth',2,'Parent',g);
hax = quiver3(0,0,0,h(1)/norm(h),h(2)/norm(h),h(3)/norm(h),0,'Color','b','LineWidth',2,'Parent',g);

%intersection point:
x = 1/sqrt(I(1) + I(2)*(w(2)/w(1))^2 + I(3)*(w(3)/w(1))^2);
y = x*w(2)/w(1);
z = x*w(3)/w(1);

set(g,'Matrix',makehgtform('yrotate',pi-ph)*makehgtform('zrotate',-th));

%plot3(x,y,z,'k.','MarkerSize',30)
% 
% allshapes = [s1,xax,yax,zax,nax,hax];
% rotate(allshapes,[1 0 0],-ph*180/pi,[0 0 0])