%script for generating correct orientations for orbit Euler angles figure

M = linspace(0,2*pi,1000);
e = 0.3;
a = 1;

[E, nu] = invKepler(M,e);
r = a*(1-e^2)./(1+ e*cos(nu));

%Cartesian coordinates in the perifocal frame
x = r.*cos(nu);
y = r.*sin(nu);
z = zeros(size(x));

rotMatE3 = @(ang) [cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
rotMatE2 = @(ang) [cos(ang) 0 -sin(ang);0 1 0;sin(ang) 0 cos(ang)];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];
rotMats = {rotMatE1,rotMatE2,rotMatE3};
a = eye(3).';

axl = {'e_1','e_2','e_3'};
mul = 1.5;

%%
figure(1)
clf
s1 = plot3(x,y,z);
hold on
ax1 = zeros(3,1);
for k=1:3
   ax1(k) = plot3([0,a(k,1)]*mul,[0,a(k,2)]*mul,[0,a(k,3)]*mul,'k');
   %text(a(k,1)*mul,a(k,2)*mul,a(k,3)*mul,axl{k});
end

axis([-1,1,-1,1,-1,1]*mul)
view(150,15)
set(gca,'Visible','off')
%%

O = 45;
rotate(s1,a(3,:),O,[0,0,0])
a = rotMats{3}(O*pi/180)*a;

for k=1:3
    set(ax1(k),'XData',[0,a(k,1)]*mul,'YData',[0,a(k,2)]*mul,'ZData',[0,a(k,3)]*mul);
end

%%

I = 30;
rotate(s1,a(1,:),I,[0,0,0])
a = rotMats{1}(I*pi/180)*a;

for k=1:3
    set(ax1(k),'XData',[0,a(k,1)]*mul,'YData',[0,a(k,2)]*mul,'ZData',[0,a(k,3)]*mul);
end

%%

w = 100;
rotate(s1,a(3,:),w,[0,0,0])
a = rotMats{3}(w*pi/180)*a;

for k=1:3
    set(ax1(k),'XData',[0,a(k,1)]*mul,'YData',[0,a(k,2)]*mul,'ZData',[0,a(k,3)]*mul);
end