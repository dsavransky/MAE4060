%Plot the basic geometry of a long. of ascending node change on a circular orbit

%Copyright (c) 2019 Dmitry Savransky (ds264@cornell.edu)

%%
%circular orbit
th = linspace(0,2*pi,1000);
r = [cos(th);sin(th);zeros(size(th))];
rotMats = DCMs();

%orient orbit
O = 30*pi/180;
I = 40*pi/180;
r1 = (rotMats{1}(I)*rotMats{3}(O)).'*r;
h1 = (rotMats{1}(I)*rotMats{3}(O)).'*[0;0;1];

%now double Omega
r2 = (rotMats{1}(I)*rotMats{3}(O*2)).'*r;
h2 = (rotMats{1}(I)*rotMats{3}(O*2)).'*[0;0;1];

%now change both Omega and I
r3 = (rotMats{1}(I*1.5)*rotMats{3}(O*2)).'*r;
h3 = (rotMats{1}(I*1.5)*rotMats{3}(O*2)).'*[0;0;1];


figure(1)
clf

plot3(r1(1,:),r1(2,:),r1(3,:),'r')

hold on
%plot3(r2(1,:),r2(2,:),r2(3,:),'g')
plot3(r3(1,:),r3(2,:),r3(3,:),'g')

oax = quiver3([0,0,0],[0,0,0],[0,0,0],[1.5,0,0],[0,1.5,0],[0,0,1.5],0,'Linewidth',2,'color','b');
nax = quiver3(0,0,0,1.5*cos(O),1.5*sin(O),0,'Linewidth',2,'color','r');
nax2 = quiver3(0,0,0,1.5*cos(O*2),1.5*sin(O*2),0,'--','Linewidth',2,'color','k');
hax1 = quiver3(0,0,0,1.5*h1(1),1.5*h1(2),1.5*h1(3),'Linewidth',2,'color','r');
%hax2 = quiver3(0,0,0,1.5*h2(1),1.5*h2(2),1.5*h2(3),'Linewidth',2,'color','g');
hax3 = quiver3(0,0,0,1.5*h3(1),1.5*h3(2),1.5*h3(3),'Linewidth',2,'color','g');


view(100,15)
axis equal

%%
dists = zeros(1000);
for j = 1:1000
    for k = 1:1000
        dists(j,k) = norm(r1(:,j) - r2(:,k));
    end
end
[v,i] = min(dists);
[vv,j] = min(v);
plot3(r2(1,j),r2(2,j),r2(3,j),'.','MarkerSize',20)
dists(i(j),j) = 1e6;
[v,i] = min(dists);
[vv,j] = min(v);
plot3(r2(1,j),r2(2,j),r2(3,j),'.','MarkerSize',20)

%%
dists = zeros(1000);
for j = 1:1000
    for k = 1:1000
        dists(j,k) = norm(r1(:,j) - r3(:,k));
    end
end
[v,i] = min(dists);
[vv,j] = min(v);
plot3(r3(1,j),r3(2,j),r3(3,j),'.','MarkerSize',20)
dists(i(j),j) = 1e6;
[v,i] = min(dists);
[vv,j] = min(v);
plot3(r3(1,j),r3(2,j),r3(3,j),'.','MarkerSize',20)
