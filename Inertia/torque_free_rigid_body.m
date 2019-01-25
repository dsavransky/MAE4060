I = [1 2 3]; %principal moments of inertia
w = [1,0.5,0.25];
t = linspace(0,100,1000);

[to,ws] = torque_free_motion(t,w(:),I);

aCb0 = eye(3);
z0 = [w(:);aCb0(:)];
[to2,z] = torque_free_motion2(t,z0,I);

aCbs = reshape(z(:,4:end),length(z),3,3);
aCbs(:,4,4) = 1;
%%

figure(1)
clf
g1 = hgtransform;
[xe,ye,ze]=ellipsoid(0,0,0,sqrt(1/I(1)),sqrt(1/I(2)),...
    sqrt(1/I(3)),100);
s1 = surface(xe,ye,ze,'FaceAlpha',0.75,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b','Parent',g1);
hold on
%set up lighting
l = light('Position',[0 -1.5 0]);
lighting phong
shading interp

view(3)
axis equal
axmx = max(max(1./sqrt(I))*1.5,1);
axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])

%%
for j = 2:length(to)
    set(g1,'Matrix',squeeze(aCbs(j,:,:)));
    pause((to(j)-to(j-1)));
end