I = [1 2 3]; %principal moments of inertia

w = [0.5,1.5,0.1]; %angular velocity
%w = [1,0.5,0.25];
n = w/norm(w); %axis of rotation

[t,ws] = torque_free_motion([0,100],w(:),I);

%intersection point of omega and ellipse
xs = 1./sqrt(I(1) + I(2)*(ws(:,2)./ws(:,1)).^2 + I(3)*(ws(:,3)./ws(:,1)).^2);
ys = xs.*ws(:,2)./ws(:,1);
zs = xs.*ws(:,3)./ws(:,1);
if max(diff(zs))/max(diff(xs)) > 10
    xs(zs < 0) = -xs(zs < 0);
    ys(zs < 0) = -ys(zs < 0);
    zs(zs < 0) = -zs(zs < 0);
end

h = diag(I)*w.'; %angular momentum
th = atan2(h(2),h(1));
ph = acos(h(3)/norm(h));

r1 = rotMatE1(-ph);

if ishandle(1), close(1); end
figure(1)
set(1,'Position',[70,100,1200,1200])
clf
g1 = hgtransform;

[xe,ye,ze]=ellipsoid(0,0,0,sqrt(1/I(1)),sqrt(1/I(2)),...
    sqrt(1/I(3)),100);
s1 = surface(xe,ye,ze,'FaceAlpha',0.2,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b','Parent',g1);
hold on
%set up lighting
l(1) = light('Position',[0 -1.5 0]);
%l(2) = light('Position',[1.5 0.5 -0.5]);
lighting phong
shading interp

view(3)
axis equal
axmx = max(max(1./sqrt(I))*1.5,1);
axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])

pc = [0.3,0.3,0.3];
xax = quiver3(0,0,0, sqrt(1/I(1))*1.5,0,0,'Linewidth',1,'Color','k','Parent',g1);
yax = quiver3(0,0,0, 0, sqrt(1/I(2))*1.5,0,'Linewidth',1,'Color','k','Parent',g1);
zax = quiver3(0,0,0,0,0,sqrt(1/I(3))*1.5,'Linewidth',1,'Color','k','Parent',g1);

nax = quiver3(0,0,0,n(1),n(2),n(3),0,'Color','r','LineWidth',2,'Parent',g1);
hax = quiver3(0,0,0,h(1)/norm(h),h(2)/norm(h),h(3)/norm(h),0,'Color','b','LineWidth',2,'Parent',g1);

polhode =  plot3(xs,ys,zs,'y','Linewidth',2,'Parent',g1);

set(g1,'Matrix',makehgtform('yrotate',pi-ph)*makehgtform('zrotate',-th));

tmp = g1.Matrix*[xs(1),ys(1),zs(1),0].';
invplane = fill3([-1,-1,1,1,-1]*axmx,[-1,1,1,-1,-1]*axmx,[1,1,1,1,1]*tmp(3),'c','FaceAlpha',0.2,'EdgeColor','None');

herpolhode = plot3(tmp(1),tmp(2),tmp(3),'k','LineWidth',2);
%%
for j = 2:length(t)
    
    hj = diag(I)*ws(j,:).'; %angular momentum
    thj = atan2(hj(2),hj(1));
    phj = acos(hj(3)/norm(hj));
    nj = ws(j,:)/norm(ws(j,:)); %axis of rotation

    set(nax,'UData',nj(1),'VData',nj(2),'WData',nj(3))
    set(hax,'UData',hj(1)/norm(hj),'VData',hj(2)/norm(hj),'WData',hj(3)/norm(hj))
    set(g1,'Matrix',makehgtform('yrotate',pi-phj)*makehgtform('zrotate',-thj));
    
    set(polhode,'XData',xs(1:j),'YData',ys(1:j),'ZData',zs(1:j))
    
    tmp = g1.Matrix*[xs(j),ys(j),zs(j),0].';
    if tmp(3) > 0
        tmp = -tmp;
    end
    
    set(herpolhode,'XData',[get(herpolhode,'XData'),tmp(1)],...
        'YData',[get(herpolhode,'YData'),tmp(2)],...
        'ZData',[get(herpolhode,'ZData'),tmp(3)])
    
    pause((t(j)-t(j-1))/2);
end
