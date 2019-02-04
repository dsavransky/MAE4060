bCi313 = @(psi,theta,phi) [-sin(phi).*sin(psi).*cos(theta) + cos(phi).*cos(psi) sin(phi).*cos(psi).*cos(theta) + sin(psi).*cos(phi) sin(phi).*sin(theta); -sin(phi).*cos(psi) - sin(psi).*cos(phi).*cos(theta) -sin(phi).*sin(psi) + cos(phi).*cos(psi).*cos(theta) sin(theta).*cos(phi); sin(psi).*sin(theta) -sin(theta).*cos(psi) cos(theta)];

omega2dang = @(psi,theta,phi,omega_1,omega_2,omega_3) [(omega_1.*sin(phi) + omega_2.*cos(phi))./sin(theta);omega_1.*cos(phi) - omega_2.*sin(phi);-omega_1.*sin(phi)./tan(theta) - omega_2.*cos(phi)./tan(theta) + omega_3];

dang2omega = @(psi,theta,phi,psidot,thetadot,phidot) [psidot.*sin(phi).*sin(theta) + thetadot.*cos(phi), psidot.*sin(theta).*cos(phi) - thetadot.*sin(phi), phidot + psidot.*cos(theta)];


%%
I = [1 2 3]; %principal moments of inertia
  w0 = [1,0.5,0.25];
% w0 = [1,0,0];
% w0 = [0,0,1];
% w0 = [0,1,0];
% w0 = [0,2,0]+1e-3;
t = linspace(0,100,1000);

%%
[to,ws] = torque_free_motion(t,w0(:),I);

%%
% phi0 = 30*pi/180;
% theta0 = 40*pi/180;
% psi0 = 50*pi/180;

h = sqrt(I(1)^2*w0(1)^2 + I(2)^2*w0(2)^2 + I(3)^2*w0(3)^2);
a = I(1)*w0(1) + I(2)*w0(2);
psi0 = 0;
theta0 = atan2(-sqrt(I(1)^2*w0(1)^2 + I(2)^2*w0(2)^2),-I(3)*w0(3));
phi0 = atan2(I(1)*w0(1),I(2)*w0(2));
bCi0 = EulerAngs2DCM([psi0,theta0,phi0],[3,1,3]);
bCi0.'*diag(I)*w0(:)


dangs0 = omega2dang(psi0,theta0,phi0,w0(1),w0(2),w0(3));

z0 = [psi0,theta0,phi0,dangs0.'];

[~,z2] = torque_free_motion3(t,z0,I);

ws2 = dang2omega(z2(:,1),z2(:,2),z2(:,3),z2(:,4),z2(:,5),z2(:,6));
ws_I = zeros(size(ws2));
hs_I = zeros(size(ws2));
ws_norm = zeros(size(ws2));
bCis = zeros(length(z2),3,3);
polhode = zeros(size(ws2));
for j = 1:length(z2)
    tmp = EulerAngs2DCM(z2(j,1:3),[3,1,3]);
    ws_norm(j,:) = ws2(j,:)/norm(ws2(j,:));
    ws_I(j,:) = (tmp.'*ws_norm(j,:).').';
    hs_I(j,:) = (tmp.'*diag(I)*ws2(j,:).').';
    bCis(j,:,:) = tmp;
    polhode(j,:) = ws2(j,:)/sqrt((ws2(j,:)*diag(I)*ws2(j,:).'));
end

mats = bCis;
mats(:,4,4) = 1;
bCa0 = squeeze(bCis(1,:,:));
 %% just a sanity check - not needed
% bCa0 = squeeze(bCis(1,:,:));
% z0 = [w0(:);bCa0(:)];
% [~,z] = torque_free_motion2(t,z0,I);
% 
% bCas = reshape(z(:,4:end),length(z),3,3);
% 
% dets = zeros(length(bCas),1);
% angs = zeros(length(bCas),3);
% for j = 1:length(bCas)
%     dets(j) = det(squeeze(bCas(j,:,:)));
%     angs(j,:) = DCM2EulerAngs(squeeze(bCas(j,:,:)),[3,1,3]);
% end
% 
% % mats = bCas;
% % mats(:,4,4) = 1;

%%
figure(1)
clf
g1 = hgtransform;
axscale = 1.75;
xscale = sqrt(1/I(1))*axscale;
zscale = sqrt(1/I(3))*axscale;
oscale  = 1.1;
[xe,ye,ze]=ellipsoid(0,0,0,sqrt(1/I(1)),sqrt(1/I(2)),...
    sqrt(1/I(3)),100);
s1 = surface(xe,ye,ze,'FaceAlpha',0.8,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b','Parent',g1);
hold on
xax = quiver3(0,0,0, sqrt(1/I(1))*axscale,0,0,'Linewidth',3,...
    'Color','r','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
yax = quiver3(0,0,0, 0, sqrt(1/I(2))*axscale,0,'Linewidth',3,...
    'Color','g','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
zax = quiver3(0,0,0,0,0,sqrt(1/I(3))*axscale,'Linewidth',3,...
    'Color','b','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
omegaax = quiver3(0,0,0,ws_norm(1,1)*oscale,ws_norm(1,2)*oscale,...
    ws_norm(1,3)*oscale,'Linewidth',3,'Color','k','Parent',g1,...
    'AutoScale','off','MaxHeadSize',0.5);
hax = quiver3(0,0,0,hs_I(1,1),hs_I(1,2),hs_I(1,3),'Linewidth',3,...
    'Color',[0.5,0.5,0.5],'AutoScale','off','MaxHeadSize',0.5); %inertially fixed
xaxhist = plot3(bCa0(1,1)*xscale,bCa0(1,2)*xscale,bCa0(1,3)*xscale,'r--');
zaxhist = plot3(bCa0(3,1)*zscale,bCa0(3,2)*zscale,bCa0(3,3)*zscale,'b--');
omegaaxhist = plot3(ws_I(1,1)*oscale,ws_I(1,2)*oscale,ws_I(1,3)*oscale,'k--');

%set up lighting
l = light('Position',[0 -1.5 0]);
lighting phong
shading interp

view(3)
axis equal
axmx = max(max(1./sqrt(I))*axscale,1);
axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])
tmp = squeeze(mats(1,1:3,1:3)).'*polhode(1,:).';
invplane = fill3([-1,-1,1,1,-1]*axmx,[-1,1,1,-1,-1]*axmx,...
    [1,1,1,1,1]*tmp(3),'c','FaceAlpha',0.2,'EdgeColor','None');


set(g1,'Matrix',squeeze(mats(1,:,:)).');
set(gca,'Visible','off')
camzoom(gca,2)

figure(2)
clf
s1 = surface(xe,ye,ze,'FaceAlpha',0.8,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b');
hold on
xaxb = quiver3(0,0,0, sqrt(1/I(1))*axscale,0,0,'Linewidth',3,...
    'Color','r','AutoScale','off','MaxHeadSize',0.5);
yaxb = quiver3(0,0,0, 0, sqrt(1/I(2))*axscale,0,'Linewidth',3,...
    'Color','g','AutoScale','off','MaxHeadSize',0.5);
zaxb = quiver3(0,0,0,0,0,sqrt(1/I(3))*axscale,'Linewidth',3,...
    'Color','b','AutoScale','off','MaxHeadSize',0.5);
omegaaxb = quiver3(0,0,0,ws_norm(1,1)*axscale,ws_norm(1,2)*axscale,...
    ws_norm(1,3)*axscale,'Linewidth',3,'Color','k',...
    'AutoScale','off','MaxHeadSize',0.5);
polhodehist = plot3(polhode(1,1),polhode(1,2),polhode(1,3),'k','LineWidth',2);

lb = light('Position',[1.5 1.5 1.5]);
lighting phong
shading interp

view(120,30)
axis equal
axmx = max(max(1./sqrt(I))*axscale,1);
axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])
set(g1,'Matrix',squeeze(mats(1,:,:)).');
set(gca,'Visible','off')
camzoom(gca,2.5)
%%
for j = 2:length(to)
    set(g1,'Matrix',squeeze(mats(j,:,:)).');
    set(xaxhist,'XData',[get(xaxhist,'XData'),mats(j,1,1)*xscale],...
                'YData',[get(xaxhist,'YData'),mats(j,1,2)*xscale],...
                'ZData',[get(xaxhist,'ZData'),mats(j,1,3)*xscale])
    set(zaxhist,'XData',[get(zaxhist,'XData'),mats(j,3,1)*zscale],...
                'YData',[get(zaxhist,'YData'),mats(j,3,2)*zscale],...
                'ZData',[get(zaxhist,'ZData'),mats(j,3,3)*zscale])
    set(omegaaxhist,'XData',ws_I(1:j,1)*oscale,...
                    'YData',ws_I(1:j,2)*oscale,...
                    'ZData',ws_I(1:j,3)*oscale)
    set(omegaax,'UData',ws_norm(j,1)*oscale,'VData',ws_norm(j,2)*oscale,...
                'WData',ws_norm(j,3)*oscale)
    set(omegaaxb,'UData',ws_norm(j,1)*axscale,'VData',ws_norm(j,2)*axscale,...
                'WData',ws_norm(j,3)*axscale)
    set(polhodehist,'XData',polhode(1:j,1),...
                    'YData',polhode(1:j,2),...
                    'ZData',polhode(1:j,3))

    pause((to(j)-to(j-1)));
end