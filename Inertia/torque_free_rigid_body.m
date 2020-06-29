function torque_free_rigid_body(varargin)

p = inputParser;
addOptional(p,'w0', [1,0.5,0.25], @(x) (numel(x) == 3) && isnumeric(x));
addOptional(p,'I', [1,2,3], @(x) (numel(x) == 3) && isnumeric(x));
addOptional(p,'t', linspace(0,100,1000),@(x) (numel(x) >=2) && ...
                                              isnumeric(x) && issorted(x));
addOptional(p,'angs0', [50,40,30]*pi/180, ...
    @(x) (numel(x) == 3) && isnumeric(x));  
addOptional(p,'poinsot',false, @(x) isscalar(x) && islogical(x));
addOptional(p,'xtrace',true, @(x) isscalar(x) && islogical(x));
addOptional(p,'ytrace',false, @(x) isscalar(x) && islogical(x));
addOptional(p,'ztrace',true, @(x) isscalar(x) && islogical(x));
addOptional(p,'animate',true, @(x) isscalar(x) && islogical(x));

parse(p,varargin{:});

w0 = p.Results.w0(:);
I = p.Results.I;
t = p.Results.t;
h = sqrt(I(1)^2*w0(1)^2 + I(2)^2*w0(2)^2 + I(3)^2*w0(3)^2); %ang momentum
T2 = w0.'*diag(I)*w0; %2*KE

if p.Results.poinsot
    %orient initially such that h = -e_3
    psi0 = 0;
    theta0 = atan2(-sqrt(I(1)^2*w0(1)^2 + I(2)^2*w0(2)^2),-I(3)*w0(3));
    phi0 = atan2(I(1)*w0(1),I(2)*w0(2));
    xtrace = false;
    ytrace = false;
    ztrace = false;
else
    %orient initially based on inputs/defaults
    psi0 = p.Results.angs0(1);
    theta0 = p.Results.angs0(2);
    phi0 = p.Results.angs0(3);
    xtrace = p.Results.xtrace;
    ytrace = p.Results.ytrace;
    ztrace = p.Results.ztrace;
end

omega2dang = @(psi,theta,phi,omega_1,omega_2,omega_3)...
  [(omega_1.*sin(phi) + omega_2.*cos(phi))./sin(theta);...
   omega_1.*cos(phi) - omega_2.*sin(phi);...
 -omega_1.*sin(phi)./tan(theta) - omega_2.*cos(phi)./tan(theta) + omega_3];

dang2omega = @(psi,theta,phi,psidot,thetadot,phidot)...
    [psidot.*sin(phi).*sin(theta) + thetadot.*cos(phi), ...
    psidot.*sin(theta).*cos(phi) - thetadot.*sin(phi), ...
    phidot + psidot.*cos(theta)];

%calculate initial angle derivatives
dangs0 = omega2dang(psi0,theta0,phi0,w0(1),w0(2),w0(3));

%create initial state
z0 = [psi0,theta0,phi0,dangs0.'];

%integrate!
[~,z] = torque_free_motion3(t,z0,I);

%calculate time histories
%body frame omegas:
ws = dang2omega(z(:,1),z(:,2),z(:,3),z(:,4),z(:,5),z(:,6));
ws_I = zeros(size(ws)); %inertial frame omegas
hs_norm = zeros(size(ws)); %body frame ang momentum direction
ws_norm = zeros(size(ws)); %body frame omega direction
bCis = zeros(length(z),3,3); %DCMs
polhode = zeros(size(ws)); %intersection of omega with ellipsoid in B frame
for j = 1:length(z)
    tmp = EulerAngs2DCM(z(j,1:3),[3,1,3]);
    ws_norm(j,:) = ws(j,:)/norm(ws(j,:));
    ws_I(j,:) = (tmp.'*ws_norm(j,:).').';
    hs_norm(j,:) = (diag(I)*ws(j,:).').'/h;
    bCis(j,:,:) = tmp;
    polhode(j,:) = sqrt(T2)*ws(j,:)/sqrt((ws(j,:)*diag(I)*ws(j,:).'));
end

mats = bCis; %hgtransform matrix (inertial frame)
mats(:,4,4) = 1;
bCi0 = squeeze(bCis(1,:,:)); %initial DCM
h_I = bCi0.'*diag(I)*ws(1,:).'/h; %inertia angular momentum 

%inertial frame
figure(1)
clf
set(1,'Position',[0 215 400 300])
g1 = hgtransform;
axscale = 1.75;
xscale = sqrt(T2/I(1))*axscale;
yscale = sqrt(T2/I(2))*axscale;
zscale = sqrt(T2/I(3))*axscale;
oscale  = mean([xscale,yscale,zscale]);
[xe,ye,ze]=ellipsoid(0,0,0,sqrt(T2/I(1)),sqrt(T2/I(2)),...
    sqrt(T2/I(3)),100);
s1 = surface(xe,ye,ze,'FaceAlpha',0.8,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b','Parent',g1);
hold on
xax = quiver3(0,0,0,xscale,0,0,'Linewidth',3,...
    'Color','r','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
yax = quiver3(0,0,0,0,yscale,0,'Linewidth',3,...
    'Color','g','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
zax = quiver3(0,0,0,0,0,zscale,'Linewidth',3,...
    'Color','b','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
omegaax = quiver3(0,0,0,ws_norm(1,1)*oscale,ws_norm(1,2)*oscale,...
    ws_norm(1,3)*oscale,'Linewidth',3,'Color','k','Parent',g1,...
    'AutoScale','off','MaxHeadSize',0.5);
hax = quiver3(0,0,0,h_I(1)*oscale,h_I(2)*oscale,h_I(3)*oscale,...
    'Linewidth',3,'Color',[0.5,0.5,0.5],'AutoScale','off','MaxHeadSize',0.5);
if xtrace
    xaxhist = plot3(bCi0(1,1)*xscale,bCi0(1,2)*xscale,bCi0(1,3)*xscale,'r--');
end
if ytrace
    yaxhist = plot3(bCi0(2,1)*yscale,bCi0(2,2)*yscale,bCi0(2,3)*yscale,'g--');
end
if ztrace
    zaxhist = plot3(bCi0(3,1)*zscale,bCi0(3,2)*zscale,bCi0(3,3)*zscale,'b--');
end
omegaaxhist = plot3(ws_I(1,1)*oscale,ws_I(1,2)*oscale,...
    ws_I(1,3)*oscale,'k--');

%set up lighting
l = light('Position',[0 -1.5 0]);
lighting phong
shading interp

view(3)
axis equal
axmx = max(max(sqrt(T2./I))*axscale,1);
axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])

if p.Results.poinsot
    invplanez = -T2/h; %height of ellipsoid center above invariable plane
    invplane = fill3([-1,-1,1,1,-1]*axmx,[-1,1,1,-1,-1]*axmx,...
        [1,1,1,1,1]*invplanez,'c','FaceAlpha',0.2,'EdgeColor','None');
    
    herpolhode = [invplanez./ws_I(:,3).*ws_I(:,1),...
                  invplanez./ws_I(:,3).*ws_I(:,2),...
                  zeros(length(t),1)+invplanez];
    set(omegaaxhist,'ZData',invplanez,'Linestyle','-','Linewidth',2,...
        'Color',[0.3,0.3,0.3],'XData',herpolhode(1,1),'YData',herpolhode(1,2))
    
    polhodehisti = plot3(polhode(1,1),polhode(1,2),polhode(1,3),'k',...
        'LineWidth',2,'Parent',g1);
end

set(g1,'Matrix',squeeze(mats(1,:,:)).');
set(gca,'Visible','off')
camzoom(gca,2)

%body frame
figure(2)
clf
set(2,'Position',[400 215 400 300])
s2 = surface(xe,ye,ze,'FaceAlpha',0.8,'SpecularExponent',10,...
    'SpecularStrength',0.3,'FaceColor','b');
hold on
xaxb = quiver3(0,0,0,xscale,0,0,'Linewidth',3,...
    'Color','r','AutoScale','off','MaxHeadSize',0.5);
yaxb = quiver3(0,0,0, 0,yscale,0,'Linewidth',3,...
    'Color','g','AutoScale','off','MaxHeadSize',0.5);
zaxb = quiver3(0,0,0,0,0,zscale,'Linewidth',3,...
    'Color','b','AutoScale','off','MaxHeadSize',0.5);
omegaaxb = quiver3(0,0,0,ws_norm(1,1)*oscale,ws_norm(1,2)*oscale,...
    ws_norm(1,3)*oscale,'Linewidth',3,'Color','k',...
    'AutoScale','off','MaxHeadSize',0.5);
polhodehist = plot3(polhode(1,1),polhode(1,2),polhode(1,3),'k','LineWidth',2);
haxb = quiver3(0,0,0,hs_norm(1,1)*oscale,hs_norm(1,2)*oscale,...
    hs_norm(1,3)*oscale,'Linewidth',3,'Color',[0.5,0.5,0.5],...
    'AutoScale','off','MaxHeadSize',0.5);

lb = light('Position',[1.5 1.5 1.5]);
lighting phong
shading interp

view(120,30)
axis equal
axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])
set(g1,'Matrix',squeeze(mats(1,:,:)).');
set(gca,'Visible','off')
camzoom(gca,2.5)

if p.Results.animate
    for j = 2:length(t)
        set(g1,'Matrix',squeeze(mats(j,:,:)).');
        if xtrace
            set(xaxhist,'XData',[get(xaxhist,'XData'),mats(j,1,1)*xscale],...
                'YData',[get(xaxhist,'YData'),mats(j,1,2)*xscale],...
                'ZData',[get(xaxhist,'ZData'),mats(j,1,3)*xscale])
        end
        if ytrace
            set(yaxhist,'XData',[get(yaxhist,'XData'),mats(j,2,1)*yscale],...
                'YData',[get(yaxhist,'YData'),mats(j,2,2)*yscale],...
                'ZData',[get(yaxhist,'ZData'),mats(j,2,3)*yscale])
        end
        if ztrace
            set(zaxhist,'XData',[get(zaxhist,'XData'),mats(j,3,1)*zscale],...
                'YData',[get(zaxhist,'YData'),mats(j,3,2)*zscale],...
                'ZData',[get(zaxhist,'ZData'),mats(j,3,3)*zscale])
        end
        if p.Results.poinsot
            set(omegaaxhist,'XData',herpolhode(1:j,1),...
                'YData',herpolhode(1:j,2),...
                'ZData',herpolhode(1:j,3))
            set(polhodehisti,'XData',polhode(1:j,1),...
                'YData',polhode(1:j,2),...
                'ZData',polhode(1:j,3))
        else
            set(omegaaxhist,'XData',ws_I(1:j,1)*oscale,...
                'YData',ws_I(1:j,2)*oscale,...
                'ZData',ws_I(1:j,3)*oscale)
        end
        set(omegaax,'UData',ws_norm(j,1)*oscale,'VData',ws_norm(j,2)*oscale,...
            'WData',ws_norm(j,3)*oscale)
        set(omegaaxb,'UData',ws_norm(j,1)*oscale,'VData',ws_norm(j,2)*oscale,...
            'WData',ws_norm(j,3)*oscale)
        set(haxb,'UData',hs_norm(j,1)*oscale,'VData',hs_norm(j,2)*oscale,...
            'WData',hs_norm(j,3)*oscale)
        set(polhodehist,'XData',polhode(1:j,1),...
            'YData',polhode(1:j,2),...
            'ZData',polhode(1:j,3))
        
        pause((t(j)-t(j-1)));
    end
else
    set(g1,'Matrix',squeeze(mats(end,:,:)).');
    if xtrace
        set(xaxhist,'XData',mats(:,1,1)*xscale,...
            'YData',mats(:,1,2)*xscale,...
            'ZData',mats(:,1,3)*xscale)
    end
    if ytrace
        set(yaxhist,'XData',mats(:,2,1)*yscale,...
            'YData',mats(:,2,2)*yscale,...
            'ZData',mats(:,2,3)*yscale)
    end
    if ztrace
        set(zaxhist,'XData',mats(:,3,1)*zscale,...
            'YData',mats(:,3,2)*zscale,...
            'ZData',mats(:,3,3)*zscale)
    end
    if p.Results.poinsot
        set(omegaaxhist,'XData',herpolhode(:,1),...
            'YData',herpolhode(:,2),...
            'ZData',herpolhode(:,3))
        set(polhodehisti,'XData',polhode(:,1),...
            'YData',polhode(:,2),...
            'ZData',polhode(:,3))
    else
        set(omegaaxhist,'XData',ws_I(:,1)*oscale,...
            'YData',ws_I(:,2)*oscale,...
            'ZData',ws_I(:,3)*oscale)
    end
    set(omegaax,'UData',ws_norm(end,1)*oscale,'VData',ws_norm(end,2)*oscale,...
        'WData',ws_norm(end,3)*oscale)
    set(omegaaxb,'UData',ws_norm(end,1)*oscale,'VData',ws_norm(end,2)*oscale,...
        'WData',ws_norm(end,3)*oscale)
    set(haxb,'UData',hs_norm(end,1)*oscale,'VData',hs_norm(end,2)*oscale,...
        'WData',hs_norm(end,3)*oscale)
    set(polhodehist,'XData',polhode(:,1),...
        'YData',polhode(:,2),...
        'ZData',polhode(:,3))
    
end
end