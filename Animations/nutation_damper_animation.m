function [T,Z] = nutation_damper_animation(varargin)
%nutation_damper_animation animates the action of a viscous nutation damper
%on a spinning, symmetric spacecraft. The spacecraft begins in a stable
%spin about its symmetry axis, and then (after 5 revolutions) experiences
%an angular impulse producing a spin component about one of the transverse
%axes of magnitude equal to 3/5 of the initial spin rate.  In the case of a
%major axis spinner, the nutation damper will re-stabilize the spin about
%the symmetry axis (a minor axis spinner will tumble). The function
%produces a plot of the time histories of the spacecraft and damper angular
%velocity components, and (optionally) animates the time history.
%
%Optional Parameters:
%   w0 (scalar float) Initial angular speed about the spin axis (defaults
%                     to 5 rad/s)
%   IsIt (scalar float) Ratio between the symmetry and non-symmetry moments
%                       of inertia.  Defined between 0 and 2 for a
%                       cylinder. > 1 produces a major axis spinner
%                       (stable) and > 1 produces a minor axis spinner
%                       (unstable).  Defaults to 1.75
%   IwIt (scalar float) Ratio between moment of inertia of nutation damper
%                       and spacecraft transverse moment of inertia.
%                       Defaults to 0.06.
%   DIt (scalar float) Ratio between D and spacecraft transverse moment of
%                      inertia, where the friction torque mangitude of the
%                      damper is equal to D times its spin rate. Defaults
%                      to 0.5.
%   tend (scalar float) Time to integrate past initial disturbance.
%   animate (scalar bool) Animate the result.  Defaults true.
%   saveAnimation (scalar bool) Save the animation to a movie. Defaults
%                               false.
%   filename (string) Filename for saved animation

% Copyright Dmitry Savransky (ds264@cornell.edu) 2021.

p = inputParser;
addOptional(p,'w0', 5, @(x) isscalar(x) && isnumeric(x));
addOptional(p,'IwIt', 0.06, @(x) isscalar(x) && isnumeric(x));
addOptional(p,'IsIt', 1.75, @(x) isscalar(x) && isnumeric(x) && x < 2 && ...
            x > 0);
addOptional(p,'DIt', 0.5, @(x) isscalar(x) && isnumeric(x));
addOptional(p,'tend', 35, @(x) isscalar(x) && isnumeric(x));
addOptional(p,'animate',true, @(x) isscalar(x) && islogical(x));
addOptional(p,'saveAnimation',false, @(x) isscalar(x) && islogical(x));
addOptional(p,'filename','nutation_damper', @(x) ischar(x));

parse(p,varargin{:});
doWrite = p.Results.saveAnimation;
animate = p.Results.animate;
fname = p.Results.filename;
w0 = p.Results.w0;
tend = p.Results.tend;
IwIt = p.Results.IwIt;
IsIt = p.Results.IsIt;
DIt = p.Results.DIt;
q0 = [0,0,0,1];

[T0,Z0] = ode45(@nutation_damper_eom,linspace(0,2*pi*w0/5,2*pi*w0/5*30),...
    [0,w0,0,0,q0,0].',odeset('AbsTol',1e-9,'RelTol',1e-9));
[T,Z] = ode45(@nutation_damper_eom,linspace(0,tend,tend*30),...
    [3/5*w0,w0,0,0,q0,0].',odeset('AbsTol',1e-9,'RelTol',1e-9));

Z = [Z0;Z];
T = [T0;T+T0(end)];

figure(1)
subplot(3,1,1)
plot(T,Z(:,[1,3]))
grid on
ylabel('radians')
xlabel('Time (s)')
leg1 = legend({'\omega_1','\omega_3'},'FontSize',16,'Location','best');
set(leg1,'Position',[0.8322 0.8566 0.1045 0.1369])
xlim([0,T(end)])
subplot(3,1,2)
plot(T,Z(:,2))
xlim([0,T(end)])
grid on
ylabel('radians')
xlabel('Time (s)')
legend({'\omega_2'},'FontSize',16,'Location','best')
subplot(3,1,3)
plot(T,Z(:,4))
xlim([0,T(end)])
grid on
ylabel('radians')
xlabel('Time (s)')
legend({'\Omega'},'FontSize',16,'Location','best')

qs = Z(:,5:8);
Cs = q2DCMvec(qs);
mats = Cs; %hgtransform matrix (inertial frame)
mats(4,4,:) = 1;

if animate || doWrite
    sz = get(0,'Screensize');
    f = figure(2);
    clf
    set(2,'Position',[0 sz(4)-sz(3)*3/8 sz(3)/2 sz(3)*3/8])
    g1 = hgtransform;
    [xc,yc,zc]=cylinder([1,1],100);
    scl = 1/sqrt(6/IsIt - 3);
    c1 = surface(xc,yc,zc/scl - 1/scl/2,'FaceAlpha',0.8,'SpecularExponent',10,...
        'SpecularStrength',0.3,'FaceColor','b','Parent',g1);
    hold on
    cf1 = fill3(xc(1,:),yc(1,:), zeros(1,length(xc)) + 1/scl/2,'b','Parent',g1);
    cf2 = fill3(xc(1,:),yc(1,:), zeros(1,length(xc)) - 1/scl/2,'b','Parent',g1);
    view(3)
    axis equal
    scl2 = 2;
    c2 = surface(xc/scl2/3,yc/scl2/3-1,zc/scl2/4,'FaceAlpha',0.8,'SpecularExponent',10,...
        'SpecularStrength',0.3,'FaceColor','r','Parent',g1);
    cf3 = fill3(xc(1,:)/scl2/3,yc(1,:)/scl2/3-1, zeros(1,length(xc)) + 1/scl2/4,'r','Parent',g1);
    rotate([c2,cf3],[1,0,0],90,[0,-1,0])
    rotate([c1,cf1,cf2,c2,cf3],[1,0,0],90,[0,0,0])
    xax = quiver3(0,0,0,-1.5,0,0,'Linewidth',3,...
        'Color','r','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
    yax = quiver3(0,0,0,0,-1.5,0,'Linewidth',3,...
        'Color','g','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
    zax = quiver3(0,0,0,0,0,-1.5,'Linewidth',3,...
        'Color','b','Parent',g1,'AutoScale','off','MaxHeadSize',0.5);
    
    l = light('Position',[-1 -1 0]);
    lighting phong
    %shading interp
    axmx = 1.5;
    axis([-axmx,axmx,-axmx,axmx,-axmx,axmx])
    set(gca,'Visible','off')
    
    if doWrite
        set(gca,'nextplot','replacechildren');
        set(f,'Visible','off','Renderer','zbuffer')
        vidObj = VideoWriter(fname);
        vidObj.Quality = 100;
        vidObj.FrameRate = 60;
        open(vidObj);
    end
    
    for j = 2:length(T)
        set(g1,'Matrix',mats(:,:,j).');
        if doWrite
            writeVideo(vidObj,getframe(f));
        else
            pause((T(j)-T(j-1))/2);
        end
    end
    
    if doWrite
        close(vidObj);
        set(f,'Visible','on')
    end
end

    function dz = nutation_damper_eom(~,z)
        
        %z = [omega_1, omega_2, omega_3, Omega, q, alpha];
        omega_1 = z(1);
        omega_2 = z(2);
        omega_3 = z(3);
        Omega = z(4);
        
        IwB_B = z(1:3);
        dq = qXi(z(5:8))*IwB_B/2;
        
        dz = [omega_2.*(IsIt.*omega_3 - IwIt.*Omega - omega_3);...
            (IwIt/IsIt).*Omega.*omega_1;...
            (DIt.*Omega - omega_1.*omega_2.*(IsIt - 1))./(1 - IwIt);...
            (-DIt.*Omega./IwIt + omega_1.*omega_2.*(IsIt - 1))./(1 - IwIt);...
            dq; Omega];
    end


end