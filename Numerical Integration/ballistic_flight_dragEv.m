function [t,r,v,te] = ballistic_flight_dragEv()
% ballistic_flight_dragEv numerically integrates the equations of motion 
% for a projectile launched from the ground and influenced by gravity and a
% simplified drag force that acts at constant magnitude in the direction 
% oppostie of the projectile's velocity.  Plots are generated of the 
% projectile's 3D trajectory and height above ground as a function of time. 
% A terminal event stops the numerical integration at the time when the
% projectile returns to a zero height.
%      
% Inputs:
%   None
%
% Outputs:
%   t (nx1 float): Array of times
%   r (nx3 float): Array of position values corresponding to each time
%   v (nx3 float): Array of position values corresponding to each time
%   te (scalar): Time of terminal event.

% Copyright (c) 2020 Dmitry Savransky (ds264@cornell.edu)

g = 9.81; %m/s^2
Fdm = 5; %m/s^2

r0 = [0,0,0]; %m
v0 = [100,100,500]; %m/s

z0 = [r0.';v0.'];

[t,z,te] = ode45(@ballistic_flight_drag_eom,linspace(0,150,1000),z0,...
    odeset('Events',@ballistic_flight_drag_event));

r = z(:,1:3);
v = z(:,4:6);
sz = get(0,'Screensize');
figure(3)
set(3,'Position',[0 sz(4)-sz(3)*3/8 sz(3)/2 sz(3)*3/8])
clf
plot(t,r(:,3))
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (s)')
ylabel('Height (m)')
ylim([0,max(r(:,3))])

figure(4)
set(4,'Position',[sz(3)/2 sz(4)-sz(3)*3/8 sz(3)/2 sz(3)*3/8])
clf
plot3(r(:,1),r(:,2),r(:,3))

    function dz = ballistic_flight_drag_eom(~,z)
        %state vector is z = [x,y,z,dx,dy,dz]
        
        r = z(1:3); %position vector (m)
        v = z(4:6); %velocity vector (m/s)
        vmag = norm(v); %m/s
        vhat = v/vmag;
        
        %acceleration is:
        a = [0;0;-g] - Fdm*vhat;
        dz = [v; a];
    end


    %0 height crossing event function
    function [pos,isterm,dir] = ballistic_flight_drag_event(~,z)
        h = z(3); %height above ground
        pos = h; %zero crossing
        isterm = 1; %halt
        dir = -1;    %look for decreasing events only
    end


end