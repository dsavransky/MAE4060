function [t,r,v,te,ze,ie] = ballistic_flight_dragEv2()
% ballistic_flight_dragEv2 numerically integrates the equations of motion 
% for a projectile launched from the ground and influenced by gravity and a
% simplified drag force that acts at constant magnitude in the direction 
% oppostie of the projectile's velocity.  Plots are generated of the 
% projectile's 3D trajectory and height above ground as a function of time. 
% A terminal event stops the numerical integration at the time when the
% projectile returns to a zero height. A non-terminal event records all
% times when the projectile passes 5000 m height above ground and marks the
% corresponding points in the plots. 
%      
% Inputs:
%   None
%
% Outputs:
%   t (nx1 float): Array of times
%   r (nx3 float): Array of position values corresponding to each time
%   v (nx3 float): Array of position values corresponding to each time
%   te (mx1 float): Array of event times
%   ze (mx6 float): Array of state values at event times
%   ie (mx1 int): Index of event. 1 is the terminal ground crossing event
%                 and 2 is the non-terminal 5000 m crossing event.

% Copyright (c) 2020 Dmitry Savransky (ds264@cornell.edu)

g = 9.81; %m/s^2
Fdm = 5; %m/s^2

r0 = [0,0,0]; %m
v0 = [100,100,500]; %m/s

z0 = [r0.';v0.'];

[t,z,te,ze,ie] = ode45(@ballistic_flight_drag_eom,linspace(0,150,1000),z0,...
    odeset('Events',@ballistic_flight_drag_event));

r = z(:,1:3);
v = z(:,4:6);
r5k = ze(ie==2,:);
t5k = te(ie==2);

figure(3)
clf
plot(t,r(:,3))
hold on
plot(t5k,r5k(:,3),'.','MarkerSize',20)
hold off
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (s)')
ylabel('Height (m)')
ylim([0,max(r(:,3))])

figure(4)
clf
plot3(r(:,1),r(:,2),r(:,3))
hold on
plot3(r5k(:,1),r5k(:,2),r5k(:,3),'.','MarkerSize',20)
hold off

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


    %0 and 5000 m height crossing event function
    function [pos,isterm,dir] = ballistic_flight_drag_event(~,z)
        h = z(3); %height above ground
        pos = [h,h-5000]; %zero crossing
        isterm = [1,0]; %halt
        dir = [-1,0];    %look for decreasing events only
    end


end