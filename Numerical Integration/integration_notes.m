% This script demonstrates several approaches to MATLAB numerical
% integration using the equations of motion for a spherical pendulum.
% For the derivation of the equations of motion, see the
% spherical_pendulum_deriv.mlx live scrip in the Notebooks folder. 

% Copyright (c) 2020 Dmitry Savransky (ds264@cornell.edu)

%% define some constants
g = 9.81; %m/s
l = 1; %m

%this is our integrator function (equations of motion of a spherical
%pendulum).  phi and theta are the The state variable is defined as:
%z = [phi, phidot, theta, thetadot].'
dz1 = @(t,z) [z(2);
            z(4)^2*sin(z(1))*cos(z(1)) - g/l*sin(z(1));
            z(4);
            -2*z(4)*z(2)/tan(z(1))];
%we also define the exact same thing in helper function
%sphericalPendulum_eom
 
%initial conditions (must match state order set in integrator function
z0 = [30*pi/180,0,0,1];   

%integration time span options
%this will allow the integrator to take whatever steps it
%wants:
%tspan = [0,30]; 
%this will still allow the integrator to take whatever steps it wants, but
%the ouptut will be exactly at the times specified:
tspan = linspace(0,30,1000);

%these two should produce exactly the same results:
[t1,res1] = ode45(@sphericalPendulum_eom,tspan,z0);
[t2,res2] = ode45(dz1,tspan,z0);

sum(res1 - res2) %all zeros!

%%
%another option - leave g and l as parameters, and wrap two nested lambda
%functions
dz2 = @(t,z,g,l) [z(2);
            z(4)^2*sin(z(1))*cos(z(1)) - g/l*sin(z(1));
            z(4);
            -2*z(4)*z(2)/tan(z(1))];
dz3 = @(t,z) dz2(t,z,9.81,1);
[t3,res3] = ode45(dz3,tspan,z0);

%or, write a helper function with a nest integrator subfunction:
[t3a,res3a] = sphericalPendulum(l,z0);

%again the same result
sum(res3 - res1) %all zero again!
sum(res3a - res1) %all zero again!


%%
%now with a different integrator
[t4,res4] = ode113(dz1,tspan,z0);

sum(res4 - res1) %not zero!

mean(abs( (res4(2:end,:)-res3(2:end,:))./res3(2:end,:) )*100) %percent error, pretty high actually


%%
%odeset() - to set options

[t6,res6] = ode45(dz1,tspan,z0,odeset('RelTol',1e-12,'AbsTol',1e-12));
[t7,res7] = ode45(dz1,tspan,z0,odeset('RelTol',1e-12));
[t8,res8] = ode45(dz1,tspan,z0,odeset('AbsTol',1e-12));

%%
En = @(ph,phd,thd,l,g) -g*l*cos(ph) + l^2*phd.^2/2 + l^2*thd.^2.*sin(ph).^2/2;

En1 = En(res1(:,1),res1(:,2),res1(:,4),l,g); %defaults ode45
En4 = En(res4(:,1),res4(:,2),res4(:,4),l,g); %defaults ode113
En6 = En(res6(:,1),res6(:,2),res6(:,4),l,g); %reltol = abstol = 1e-12
En7 = En(res7(:,1),res7(:,2),res7(:,4),l,g); %reltol=1e-12, default abstol
En8 = En(res8(:,1),res8(:,2),res8(:,4),l,g); %abstol=1e-12, defatul reltol

plot(tspan,En1,tspan,En4,tspan,En6,'--',tspan,En7,tspan,En8,'LineWidth',2)
xlabel('Time (s)')
ylabel('Energy (N m)')
set(gca,'FontName','Times','FontSize',16)
legend({'ODE45 Default', 'ODE113 Default', 'RelTol=AbsTol=1e-12',...
    'RelTol=1e-12, AbsTol=1e-6', 'RelTol=1e-3, AbsTol=1e-12'},...
    'Location','SouthWest')

%% limit to 2nd best energy conservation
ylim([min(En7),max(En7)])

%% limit to best energy conservation
ylim([min(En6),max(En6)])