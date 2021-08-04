function [t,res] = sphericalPendulum(l,z0)
% sphericalPendulum numerically integrates the equations of motion for a
% spherica pendulum over a fixed time span.
%      
% Inputs:
%   l (float): Pendulum arm length (m)
%   z0 (1x4 float): Initial conditions in order: phi, phidot, theta,
%                   thetadot where phi and theta are the zenith and azimuth
%                   angles, respectively (and dot represents derivative in
%                   time).
% Outputs:
%   t (1000x1 float): Array of times
%   res (1000x4 float): Array of state values corresponding to each time

% Copyright (c) 2020 Dmitry Savransky (ds264@cornell.edu)

%constants
g = 9.81; %m/s

[t,res] = ode45(@sphericalPendulum_eq,linspace(0,30,1000),z0);

    function dz = sphericalPendulum_eq(~,z)
        
        %z = [phi, phidot, theta, thetadot]
        ph = z(1);
        phd = z(2);
        thd = z(4);
        
        dz = [phd;
            thd^2*sin(ph)*cos(ph) - g/l*sin(ph);
            thd;
            -2*thd*phd/tan(ph)];
    end

end
