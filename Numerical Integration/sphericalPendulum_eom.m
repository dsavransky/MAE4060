function dz = sphericalPendulum_eom(~,z)
% sphericalPendulum_eom provides the equations of motion for a spherical
% pendulum suitalbe for integration with ode45 or any other MATLAB built-in
% integrator. The pendulum arm length is fixed to 1 m.
%      
% Inputs:
%   z (1x4 float): State vector: phi, phidot, theta, thetadot where phi and
%                  theta are the zenith and azimuth angles, respectively
%                  (and dot represents derivative in time).
% Outputs:
%   dz (1x4 float): State vector derivative. 

g = 9.81; %m/s
l = 1; %m

%z = [phi, phidot, theta, thetadot]
ph = z(1);
phd = z(2);
thd = z(4);

dz = [phd;
    thd^2*sin(ph)*cos(ph) - g/l*sin(ph);
    thd;
    -2*thd*phd/tan(ph)];

end
