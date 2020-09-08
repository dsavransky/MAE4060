function dz = sphericalPendulum_eom(~,z)
%integrator function for spherical pendulum equations of motion
%constants are hard-coded

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
