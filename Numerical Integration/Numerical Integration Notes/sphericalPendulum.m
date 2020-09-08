function [t,res] = sphericalPendulum(l,z0)
%numerical integration of spherical pendulum allowing for variable length
%as input

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
