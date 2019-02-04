function [t,z] = torque_free_motion3(t0,z0,I)
%TORQUE_FREE_MOTION2 integrates Euler's equations in the absence of 
%external moments
%
%   [t,w] = torque_free_motion(t0,w0,I) returns the time and angular
%   velocities as a function of time for the torque free motion of a rigid
%   body with principal axis moments of inertia given by the elements of I
%   (a 3 element array) over a time span t0 (a monotonically increasing
%   array of at least 2 elements) with initial angular velocities w0 (a 3  
%   element array).
%
%SEE ALSO: ode45

% Copyright (c) 2017 Dmitry Savransky (ds264@cornell.edu)
    I_1 = I(1);
    I_2 = I(2);
    I_3 = I(3);

    %state  - [psi,theta,phi,psidot,thetadot,phidot]
    [t,z] = ode113(@eulers_eqs,t0,z0,odeset('RelTol',eps(1)*100,'AbsTol',eps(1)*100));

    function dz = eulers_eqs(~,z)
        psi = z(1);
        theta = z(2);
        phi = z(3);
        psidot = z(4);
        thetadot = z(5);
        phidot = z(6);
        
        psiddot = (-I_1.^2.*phidot.*psidot.*cos(2*phi - theta)/4 + I_1.^2.*phidot.*psidot.*cos(2*phi + theta)/4 - I_1.^2.*phidot.*thetadot.*cos(2*phi)/2 - I_1.^2.*phidot.*thetadot/2 - I_1.^2.*psidot.^2.*cos(2*phi - 2*theta)/8 + I_1.^2.*psidot.^2.*cos(2*phi + 2*theta)/8 - I_1.^2.*psidot.*thetadot.*cos(theta)/2 - I_1.^2.*psidot.*thetadot.*cos(2*phi - theta)/4 - I_1.^2.*psidot.*thetadot.*cos(2*phi + theta)/4 + I_1.*I_2.*phidot.*thetadot - I_1.*I_2.*psidot.*thetadot.*cos(theta) + I_1.*I_3.*phidot.*psidot.*cos(2*phi - theta)/4 - I_1.*I_3.*phidot.*psidot.*cos(2*phi + theta)/4 + I_1.*I_3.*phidot.*thetadot.*cos(2*phi)/2 + I_1.*I_3.*phidot.*thetadot/2 + I_1.*I_3.*psidot.^2.*cos(2*phi - 2*theta)/8 - I_1.*I_3.*psidot.^2.*cos(2*phi + 2*theta)/8 + I_1.*I_3.*psidot.*thetadot.*cos(theta)/2 + I_1.*I_3.*psidot.*thetadot.*cos(2*phi - theta)/4 + I_1.*I_3.*psidot.*thetadot.*cos(2*phi + theta)/4 + I_2.^2.*phidot.*psidot.*cos(2*phi - theta)/4 - I_2.^2.*phidot.*psidot.*cos(2*phi + theta)/4 + I_2.^2.*phidot.*thetadot.*cos(2*phi)/2 - I_2.^2.*phidot.*thetadot/2 + I_2.^2.*psidot.^2.*cos(2*phi - 2*theta)/8 - I_2.^2.*psidot.^2.*cos(2*phi + 2*theta)/8 - I_2.^2.*psidot.*thetadot.*cos(theta)/2 + I_2.^2.*psidot.*thetadot.*cos(2*phi - theta)/4 + I_2.^2.*psidot.*thetadot.*cos(2*phi + theta)/4 - I_2.*I_3.*phidot.*psidot.*cos(2*phi - theta)/4 + I_2.*I_3.*phidot.*psidot.*cos(2*phi + theta)/4 - I_2.*I_3.*phidot.*thetadot.*cos(2*phi)/2 + I_2.*I_3.*phidot.*thetadot/2 - I_2.*I_3.*psidot.^2.*cos(2*phi - 2*theta)/8 + I_2.*I_3.*psidot.^2.*cos(2*phi + 2*theta)/8 + I_2.*I_3.*psidot.*thetadot.*cos(theta)/2 - I_2.*I_3.*psidot.*thetadot.*cos(2*phi - theta)/4 - I_2.*I_3.*psidot.*thetadot.*cos(2*phi + theta)/4)./(I_1.*I_2.*sin(theta));
        thetaddot = (I_1.^2.*phidot.*psidot.*sin(phi).^2.*sin(theta) + I_1.^2.*phidot.*thetadot.*sin(2*phi)/2 + I_1.^2.*psidot.^2.*sin(phi).^2.*sin(theta).*cos(theta) + I_1.^2.*psidot.*thetadot.*(sin(2*phi - theta) + sin(2*phi + theta))/4 - I_1.*I_2.*phidot.*psidot.*sin(phi).^2.*sin(theta) - I_1.*I_2.*phidot.*psidot.*sin(theta).*cos(phi).^2 - I_1.*I_3.*phidot.*psidot.*sin(phi).^2.*sin(theta) - I_1.*I_3.*phidot.*thetadot.*sin(2*phi)/2 - I_1.*I_3.*psidot.^2.*sin(phi).^2.*sin(theta).*cos(theta) - I_1.*I_3.*psidot.*thetadot.*(sin(2*phi - theta) + sin(2*phi + theta))/4 + I_2.^2.*phidot.*psidot.*sin(theta).*cos(phi).^2 - I_2.^2.*phidot.*thetadot.*sin(2*phi)/2 + I_2.^2.*psidot.^2.*sin(theta).*cos(phi).^2.*cos(theta) - I_2.^2.*psidot.*thetadot.*(sin(2*phi - theta) + sin(2*phi + theta))/4 - I_2.*I_3.*phidot.*psidot.*sin(theta).*cos(phi).^2 + I_2.*I_3.*phidot.*thetadot.*sin(2*phi)/2 - I_2.*I_3.*psidot.^2.*sin(theta).*cos(phi).^2.*cos(theta) + I_2.*I_3.*psidot.*thetadot.*(sin(2*phi - theta) + sin(2*phi + theta))/4)./(I_1.*I_2);
        phiddot = (I_1.*psidot.^2.*sin(phi).*sin(theta).^2.*cos(phi) - 2*I_1.*psidot.*thetadot.*sin(phi).^2.*sin(theta) + I_1.*psidot.*thetadot.*sin(theta) - I_1.*thetadot.^2.*sin(phi).*cos(phi) - I_2.*psidot.^2.*sin(phi).*sin(theta).^2.*cos(phi) + 2*I_2.*psidot.*thetadot.*sin(phi).^2.*sin(theta) - I_2.*psidot.*thetadot.*sin(theta) + I_2.*thetadot.^2.*sin(phi).*cos(phi) - I_3.*psiddot.*cos(theta) + I_3.*psidot.*thetadot.*sin(theta))./I_3;
        
        dz = [psidot;thetadot;phidot;psiddot;thetaddot;phiddot];
    end

end