function [t,z] = torque_free_motion2(t0,z0,I)
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

    %state  - [w1,w2,w3, ^B{C}^A(:)]
    [t,z] = ode113(@eulers_eqs,t0,z0,odeset('RelTol',1e-16,'AbsTol',1e-16));

    function dz = eulers_eqs(~,z)
        dw = [z(2)*z(3)*(I(2)-I(3))/I(1);...
              z(3)*z(1)*(I(3)-I(1))/I(2);...
              z(1)*z(2)*(I(1)-I(2))/I(3)];
        bCa = reshape(z(4:end),3,3);
        w = skew(z(1:3));
        dbCa = -w*bCa;
        dz = [dw;dbCa(:)];
    end

end