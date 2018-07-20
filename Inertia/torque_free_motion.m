function [t,w] = torque_free_motion(t0,w0,I)
%TORQUE_FREE_MOTION integrates Euler's equations in the absence of external
%moments
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


    %state  - [w1,w2,w3]
    [t,w] = ode45(@eulers_eqs,t0,w0,odeset('RelTol',1e-12,'AbsTol',1e-12));

    function dw = eulers_eqs(t,w)
        dw = [w(2)*w(3)*(I(2)-I(3))/I(1);...
              w(3)*w(1)*(I(3)-I(1))/I(2);...
              w(1)*w(2)*(I(1)-I(2))/I(3)];
        
    end

end