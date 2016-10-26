function [t,w] = torque_free_motion(t0,w0,I)

    %state  - [w1,w2,w3]
    [t,w] = ode45(@eulers_eqs,t0,w0,odeset('RelTol',1e-12,'AbsTol',1e-12));

    function dw = eulers_eqs(t,w)
        dw = [w(2)*w(3)*(I(2)-I(3))/I(1);...
              w(3)*w(1)*(I(3)-I(1))/I(2);...
              w(1)*w(2)*(I(1)-I(2))/I(3)];
        
    end

end