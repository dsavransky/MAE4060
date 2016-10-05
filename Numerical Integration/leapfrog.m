function [tout,yout,dyout] = leapfrog(odefun,tspan,y0,dy0,h)
% leapfrog solves ODEs using the leapfrog method.
%     [TOUT,YOUT,DYOUT] = leapfrog(ODEFUN,TSPAN,Y0,DY0,H) with TSPAN = [T0
%     TFINAL] integrates the system of differential equations y'' = f(t,y)
%     from time T0 to TFINAL with initial conditions Y0 and time step H.
%     ODEFUN is a function handle. For a scalar T and a vector Y,
%     ODEFUN(T,Y) must return a column vector corresponding to f(t,y). Each
%     row in the solution array YOUT corresponds to a time returned in the
%     column vector TOUT.  To obtain solutions at specific times
%     T0,T1,...,TFINAL (all increasing or all decreasing), use TSPAN = [T0
%     T1 ... TFINAL]. Modeled on MATLAB's built-in ODE functions.

% input error checks

%require 4 inputs
narginchk(5,5);
assert(isa(odefun, 'function_handle'),...
    'leapfrog:odefun_not_fun', ...
    'Derivative f(t,y) must be given as function handle.');
assert(~isempty(tspan) && numel(tspan) >= 2,...
    'leapfrog:tspan_empty', ...
    'Time interval [tspan] must contain at least 2 values.');
assert(all(diff(tspan) ~= 0),...
    'leapfrog:tspan_err',...
    'Values in [tspan] must all span a non-zero interval.');
assert(all(diff(tspan) > 0) || all(diff(tspan) < 0), ...
    'leapfrog:tspan_monotonic',...
    'The entries in tspan must monotonically increase or decrease.');
assert(h <= abs(tspan(end) - tspan(1)), ...
    'leapfrog:time_step_large',...
    'The time step cannot be larger than the integration interval.');
assert(numel(y0) == numel(dy0), ...
    'leapfrog:time_step_large',...
    'The time step cannot be larger than the integration interval.');

%determine direction of time:
if tspan(end) > tspan(1), direction = 1; else direction = -1; end

%make sure time step is in the right direction:
h = abs(h)*direction;

%allocate outputs
if numel(tspan) > 2
    tout = tspan(:);
else
    tout = (tspan(1):h:tspan(end)).';
    if tout(end) ~= tspan(end), tout = [tout;tspan(end)]; end
end
yout = zeros(numel(tout),numel(y0));
dyout = zeros(numel(tout),numel(dy0));
yout(1,:) = y0(:).';
dyout(1,:) = dy0(:).';

%initialize
t = tout(1);
y = y0(:);
dy = dy0(:);
counter = 2; %index to current place in output

%attempt a first step to make sure that the given function gives you the
%expected result
try
    tmp = feval(odefun, t, y);
catch ME
    error('leapfrog:bad_odefun','Supplied ODEFUN does not work.')
end
if numel(tmp) ~= numel(y0)
    error('leapfrog:bad_odefun',...
        'ODEFUN should return a %3.0f-element column vector.', numel(y0))
end

%now integrate
while counter <= numel(tout)
   %figure out if you need to modify the current time step to hit the next
   %requested output:
   if abs(tout(counter) - t) < abs(h)
       dt = tout(counter) - t;
   else
       dt = h;
   end
   
   %update state
   fi = odefun(t,y);
   y = y + dt*dy + dt^2*fi/2;
   dy = dy + (fi + odefun(t,y))/2*dt;
   t = t + dt;
   
   %write output as needed
   if tout(counter) == t
       yout(counter,:) = y.';
       dyout(counter,:) = dy.';
       counter = counter + 1;
   end
end

end