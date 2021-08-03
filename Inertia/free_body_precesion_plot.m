% Generate a plot of angular velocity components of a torque free rigid
% body using numerical integration and the analytical solution.

% Copyright (c) 2017 Dmitry Savransky (ds264@cornell.edu)

%set up constants and initial conditions
It = 1;
Is = 2;
I = [It,It,Is];
ws = 1;
w0 = 0.1;

%numerically integrate
[t,w] = torque_free_motion([0,10],[w0,0,ws],[It,It,Is]);

figure(1)
clf
plot(t,w(:,1:2)) %third components is constant, so omit it in plot

wn = (Is - It)/It*ws; %find natural frequency of system

%overplot analytical solution
hold on
t1 = linspace(t(1),t(end),50);
plot(t1,w0*cos(wn*t1),'o',t1,w0*sin(wn*t1),'o')
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
legend({'$\omega_1$','$\omega_2$'},'Interpreter','LaTeX','FontSize',16)