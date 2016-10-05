%% Earth's orbit (two-body)

%define the orbit:
a = 1; %AU
e = 0.0167; 
Tp = 1; %year

%calculate \mu and n from given info:
mu = a^3*(Tp/2/pi)^(-2); %AU^3/year^2
n = sqrt(mu/a^3); %yr^(-1)
h = sqrt(mu*a*(1-e^2));

%define 1 year of true anomaly:
np = 10;
M = linspace(0,np*2*pi,np*100);
%get the eccentric anomaly:
[E,nu] = invKepler(mod(M,2*pi),e);

%the (heliocentric) position and velocity in the perifocal frame
%(only need two coordinates as we're in the plane of the orbit
r = a*(1-e^2)./(1+ e*cos(nu));

r_ps = [r.*cos(nu),r.*sin(nu)];
v_ps = [-mu/h*sin(nu), mu/h*(e+cos(nu))];

%get the intial position and velocity
r_ps_0 = r_ps(1,:);
v_ps_0 = v_ps(1,:);
%figure out the time (M = n(t - t_0) - choose t_0 = 0)
t = M/n;
%% RK45
%state space x = [r_1,r_2,v_1,v_2].';
%dr^2/dt^2 = -mu r/|r|^3
newton_gravity = @(t,x) [x(3);x(4);-mu*x(1:2)/norm(x(1:2))^3];
x_0 = [r_ps_0.';v_ps_0.'];

%integrate
[tout,xout] = ode45(newton_gravity,t,x_0);
r_ps_int = xout(:,1:2);

%compare
figure(1)
clf
subplot(2,1,1)
plot(r_ps(:,1),r_ps(:,2),r_ps_int(:,1),r_ps_int(:,2))
set(gca,'FontName','Times','FontSize',14)
xlabel('$\hat{\mathbf{e}}$ (AU)','Interpreter','Latex')
ylabel('$\hat{\mathbf{q}}$ (AU)','Interpreter','Latex')
legend({'Analytical','Numerical'},'Location','NorthWest')
subplot(2,1,2)
plot(t,r_ps(:,1)-r_ps_int(:,1),t,r_ps(:,2)-r_ps_int(:,2))
set(gca,'FontName','Times','FontSize',14)
xlabel('Time (years)')
ylabel('Error (AU)')
legend({'$\hat{\mathbf{e}}$','$\hat{\mathbf{q}}$ (AU)'},...
    'Interpreter','Latex','Location','NorthWest')

%% leapfrog
%state space x = [r_1,r_2].';
%dr^2/dt^2 = -mu r/|r|^3
newton_gravity2 = @(t,x) -mu*x(1:2)/norm(x(1:2))^3;

%integrate
dt = mean(diff(t));
[tout2,xout2,dxout2] = leapfrog(newton_gravity2,t,r_ps_0.',v_ps_0.',dt);
r_ps_int = xout2(:,1:2);

%compare
figure(2)
clf
subplot(2,1,1)
plot(r_ps(:,1),r_ps(:,2),r_ps_int(:,1),r_ps_int(:,2))
set(gca,'FontName','Times','FontSize',14)
xlabel('$\hat{\mathbf{e}}$ (AU)','Interpreter','Latex')
ylabel('$\hat{\mathbf{q}}$ (AU)','Interpreter','Latex')
legend({'Analytical','Numerical'},'Location','NorthWest')
subplot(2,1,2)
plot(t,r_ps(:,1)-r_ps_int(:,1),t,r_ps(:,2)-r_ps_int(:,2))
set(gca,'FontName','Times','FontSize',14)
xlabel('Time (years)')
ylabel('Error (AU)')
legend({'$\hat{\mathbf{e}}$','$\hat{\mathbf{q}}$ (AU)'},...
    'Interpreter','Latex','Location','NorthWest')


%% Total energy
%energy = 1/2 v^2 + sum_j>i mu_j/|r_ij|


E_rk = sum(xout(:,3:4).^2,2)/2 - mu./sqrt(sum(xout(:,1:2).^2,2));
E_leapfrog = sum(dxout2.^2,2)/2 - mu./sqrt(sum(xout2.^2,2));

figure(3)
clf
plot(tout,E_rk,tout,E_leapfrog)
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (days)','Interpreter','Latex')
ylabel('Specific Energy (kg AU$^2$/day$^2$','Interpreter','Latex')
legend({'ODE45','Leapfrog'},'Interpreter','Latex','Location','SouthWest')
