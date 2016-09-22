%% case 1 - s = 0
n = 1;
r1 = 1;
x0 = 0.01;
y0 = 0;
yd0 = -3*n*x0/2;
t = linspace(0,4*pi,1000);
x = x0*ones(size(t));
y = -3*n*t*x0/2 + y0;
M = mod(n*t,2*pi);
reforb = [r1*cos(M);r1*sin(M)];
hillorb_rot = [x;y];
dt = 15/length(M);

anim_hill_orbits(reforb,hillorb_rot,M,dt)

%% case 2 - s = \pm ni
x0 = 0.05;
y0 = 0.05;
yd0 = -2*n*x0;

x = x0*cos(n*t);
y = -2*x0*sin(n*t) + y0;
hillorb_rot = [x;y];

anim_hill_orbits(reforb,hillorb_rot,M,dt)

%% case 3 - drift - x0=y0=xd0=0, yd0 >0
n = 3;
r1 = 1;
x0 = 0;
y0 = 0;
yd0 = 0.01;
t = linspace(0,4*pi,1000);
x = 2*(1-cos(n*t))*yd0/n;
y = yd0/n*(4*sin(n*t) - 3*n*t);

M = mod(n*t,2*pi);
reforb = [r1*cos(M);r1*sin(M)];
hillorb_rot = [x;y];
dt = 15/length(M);

anim_hill_orbits(reforb,hillorb_rot,M,dt)