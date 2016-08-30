%grid
load('topo.mat','topo','topomap1');
phi = linspace(0,pi,50);
theta = linspace(0,2*pi,50);

%% base
%l = 0,m=0
[x0,y0,z0] = calc_spherical_harmonic(0,0,phi,theta);
plot_spherical_harmonic(x0,y0,z0,1)
plot_spherical_harmonic(-x0,-y0,-z0,2,topo)

pause 

%% zonal
%l = 2,m=0
[x1,y1,z1] = calc_spherical_harmonic(2,0,phi,theta);
plot_spherical_harmonic(x1,y1,z1,1)
plot_spherical_harmonic(x1/10-x0,y1/10-y0,z1/10-z0,2,topo)

pause

%% zonal
%l = 4,m=0
[x1,y1,z1] = calc_spherical_harmonic(4,0,phi,theta);
plot_spherical_harmonic(x1,y1,z1,1)
plot_spherical_harmonic(x1/10-x0,y1/10-y0,z1/10-z0,2,topo)

pause

%% sectoral
%l = 2,m=2
[x1,y1,z1] = calc_spherical_harmonic(2,2,phi,theta);
plot_spherical_harmonic(x1,y1,z1,1)
plot_spherical_harmonic(x1/10-x0,y1/10-y0,z1/10-z0,2,topo)

pause
%% sectoral
%l = 3,m=3
[x1,y1,z1] = calc_spherical_harmonic(3,3,phi,theta);
plot_spherical_harmonic(x1,y1,z1,1)
plot_spherical_harmonic(x1/10-x0,y1/10-y0,z1/10-z0,2,topo)

pause
%% tesseral

%l = 6,m=4
[x1,y1,z1] = calc_spherical_harmonic(6,4,phi,theta);
plot_spherical_harmonic(x1,y1,z1,1)
plot_spherical_harmonic(x1/10-x0,y1/10-y0,z1/10-z0,2,topo)

return 
