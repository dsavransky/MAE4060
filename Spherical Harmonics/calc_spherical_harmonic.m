function [x,y,z] = calc_spherical_harmonic(l,m,phi,theta)
% calc_spherical_harmonic evaluates the spherical harmonic given by indices
% l,m over the ranges of spherical angles phi,theta (where theta is azimuth
% and phi is zenith angle).
%
% [x,y,z] =  calc_spherical_harmonic(l,m,phi,theta) returns matrices of
% Cartesian coordinates x,y,z equal in size to the outputs of
% meshgrid(theta,phi)

% Copyright (c) 2014 Dmitry Savransky (ds264@cornell.edu)

%setup grid
[Theta,Phi] = meshgrid(theta,phi);

% Calculate the harmonic
Pmls = legendre(l,cos(phi));
Pml = Pmls(m+1,:).';

U = Pml*(cos(m*theta)+sin(m*theta));
r = U/max(abs(U(:)));
x = r.*sin(Phi).*cos(Theta);
y = r.*sin(Phi).*sin(Theta);
z = r.*cos(Phi);