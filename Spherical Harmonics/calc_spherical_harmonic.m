function [x,y,z] = calc_spherical_harmonic(l,m,phi,theta)

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