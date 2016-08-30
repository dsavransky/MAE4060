function plot_spherical_harmonic(x,y,z,f,topo)
% plot_spherical_harmonic creates a surface plot of a spherical harmonic
%
% plot_spherical_harmonic(x,y,z) generates the surface described by
% Cartesian coordinates x,y,x
% plot_spherical_harmonic(...,f) plots in figure with handle f
% plot_spherical_harmonic(...,topo) sets the surface Cdata to topo

% Copyright (c) 2014 Dmitry Savransky (ds264@cornell.edu)

if ~exist('f','var'), f = 1; end

props.EdgeColor = 'none';
props.FaceLighting = 'phong';
if exist('topo','var')
    props.AmbientStrength = 0.1;
    props.DiffuseStrength = 1;
    props.SpecularColorReflectance = .5;
    props.SpecularExponent = 20;
    props.SpecularStrength = 1;
    props.FaceColor= 'texture';
    props.Cdata = topo;
end

figure(f)
clf
surface(x,y,z,props)
axis tight equal off
view(3)
camzoom(1.5)