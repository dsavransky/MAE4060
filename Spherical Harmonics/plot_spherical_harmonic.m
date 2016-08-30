function plot_spherical_harmonic(x,y,z,f,topo)

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