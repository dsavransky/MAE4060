function [path1,path2] = vacFocus2transfer(r1v,r2v,Fvac,col,hyp)

if ~exist('col','var') || isempty(col)
    col = 'k';
end

if ~exist('hyp','var')
    hyp = false;
end

rotMats = DCMs();
r1 = norm(r1v); %magniutde of r1
Fshat = Fvac/norm(Fvac); %unit direction of Focus to vacant focus (-eccen vec)
rf = dot(r1v(1:2),Fshat); %projection of r1 onto Fshat
d = norm(Fvac); %distance between foci = 2*abs(a)*e

if hyp
    a = r1/2 - sqrt(d^2 - 2*d*rf + r1^2)/2;
    e = (r1 + sqrt(d^2 - 2*d*rf + r1^2))/(d - 2*rf);
    w = atan2(Fvac(2),Fvac(1));
else
    a = r1/2 + sqrt(d^2 - 2*d*rf + r1^2)/2;
    e = (-r1 + sqrt(d^2 - 2*d*rf + r1^2))/(d - 2*rf);
    w = atan2(Fvac(2),Fvac(1)) - pi;
end
ell = a*(1-e^2);

tmp = rotMats{3}(w)*r1v;
nu1 = atan2(tmp(2),tmp(1));
tmp = rotMats{3}(w)*r2v;
nu2 = atan2(tmp(2),tmp(1));
nu = linspace(nu1,nu2,100);
rmag = ell./(1 + e*cos(nu));
r = [rmag.*cos(nu);rmag.*sin(nu);zeros(size(rmag))];
path1 = rotMats{3}(-w)*r;

if ~hyp
    nu = linspace(nu2,nu1+2*pi,500);
    rmag = ell./(1 + e*cos(nu));
    r = [rmag.*cos(nu);rmag.*sin(nu);zeros(size(rmag))];
    path2 = rotMats{3}(-w)*r;
else
    path2 = [];
end

hold on
plot(Fvac(1),Fvac(2),'.','Color',col,'MarkerSize',20)
plot(path1(1,:),path1(2,:),'Color',col,'Linewidth',2)
if ~hyp,plot(path2(1,:),path2(2,:),'Color',col,'Linewidth',2);end
hold off
