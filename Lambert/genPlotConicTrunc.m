function rE = genPlotConicTrunc(Fs,a,e,f,r1v,r2v,col,branch)

w = atan2(Fs(2),Fs(1));
if a > 0
    %elliptical trajectory
    E = linspace(0,2*pi,1000);
    rE = a*[(cos(E) - e);sqrt(1 - e^2)*sin(E)];
    rE = [cos(w), -sin(w); sin(w), cos(w)]*rE+repmat(Fs.',1,length(E));
else
    %hyperbolic trajectory
    Eh = linspace(-pi,pi,1000);
    rE = -a*[(e - cosh(Eh));sqrt(e^2 - 1)*sinh(Eh)];
    rE = [cos(w), -sin(w); sin(w), cos(w)]*rE;
end

[~,ind1] = min(sqrt((rE(1,:) - r1v(1)).^2 + (rE(2,:) - r1v(2)).^2));
[~,ind2] = min(sqrt((rE(1,:) - r2v(1)).^2 + (rE(2,:) - r2v(2)).^2));

if branch == 1
    inds = [min([ind1,ind2]):max([ind1,ind2])];
else
    inds = [max([ind1,ind2]):length(rE),1:min([ind1,ind2])];
end
rE = rE(:,inds);

figure(f)
hold on
plot(rE(1,:),rE(2,:),'Color',col,'LineWidth',2)
hold off
axis equal