function genPlotConic(Fs,a,e,f,subsc,subscPos,col,stickax)

if ~exist('subsc','var') || isempty(subsc)
    subsc = '';
else
    subsc = ['_',subsc];
end
if ~exist('subscPos','var') || isempty(subscPos)
    subscPos = 'top';
end
if ~exist('col','var') || isempty(col)
    col = 'b';
end
if ~exist('stickax','var') || isempty(stickax)
    stickax = false;
end

w = atan2(Fs(2),Fs(1));
if a > 0
    %elliptical trajectory
    Fsb = Fs/norm(Fs);
    P = -Fsb*a*(1 + e) + Fs;
    A = Fsb*a*(1 - e) + Fs;
    
    E = linspace(0,2*pi,1000);
    rE = a*[(cos(E) - e);sqrt(1 - e^2)*sin(E)];
    rE = [cos(w), -sin(w); sin(w), cos(w)]*rE+repmat(Fs.',1,length(E));
else
    %hyperbolic trajectory
    Eh = linspace(-pi,pi,1000);
    rE = -a*[(e - cosh(Eh));sqrt(e^2 - 1)*sinh(Eh)];
    rE = [cos(w), -sin(w); sin(w), cos(w)]*rE;
end


figure(f)
hold on
if stickax
    set(gca,'XlimMode','manual','YlimMode','manual');
end
plot(Fs(1),Fs(2),[col,'.'])
plot(rE(1,:),rE(2,:),col)
if a > 0
    plot(A(1),A(2),'k.',P(1),P(2),'k.')
    plot([A(1),P(1)],[A(2),P(2)],'k--')
else
    plot([Fs(1),0],[Fs(2),0],'k--')
end
hold off
axis equal

shim = max(diff(axis))*0.02/3;
subxy = Fs;
switch lower(subscPos)
    case 'top'
        subxy = subxy + [0,shim*5];
    case 'bottom'
        subxy = subxy - [0,shim*5];
    otherwise
        subxy = subxy - [shim*5,shim*5];
end

text(subxy(1),subxy(2),['$$F^\star',subsc,'$$'],'HorizontalAlignment','center')