function anim_hill_orbits(reforb,hillorb_rot,M,dt)
%anim_hill_orbits animates relative motion orbits in the inertial and
%rotation (Euler-Hill) referece frames
%
%Inputs:
%   reforb      (2xn float) Reference orbit x and y positions
%   hillorb_rot (2xn float) Rotating frame orbit x and y positions
%   M           (1xn float) Ref. orbit mean anomalies
%   dt          (scalar float) Time step

% Copyright Dmitry Savransky (ds264@cornell.edu) 2016

hillorb = hillorb_rot;
for j = 1:length(M)
    hillorb(:,j) = [cos(M(j)), -sin(M(j));...
        sin(M(j)) ,cos(M(j))]*hillorb_rot(:,j);
end
hillorb = hillorb + reforb;

f1 = figure(1);
f1.Position = [0,f1.Position(2:4)];
clf
hold on
axis equal
plot(0,0,'k.','MarkerSize',50);
msr = plot(reforb(1,1),reforb(2,1),'b.','MarkerSize',20);
ms = plot(hillorb(1,1),hillorb(2,1),'r.','MarkerSize',40);
ps = plot(reforb(1,1),reforb(2,1),'b--',hillorb(1,1),hillorb(2,1),'r--');
hold off
set(gca,'Xlim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],...
    'YLim',[-max(abs(hillorb(:))),max(abs(hillorb(:)))],...
    'YTick',[],'Xtick',[])


f2 = figure(2);
f2.Position = [f1.Position(1)+f1.Position(3)+1,...
    f1.Position(2),f2.Position(3:4)];
clf
hold on
plot(0,0,'b.','MarkerSize',50);
ms2 = plot(hillorb_rot(1,1),hillorb_rot(2,1),'r.','MarkerSize',40);
ps2 = plot(hillorb_rot(1,1),hillorb_rot(2,1),'r--');
hold off
xlm = [min(hillorb_rot(1,:)),max(hillorb_rot(1,:))];
ylm = [min(hillorb_rot(2,:)),max(hillorb_rot(2,:))];
if xlm(1) == xlm(2)
    xlm = [xlm(1)-xlm(1)/2,xlm(1)+xlm(1)/2];
end
if ylm(1) == ylm(2)
    ylm = [ylm(1)-ylm(1)/2,ylm(1)+ylm(1)/2];
end
set(gca,'Xlim',xlm,...
    'YLim',ylm,...
    'YTick',[],'Xtick',[])


for j = 2:length(M)
    set(ms,'XData',hillorb(1,j),'YData',hillorb(2,j));
    set(msr,'XData',reforb(1,j),'YData',reforb(2,j));
    set(ps(1),'XData',reforb(1,1:j),'YData',reforb(2,1:j));
    set(ps(2),'XData',hillorb(1,1:j),'YData',hillorb(2,1:j));
    set(ms2,'XData',hillorb_rot(1,j),'YData',hillorb_rot(2,j));
    set(ps2,'XData',hillorb_rot(1,1:j),'YData',hillorb_rot(2,1:j));
    pause(dt);
end