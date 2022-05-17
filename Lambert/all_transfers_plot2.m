% Draw figure with all departure regions and transfer orbit types for the
% Lambert Problem.

% Copyright (c) 2009 Dmitry Savransky (ds264@cornell.edu)
% Updated 2022-05-10

% define the vectors
r1 = 1;
r2 = 1.8;
dth = 120*pi/180;

%calculate the rest of the triangle
c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dth)); %\Vert \mf r_{1/2} \Vert
th1 = asin(r1*sin(dth)/c);
th2 = asin(r2*sin(dth)/c);
%r1 and r2 vectors
r1v = r1*[-cos(th2), sin(th2)];
r2v = r2*[ cos(th1), sin(th1)];

%define rotation matrix:
rotMat = @(ang) [cos(ang), -sin(ang); sin(ang), cos(ang)];

%define locus about P1 for plotting departure regions:
tmpth =  linspace(0,2*pi,500);
depcircrad = c/3.5;
depcirc = depcircrad*[cos(tmpth);sin(tmpth)]+repmat(r1v.',1,length(tmpth));

bottom = -1.2;
top = 2;
%% minimum energy orbit
am = (r1 + r2 + c)/4;
s = 2*am; %equivalent to (r1 + r2 + c)/2
Fsm = r1v + [s - r1, 0]; %vacant focus location of min E transfer
ellm = 2/c*(s - r1)*(s - r2); %semi-param of min E transfer
em = sqrt(1 - 2*ellm/s); %eccentricity of min E transfer
wm = atan2(Fsm(2),Fsm(1)) - pi; %arg of periapsis of min E transfer

% generate full orbit
nu = linspace(0,2*pi,500);
rmag = am*(1 - em^2)./(1 + em*cos(nu));
rm = rotMat(wm)*[rmag.*cos(nu);rmag.*sin(nu)];

%find intersections with departure circle
[rm1, rmcirc1, rmcirc2, rmcinds] = ...
    findLambertIntersections(r1v,r2v,depcirc,rm);
inds = sort([find(rm(1,:) == rm1(1,1) & rm(2,:) == rm1(2,1)),...
    find(rm(1,:) == rm1(1,end) & rm(2,:) == rm1(2,end))]);
inds = inds + [-1,1];
rm2 = rm(:,[inds(2):end,1:inds(1)]);

%% all vacant foci lie on this hyperbola:
ah = -abs((r1 - r2))/2;
eh = c/abs(r2 - r1);
Eh = linspace(-pi/2.5,pi/3,1000);
rE = -ah*[(eh - cosh(Eh));sqrt(eh^2 - 1)*sinh(Eh)];

%split vacant foci hyperbolae:
if r2 > r1
    lhyp = [rE(1,:)+r1v(1);rE(2,:)+r1v(2)];
    rhyp = [-rE(1,:)+r2v(1);rE(2,:)+r2v(2)];
else
    rhyp = [rE(1,:)+r1v(1);rE(2,:)+r1v(2)];
    lhyp = [-rE(1,:)+r2v(1);rE(2,:)+r2v(2)];
end

%E1/E3 and E2/E4 branches split at F^\star_min
[~,Fmind] = min(sqrt((rhyp(1,:) - Fsm(1)).^2 + (rhyp(2,:) - Fsm(2)).^2));

% H_1 branch starts at F_0 focus (a = 0)
a1 = (r1^2 - r2^2 + c^2)/(2*c);
h1 = sqrt(r1^2 - a1^2);
G = r1v + a1/c*(r2v - r1v);
Fp  = [h1*(r2v(2) - r1v(2))/c,-h1*(r2v(1) - r1v(1))/c];
F0 = G - Fp;
[~,F0ind] = min(sqrt((lhyp(1,:) - F0(1)).^2 + (lhyp(2,:) - F0(2)).^2));

% H_2 branch starts at Filled focus (0,0)
[~,Find] = min(sqrt(sum(lhyp.^2))); %index of filled Focus (0,0)

%% parabolic trajectories (escape velocity curves)
ellp1 = 4*(s-r1)*(s-r2)/c^2*(sqrt(s/2) + sqrt((s-c)/2))^2;
ellp2 = 4*(s-r1)*(s-r2)/c^2*(sqrt(s/2) - sqrt((s-c)/2))^2;
wp = acos((r2-r1)/c) - pi; %rotation of parabolae

nu = linspace(-126*pi/180,105*pi/180,1000);
B2 = tan(nu/2).^2;
rmag1 = ellp1/2*(1+B2);
rp1 = rotMat(-wp)*[rmag1.*cos(nu);rmag1.*sin(nu)];
[~,ind1] = min(sqrt((rp1(1,:) - r1v(1)).^2 + (rp1(2,:) - r1v(2)).^2));
[~,ind2] = min(sqrt((rp1(1,:) - r2v(1)).^2 + (rp1(2,:) - r2v(2)).^2));
inds = sort([ind1,ind2]);
rp1r = rp1(:,1:inds(1)-1);
rp1l = rp1(:,inds(2)+1:end);
rp1mid = rp1(:,inds(1):inds(2));

nu = linspace(-135*pi/180,145*pi/180,1000);
B2 = tan(nu/2).^2;
rmag2 = ellp2/2*(1+B2);
rp2 = rotMat(wp)*[rmag2.*cos(nu);rmag2.*sin(nu)];

[~,ind1] = min(sqrt((rp2(1,:) - r1v(1)).^2 + (rp2(2,:) - r1v(2)).^2));
[~,ind2] = min(sqrt((rp2(1,:) - r2v(1)).^2 + (rp2(2,:) - r2v(2)).^2));
inds = sort([ind1,ind2]);
rp2l = rp2(:,1:inds(1)-1);
rp2r = rp2(:,inds(2)+1:end);
rp2mid = rp2(:,inds(1):inds(2));

%% departure regions (only shaded within local circl aroudn P1)

% % H1 departure region is bounded by r_1/2 and rp1
% [~,ind1] = min(rp1cinds); %looking for intersection closest to circle start
% if ind1 == 1
%     H1rp = rp1circ1;
%     E3rp = rp1circ2;
% else
%     H1rp = rp1circ2;
%     E3rp = rp1circ1;
% end
% H1 = [r1v.', r1v.' + [depcircrad; 0], depcirc(:,1:rp1cinds(ind1)),H1rp];
% 
% 
% [~,ind2] = min(rmcinds); %again looking for closer one
% if ind2 == 1
%     E1rm = rmcirc1;
%     E3rm = rmcirc2;
% else
%     E1rm = rmcirc2;
%     E3rm = rmcirc1;
% end
% E1 = [r1v.', fliplr(H1rp), depcirc(:,rp1cinds(ind1):rmcinds(ind2)),E1rm];
% 
% % E2 is bounded by vmin and rp2
% [~,ind3] = min(rp2cinds); %again looking for closer one
% if ind3 == 1
%     E2rp = rp2circ1;
%     E4rp = rp2circ2;
% else
%     E2rp = rp2circ2;
%     E4rp = rp2circ1;
% end
% E2 = [r1v.', fliplr(E1rm), depcirc(:,rmcinds(ind2):rp2cinds(ind3)), E2rp];
% 
% % E3 is bounded by rp1 and vmin
% E3 = [r1v.', E3rp, depcirc(:,rp1cinds(3-ind1):rmcinds(3-ind2)), fliplr(E3rm)];
% 
% % E4 is bounded by vmin and rp2
% E4 = [r1v.', E3rm, depcirc(:,rmcinds(3-ind2):rp2cinds(3-ind3)), fliplr(E4rp)];
% 
% % finally, H2 is bounded by rp2 and r1
% [x0,y0,~,cind] = intersections([0,r1v(1)],[0,r1v(2)],depcirc(1,:),depcirc(2,:));
% H2 = [r1v.', E4rp, depcirc(:,rp2cinds(3-ind3):round(cind)), r1v.'];

%% alt departure regions  (shaded throuhgout whole trajectory region)

% H1 departure region is bounded by r_1/2 and rp1
H1 = [r1v.', r2v.' , rp1mid];
% E1 is bounded by rp1 and vmin
E1 = [fliplr(rp1mid),rm1];
% E2 is bounded by vmin and rp2
E2 = [fliplr(rp2l),fliplr(rp2r),rm1];
% E3 is bounded by rp1 and vmin
E3 = [rp1l, rp1r, fliplr(rm2)];
% E4 is bounded by vmin and rp2
E4 = [rm2,fliplr(rp2mid)];
% finally, H2 is bounded by rp2 and r1 and r2
H2 = [r1v.',[0;0],r2v.',rp2mid];
%% plot 
H1color = [255,140,0]/255;
H2color = [85,107,47]/255;
rpcolor = [255, 0, 189]/255;

f = figure(2);
clf
f.Position = [ -428        1011        1546        1286];
genInitTriangleFig(r1v,r2v,2) %system steup
hold on

%foci hyporbolae:
p1 = plot(lhyp(1,1:Find),lhyp(2,1:Find),'b',...
    lhyp(1,Find:F0ind),lhyp(2,Find:F0ind),'k--',...
    lhyp(1,F0ind:end),lhyp(2,F0ind:end),'r','LineWidth',2);
set(p1(1),'Color',[85,107,47]/255);
set(p1(3),'Color',[255,140,0]/255);
p2 = plot(rhyp(1,1:Fmind),rhyp(2,1:Fmind),'b',...
    rhyp(1,Fmind:end),rhyp(2,Fmind:end),'r','LineWidth',2);
plot(Fsm(1),Fsm(2),'b.','MarkerSize',30)
plot(F0(1),F0(2),'.','Color',H1color,'MarkerSize',30)

% min energy transfer
p3 = plot(rm(1,:),rm(2,:),'b','LineWidth',2);

%parabolae:
p4 = plot(rp1(1,:),rp1(2,:),rp2(1,:),rp2(2,:),'LineWidth',2,'Color',rpcolor);

% departure regions
deps = fill(H1(1,:),H1(2,:),H1color, E1(1,:),E1(2,:),'b',...
    E2(1,:),E2(2,:),'r', E3(1,:),E3(2,:),'b', E4(1,:),E4(2,:),'r',...
    H2(1,:),H2(2,:),H2color,'Edgecolor','none','FaceAlpha',0.3);
uistack([p1;p2;p3;p4;deps],'bottom') 
axis([-2,4.25,-1.2,2])
set(gca,'position',[0 0 1 1])
hold off

%% annotate
shim = max(diff(axis))*0.02/3;
text(Fsm(1)+shim,Fsm(2)+shim,'$$F^\star_m$$','HorizontalAlignment','left',...
    'VerticalAlignment','bottom')
text(F0(1)+shim,F0(2)+shim,'$$F^\star_0$$','HorizontalAlignment','left',...
    'VerticalAlignment','top')
text(rhyp(1,50)+shim,rhyp(2,50),'$$E_1,E_3$$',...
    'HorizontalAlignment','left','Color','b')
text(rhyp(1,end-50)+shim,rhyp(2,end-50),'$$E_2, E_4$$',...
    'HorizontalAlignment','left','Color','r')
text(lhyp(1,40)+shim,lhyp(2,40),'$$H_2$$',...
    'HorizontalAlignment','left','Color',H2color)
text(lhyp(1,end-50)+shim, lhyp(2,end-50),'$$H_1$$',...
    'HorizontalAlignment','left','Color',H1color)
text((r1v(1) + r2v(1))/2,(max(rp1(2,:))+r1v(2))/2,'$$H_1$$',...
    'HorizontalAlignment','center','VerticalAlignment','middle','Color',H1color)
text((r1v(1) + r2v(1))/2,(max(rm1(2,:))+max(rp1(2,:)))/2,'$$E_1$$',...
    'HorizontalAlignment','left','VerticalAlignment','middle','Color','b')
text(rp2(1,end-40),(rp2(2,end-40)+top)/2,'$$E_2$$',...
    'HorizontalAlignment','left','VerticalAlignment','middle','Color','r')
text(rp1(1,125),(rp1(2,125)+bottom)/2,'$$E_3$$',...
    'HorizontalAlignment','left','VerticalAlignment','middle','Color','b')
text((r1v(1) + r2v(1))/2,(min(rm2(2,:))+min(rp2(2,:)))/2,'$$E_4$$',...
    'HorizontalAlignment','left','VerticalAlignment','middle','Color','r')
text((r1v(1) + r2v(1))/2,min(rp2(2,:))/2,'$$H_2$$',...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Color',H2color)
tmp = (r1v+r2v)/2;
text(tmp(1),tmp(2),'$$v = \infty$$','HorizontalAlignment',...
    'center','VerticalAlignment','top','Color',H1color)
[~,ind] = min(rp2(2,:));
text(rp2(1,ind)+shim*5,rp2(2,ind),'$$v = v_{esc}$$','HorizontalAlignment',...
    'center','VerticalAlignment','top','Color',rpcolor)
[~,ind] = max(rp1(2,:));
text(rp1(1,ind)-shim,rp1(2,ind),'$$v = v_{esc}$$','HorizontalAlignment',...
    'left','VerticalAlignment','bottom','Color',rpcolor)
[~,ind] = max(rm(2,:));
text(rm(1,ind),rm(2,ind),'$$v = v_{min}$$','HorizontalAlignment',...
    'center','VerticalAlignment','bottom','Color','b')
[~,ind] = min(rm(2,:));
text(rm(1,ind),rm(2,ind),'$$v = v_{min}$$','HorizontalAlignment',...
    'center','VerticalAlignment','top','Color','b')
 