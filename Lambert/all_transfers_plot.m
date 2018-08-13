%Draw figure with all departure regions and transfer orbit types

% define the vectors
r1 = 1;
r2 = 1.8;
dth = 120*pi/180;

%calculate the rest of the triangle
c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dth));
th1 = asin(r1*sin(dth)/c);
th2 = asin(r2*sin(dth)/c);
%r1 and r2 vectors
r1v = r1*[-cos(th2), sin(th2)];
r2v = r2*[cos(th1),sin(th1)];

%minimum energy orbit
am = (r1 + r2 + c)/4;
s = 2*am;
Fs = r1v + [s - r1, 0];
pm = 2/c*(s - r1)*(s - r2);
em = sqrt(1 - 2*pm/s);

%define rotation matrix:
rotMat = @(ang) [cos(ang), -sin(ang); sin(ang), cos(ang)];

%all vacant foci lie on this hyperbola:
ah = -abs((r1 - r2))/2;
eh = c/abs(r2 - r1);
bh = ah*sqrt(eh^2 - 1);
Eh = linspace(-pi/2.5,pi/3,1000);
rE = -ah*[(eh - cosh(Eh));sqrt(eh^2 - 1)*sinh(Eh)];

%E1, E3 branch starts at intersection of vacant foci hyperbola and segment
%OP_2 
x0 = r2v(1);
y0 = r2v(2);
if r2 > r1
    xi = (x0*(bh^2*x0*(ah*eh + x0) - ah^2*y0^2) + ...
        sqrt(ah^2*bh^4*x0^4 + ah^4*bh^2*(-1 + eh^2)*x0^2*y0^2))/...
        ((bh*x0 - ah*y0)*(bh*x0 + ah*y0));
    yi = y0/x0*xi;
else
    x1 = r1v(1);
    y1 = r1v(2);
    xi = (x1*(bh^2*(ah*eh + x0)*x1 - ah^2*y0*y1) - ...
        sqrt(ah^2*bh^2*x1^2*(bh^2*x1^2 + (x1*y0 - ...
        (ah*(-1 + eh)+x0)*y1)*(x1*y0 -(ah + ah*eh + x0)*y1))))/...
        ((bh*x1 - ah*y1)*(bh*x1 + ah*y1));
    yi = y1/x1*xi;
end
aE1 = (norm(r2v - [xi,yi]) + r2)/2;
eE1 = sqrt(xi^2 + yi^2)/2/aE1;
wE1 = atan2(yi,xi);
uE1 = rotMat(wE1)*[sin(wE1+th2); eE1+cos(wE1+th2)];
uE1 = uE1/norm(uE1);

% H_1 branch starts at F_0 focus (a = 0)
a1 = (r1^2 - r2^2 + c^2)/(2*c);
h1 = sqrt(r1^2 - a1^2);
G = r1v + a1/c*(r2v - r1v);
Fp  = [h1*(r2v(2) - r1v(2))/c,-h1*(r2v(1) - r1v(1))/c];
F0 = G - Fp;

%parabolic trajectories (escape velocity curves)
if r2 > r1
    th = pi - acos((r2-r1)/c);
else
    th = acos((r1-r2)/c);
end
h1 = r1*(1 - cos(th - th2));
h2 = r1*(1 - cos(th + th2));
xs = linspace(-2*max([r2,r1]),2*max([r2,r1]),1000);
rp1 = [sin(th), -cos(th); cos(th), sin(th)]*[xs;xs.^2/2/h1 - h1/2];
rp2 = [-sin(th), -cos(th); cos(th), -sin(th)]*[xs;xs.^2/2/h2 - h2/2];
%truncate parabolas:
[~,ind1] = min(sqrt((rp1(1,:) - r1v(1)).^2 + (rp1(2,:) - r1v(2)).^2));
[~,ind2] = min(sqrt((rp1(1,:) - r2v(1)).^2 + (rp1(2,:) - r2v(2)).^2));
uH2 = mean(diff(rp1(:,ind1-1:ind1+1),[],2),2);
uH2 = uH2/norm(uH2);
rp1tail = rp1(:,1:ind1);
rp1 = rp1(:,ind1:ind2);
[~,ind1] = min(sqrt((rp2(1,:) - r1v(1)).^2 + (rp2(2,:) - r1v(2)).^2));
[~,ind2] = min(sqrt((rp2(1,:) - r2v(1)).^2 + (rp2(2,:) - r2v(2)).^2));
uH1 = mean(diff(rp2(:,ind1+1:-1:ind1-1),[],2),2);
uH1 = uH1/norm(uH1);
rp2tail = rp2(:,ind1:end);
rp2 = rp2(:,ind2:ind1);

%split vacant foci hyperbolae:
if r2 > r1
    lhyp = [rE(1,:)+r1v(1);rE(2,:)+r1v(2)];
    rhyp = [-rE(1,:)+r2v(1);rE(2,:)+r2v(2)];
else
    rhyp = [rE(1,:)+r1v(1);rE(2,:)+r1v(2)];
    lhyp = [-rE(1,:)+r2v(1);rE(2,:)+r2v(2)];
end
[~,Oind] = min(sqrt(sum(lhyp.^2)));
[~,F0ind] = min(sqrt((lhyp(1,:) - F0(1)).^2 + (lhyp(2,:) - F0(2)).^2));
[~,Fmind] = min(sqrt((lhyp(1,:) - Fs(1)).^2 + (lhyp(2,:) - Fs(2)).^2));
[~,F1ind] = min(sqrt((lhyp(1,:) - xi).^2 + (lhyp(2,:) - yi).^2));

%plot 
f = figure(1);
clf
f.Position = [319   448  840 630];
hold on
%foci hyporbolae:
p1 = plot(lhyp(1,1:Oind),lhyp(2,1:Oind),'b',...
    lhyp(1,Oind:F0ind),lhyp(2,Oind:F0ind),'k--',...
    lhyp(1,F0ind:end),lhyp(2,F0ind:end),'r','LineWidth',2);
set(p1(1),'Color',[85,107,47]/255);
set(p1(3),'Color',[255,140,0]/255);
p2 = plot(rhyp(1,1:F1ind),rhyp(2,1:F1ind),'b',...
    rhyp(1,F1ind:Fmind),rhyp(2,F1ind:Fmind),'k--',...
    rhyp(1,Fmind:end),rhyp(2,Fmind:end),'r','LineWidth',2);
plot(Fs(1),Fs(2),'b.',F0(1),F0(2),'b.',xi,yi,'r.','MarkerSize',20)
%parabolae:
plot(rp1(1,:),rp1(2,:),'Color',[85,107,47]/255,'LineWidth',2)
plot(rp2(1,:),rp2(2,:),'Color',[255,140,0]/255,'LineWidth',2)

genInitTriangleFig(r1v,r2v,1)
rE1 = genPlotConicTrunc(Fs,am,em,1,r1v,r2v,'b',2);
rE2 = genPlotConicTrunc([xi,yi],aE1,eE1,1,r1v,r2v,'r',2);

hold on

%departure regions:
de = max(diff(axis))*0.1;
sliceInds = @(u1,u2) [r1v.',...
    de*[cos(linspace(atan2(u1(2),u1(1)),atan2(u2(2),u2(1)),100));...
        sin(linspace(atan2(u1(2),u1(1)),atan2(u2(2),u2(1)),100))] + ...
        repmat(r1v.',1,100),r1v'];

tmpth =  linspace(0,2*pi,100);
tmpcirc = c/10*[cos(tmpth);sin(tmpth)]+repmat(r1v.',1,100);

%H1 departure region
[x0,y0] = intersections(rp2(1,:),rp2(2,:),tmpcirc(1,:),tmpcirc(2,:));
[~,pind] = min(sum((rp2 - repmat([x0;y0],1,length(rp2))).^2,1));
[~,cind] = min(sum((tmpcirc - repmat([x0;y0],1,100)).^2,1));
H1 = [r1v.', r1v.' + [c/10;0],tmpcirc(:,1:cind),rp2(:,pind:end)];

%E1 departure region
[x0,y0] = intersections(rE1(1,:),rE1(2,:),tmpcirc(1,:),tmpcirc(2,:));
[~,e1ind] = min(sum((rE1 - repmat([x0;y0],1,length(rE1))).^2,1));
[~,cind2] = min(sum((tmpcirc - repmat([x0;y0],1,100)).^2,1));
E1 = [r1v.', rp2(:,pind:end), tmpcirc(:,cind:cind2),rE1(:,e1ind:end)];

%H2 departure region
[x0a,y0a] = intersections([0,r1v(1)],[0,r1v(2)],tmpcirc(1,:),tmpcirc(2,:));
[~,cind3] = min(sum((tmpcirc - repmat([x0a;y0a],1,100)).^2,1));
[x0,y0] = intersections(rp1(1,:),rp1(2,:),tmpcirc(1,:),tmpcirc(2,:));
[~,pind2] = min(sum((rp1 - repmat([x0;y0],1,length(rp1))).^2,1));
[~,cind4] = min(sum((tmpcirc - repmat([x0;y0],1,100)).^2,1));
H2 = [r1v.',[x0a;y0a],fliplr(tmpcirc(:,cind4:cind3)),...
    fliplr(rp1(:,1:pind2))];

%E4 departure region
[x0,y0] = intersections(rE2(1,:),rE2(2,:),tmpcirc(1,:),tmpcirc(2,:));
[~,e2ind] = min(sum((rE2 - repmat([x0;y0],1,length(rE2))).^2,1));
[~,cind5] = min(sum((tmpcirc - repmat([x0;y0],1,100)).^2,1));
E4 = [r1v.', rp1(:,1:pind2), fliplr(tmpcirc(:,cind5:cind4)),...
    fliplr(rE2(:,1:e2ind))];

%E3 departure region
[x0,y0] = intersections(rp2tail(1,:),rp2tail(2,:),tmpcirc(1,:),...
    tmpcirc(2,:));
[~,pind3] = min(sum((rp2tail - repmat([x0;y0],1,length(rp2tail))).^2,1));
[~,cind6] = min(sum((tmpcirc - repmat([x0;y0],1,100)).^2,1));
E3 = [r1v.',rE2(:,1:e2ind),fliplr(tmpcirc(:,cind6:cind5)),...
    fliplr(rp2tail(:,1:pind3))];

%E2 departure region
[x0,y0] = intersections(rp1tail(1,:),rp1tail(2,:),tmpcirc(1,:),...
    tmpcirc(2,:));
[~,pind4] = min(sum((rp1tail - repmat([x0;y0],1,length(rp1tail))).^2,1));
[~,cind7] = min(sum((tmpcirc - repmat([x0;y0],1,100)).^2,1));
E2 = [r1v.',fliplr(rE1(:,e1ind:end)),tmpcirc(:,cind2:cind7),...
    rp1tail(:,pind4:end)];

fill(H1(1,:),H1(2,:),[255,140,0]/255, E1(1,:),E1(2,:),'b',...
    H2(1,:),H2(2,:),[85,107,47]/255, E4(1,:),E4(2,:),'r',...
    E3(1,:),E3(2,:),'b',E2(1,:),E2(2,:),'r','Edgecolor','none');

% H1 = sliceInds([1,0],uH1);
% E1 = sliceInds(uH1,uE1);
% E2 = sliceInds(uE1,-uH2);
% H2 = sliceInds(-r1v/norm(r1v),uH2);
% E4 = sliceInds(uH2,-uE1);
% E3 = sliceInds(-uE1,-uH1);
% fill(H1(1,:),H1(2,:),[255,140,0]/255,E1(1,:),E1(2,:),'b',...
%     E2(1,:),E2(2,:),'r',H2(1,:),H2(2,:),[85,107,47]/255,...
%     E4(1,:),E4(2,:),'r',E3(1,:),E3(2,:),'b','Edgecolor','none')

hold off


%annotate
shim = max(diff(axis))*0.02/3;
text(Fs(1)+shim,Fs(2)+shim,'$$F^\star_m$$','HorizontalAlignment','left',...
    'VerticalAlignment','bottom')
text(F0(1)+shim,F0(2)+shim,'$$F^\star_0$$','HorizontalAlignment','left',...
    'VerticalAlignment','top')
text(rhyp(1,round(F1ind/2))+shim,rhyp(2,round(F1ind/2)),'$$E_1,E_3$$',...
    'HorizontalAlignment','left','Color','b')
text(rhyp(1,round((Fmind +length(rhyp))/1.75))+shim,...
    rhyp(2,round((Fmind +length(rhyp))/1.75)),'$$E_2, E_4$$',...
    'HorizontalAlignment','left','Color','r')
text(lhyp(1,round(Oind/3))+shim,lhyp(2,round(Oind/3))+shim,'$$H_2$$',...
    'HorizontalAlignment','left','Color',[85,107,47]/255)
text(lhyp(1,round((F0ind +length(lhyp))/2))+shim,...
    lhyp(2,round((F0ind +length(lhyp))/2)),'$$H_1$$',...
    'HorizontalAlignment','left','Color',[255,140,0]/255)
tmp = round(mean([1,cind]));
text(tmpcirc(1,tmp),tmpcirc(2,tmp),'$$H_1$$','HorizontalAlignment',...
    'left','Color',[255,140,0]/255)
tmp = round(mean([cind,cind2]));
text(tmpcirc(1,tmp),tmpcirc(2,tmp)-shim,'$$E_1$$','HorizontalAlignment',...
    'left','VerticalAlignment','bottom','Color','b')
tmp = round(mean([cind2,cind7]));
text(tmpcirc(1,tmp),tmpcirc(2,tmp),'$$E_2$$','HorizontalAlignment',...
    'center','VerticalAlignment','bottom','Color','r')
tmp = round(mean([100,cind3]));
text(tmpcirc(1,tmp),tmpcirc(2,tmp),'$$H_2$$','HorizontalAlignment',...
    'center','VerticalAlignment','bottom','Color',[85,107,47]/255)

