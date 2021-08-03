% Plots and Animations for the Circular Restricted 3-Body Problem

% Copyright (c) 2014 Dmitry Savransky

%% Hill Curves and Surface
mu = 0.3;
x = linspace(-1.5,1.5,1000);
[x,y] = meshgrid(x,x);
r1 = sqrt((x+mu).^2+y.^2);
r2 = sqrt((x - (1-mu)).^2+y.^2);
U = -(x.^2+y.^2)/2 - ((1-mu)./r1 + mu./r2);
tmp = -log(abs(U)); %visualize in log stretch
%% Curves
lvls = [ceil(min(tmp(:))),ceil(min(tmp(:)))+[2,3],ceil(min(tmp(:)))+3.5:0.25:-1,-0.9:0.05:-0.75,linspace(-0.7,max(tmp(:)),12)];
figure(1)
clf
[C,H] = contour(x,y,tmp,lvls,'Linewidth',2);
caxis([-0.7,max(tmp(:))])
axis equal
set(gca,'Xtick',[],'YTick',[])

%% Surface
figure(2)
clf
surface(x,y,tmp)
shading flat; 
caxis([-0.7,max(tmp(:))])
zlim([-0.7,max(tmp(:))])
set(gca,'Visible','off')
view(72,66)
%% Animation
figure(3)
clf
frames = logspace(log10(abs(max(U(:)))),log10(abs(min(U(:)))/100),200);
for j=length(frames):-1:2
    buh = zeros(size(U)); 
    buh(U > -frames(j)) = max(frames) - frames(j)+min(frames);
    imagesc(buh,[min(frames),max(frames)]);
    set(gca,'XTick',[],'YTick',[])
    pause(0.01); 
end

%% Lagrange Points
f = @(x) x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
L3 = fsolve(f,-mu-0.1);
L2 = fsolve(f,1-mu+0.1);
L1 = fsolve(f,0.1);
C4 = -1/2*(3 - mu + mu^2);
Cl123 = @(x0) -(x0.^2)/2 - ((1-mu)./sqrt((x0+mu).^2) + mu./sqrt((x0 - (1-mu)).^2));

figure(2)
hold on
plot3(L1,0,-log(abs(Cl123(L1))),'k.','MarkerSize',20)
plot3(L2,0,-log(abs(Cl123(L2))),'k.','MarkerSize',20)
plot3(L3,0,-log(abs(Cl123(L3))),'k.','MarkerSize',20)
plot3([1/2-mu,1/2-mu],[sqrt(3)/2,-sqrt(3)/2],zeros(1,2)-log(abs(C4)),'k.','MarkerSize',20)

t1 = text(L1,0,-log(abs(Cl123(L1))),'$L_1$','FontSize',32,'Interpreter',...
    'Latex','VerticalAlignment','bottom','HorizontalAlignment','center');
t2 = text(L2,0,-log(abs(Cl123(L2))),'$L_2$','FontSize',32,'Interpreter',...
    'Latex','VerticalAlignment','top','HorizontalAlignment','left');
t3 = text(L3,0,-log(abs(Cl123(L3))),'$L_3$','FontSize',32,'Interpreter',...
    'Latex','VerticalAlignment','top','HorizontalAlignment','left');
t4 = text(1/2-mu,-sqrt(3)/2,-log(abs(C4)),'$L_4$','FontSize',32,'Interpreter',...
    'Latex','VerticalAlignment','bottom','HorizontalAlignment','right');
t5 = text(1/2-mu,sqrt(3)/2,-log(abs(C4)),'$L_5$','FontSize',32,'Interpreter',...
    'Latex','VerticalAlignment','bottom','HorizontalAlignment','left');
