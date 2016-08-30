% Generate a plot of the effective potential as a function of radius

% Copyright (c) 2016 Dmitry Savransky (ds264@cornell.edu)

%set Umin and h
Umin = 1
h = 1
mu = sqrt(2*h^2*Umin)
r0 = h^2/2/mu

% %set Umin and r0
% Umin = 1
% r0 = 1
% h = sqrt(8*r0^2*Umin)
% mu = h^2/2/r0

rmin = h^2/mu;
r = linspace(r0/10,rmin*5,1000);
ub = h^2/2./r.^2;
lb = -mu./r;
Ueff = h^2/2./r.^2-mu./r;
yax = [-Umin*2.5,Umin*3];


figure(1)
clf
plot(r,Ueff,'Linewidth',2)
hold on
plot(r,ub,'--',r,lb,'--')
m1 = plot(rmin,-Umin,'o','MarkerSize',10);
set(m1, 'MarkerFaceColor', get(m1, 'Color'));
m2 = plot(r0,0,'d','MarkerSize',10);
set(m2, 'MarkerFaceColor', get(m2, 'Color'));

plot([0,max(r)],[0,0],'k--')
plot([0,max(r)],[-Umin,-Umin],'k--')
plot([rmin,rmin],yax,'k--')
plot([r0,r0],yax,'k--')

legend({'$U_\mathrm{eff}$','$h^2/2r^2$','$-\mu/r$','Circular Orbit','Parabolic Orbit'},...
    'Interpreter','Latex','FontSize',18);

uistack(m1,'top')
uistack(m2,'top')

set(gca,'FontName','Times','FontSize',18,'XTick',[0,r0,rmin,max(r)],'XTickLabel',...
    {'$0$','$\frac{h^2}{2\mu}$','$\frac{h^2}{\mu}$','$\infty\rightarrow$'},...
    'YTick',[min(yax),-Umin,0,max(yax)],'YTickLabel',...
    {'$-\infty\downarrow$','$\min U_\textrm{eff}$','$0$','$\infty \uparrow$'},...
    'TickLabelInterpreter', 'Latex')
xlabel('$r$','Interpreter', 'Latex')

text(rmin*1.1,-Umin/4,'Elliptical Orbits (2 turning points)','Interpreter','Latex','FontSize',18)
text(rmin*1.1,max(yax)*0.9,{'Hyperbolic Orbits','(1 turning point)'},...
    'Interpreter','Latex','FontSize',18,'VerticalAlignment','Top')
xlim([0,max(r)])
ylim(yax)

ax1 = gca;
ax2 = axes('Position', ax1.Position,...
    'YAxisLocation','right',...
    'Color','none','XTick',[],'TickLabelInterpreter', 'Latex',...
    'FontName','Times','FontSize',18,'YLim',yax,'YTick',[0],...
    'YTickLabelRotation',90,'YTickLabel',{'$\leftarrow$ Closed $\,\,$   Open $\rightarrow$'});