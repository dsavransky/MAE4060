function genInitTriangleFig(r1v,r2v,f)
%Plot the basic geometry of the Lambert problem
%
%INPUT
%   r1v (1x2 float) x and y components of vector from focus (central body)
%                   to origin of trajectory 
%   r2v (1x2 float) x and y components of vector from focus (central body)
%                   to destination of trajectory 
%   f   (int)       Figure number to plot in

%Copyright (c) 2009 Dmitry Savransky (ds264@cornell.edu)
%Updated 2022-05-10

msize = 30;
figure(f)
hold on
plot(0,0,'k.','MarkerSize',msize)
set(gca,'FontName','Times','FontSize',18,'Visible','Off')
set(gca,'DefaulttextFontName','Times','DefaulttextFontSize',...
    18,'DefaulttextInterpreter','LaTeX','Box','on')
hold on
plot(r1v(1),r1v(2),'k.','MarkerSize',msize)
%plot([0,r1v(1)],[0,r1v(2)],'k')
quiver(0,0,r1v(1),r1v(2),0,'k','MaxHeadSize',0.25/norm(r1v),'Linewidth',2)
plot(r2v(1),r2v(2),'k.','MarkerSize',msize)
%plot([0,r2v(1)],[0,r2v(2)],'k')
quiver(0,0,r2v(1),r2v(2),0,'k','MaxHeadSize',0.25/norm(r2v),'Linewidth',2)
plot([r1v(1),r2v(1)],[r2v(2),r2v(2)],'k','Linewidth',2)
axis equal
hold off

%annotate
shim = max(diff(axis))*0.02/3;
text(0, 0,'$$F$$','HorizontalAlignment','left','VerticalAlignment','top')
text(r1v(1)-shim,r1v(2),'$$P_1$$','HorizontalAlignment','right',...
    'VerticalAlignment','bottom')
text(r2v(1),r2v(2),'$$P_2$$','VerticalAlignment','bottom')