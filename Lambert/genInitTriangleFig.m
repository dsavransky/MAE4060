function genInitTriangleFig(r1v,r2v,f)

msize = 20;
figure(f)
hold on
plot(0,0,'k.','MarkerSize',msize)
set(gca,'FontName','Times','FontSize',18,'Visible','Off')
set(gca,'DefaulttextFontName','Times','DefaulttextFontSize',...
    18,'DefaulttextInterpreter','LaTeX','Box','on')
hold on
plot(r1v(1),r1v(2),'k.','MarkerSize',msize)
plot([0,r1v(1)],[0,r1v(2)],'k')
plot(r2v(1),r2v(2),'k.','MarkerSize',msize)
plot([0,r2v(1)],[0,r2v(2)],'k')
plot([r1v(1),r2v(1)],[r2v(2),r2v(2)],'k')
axis equal
hold off

%annotate
shim = max(diff(axis))*0.02/3;
text(0, 0,'$$F$$','HorizontalAlignment','left','VerticalAlignment','top')
text(r1v(1)-shim,r1v(2),'$$P_1$$','HorizontalAlignment','right','VerticalAlignment','bottom')
text(r2v(1),r2v(2),'$$P_2$$','VerticalAlignment','bottom')