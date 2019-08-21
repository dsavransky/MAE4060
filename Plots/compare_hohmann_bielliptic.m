% Compare Hohmann and Bi-elliptic transfers

% Copyright (c) 2016 Dmitry Savransky (ds264@cornell.edu)

eta = logspace(-1.5,2,1000);

heff = abs(sqrt(2*eta./(1+eta)) - 1 + sqrt(1./eta) - sqrt(2./eta./(1+eta)));
xiinf = sqrt(2) - 1+abs(sqrt(1./eta)- sqrt(2./eta));
xicalc = @(xi) abs(sqrt(2*xi./(1+xi)) - 1) + ...
    abs(sqrt(2*eta./xi./(eta+xi)) - sqrt(2./xi./(1+xi)))+...
    abs(sqrt(1./eta) - sqrt(2*xi./eta./(eta+xi)));
figure(1)
clf
plot(eta,heff,'k','LineWidth',2)
hold on
plot(eta,xiinf,'k--','LineWidth',3)

xis = [2,15.58171,50,100];
legvals = {'Hohmann','$\xi=\infty$'};
for j=1:length(xis)
    plot(eta,xicalc(xis(j)),'LineWidth',2)
    legvals{j+2} = ['$\xi=',num2str(xis(j)),'$'];
end

plot([11.93876,11.93876],[0.2,10],'k--')
plot(1./[11.93876,11.93876],[0.2,10],'k--')
plot([15.58171,15.58171],[0.2,10],'k--')

set(gca,'FontName','Times','FontSize',18,'XScale','log','YScale','log')
xlim([min(eta),max(eta)])
ylim([0.25,max(heff)])
xlabel('$\eta$','Interpreter','Latex')
ylabel('$\Delta v/v_i$','Interpreter','Latex')
legend(legvals,'Interpreter','Latex','Location','Best')



