%function to generate nutation damper plots
function nutationDamper(f,IsIt)
% generate plots of spin for perturbed spacecraft with nutation damper
%      
% Inputs:
%   f (int): Figure to plot in
%   IsIt (float) Ratio of symmetry to transverse axis moi

% Copyright (c) 2016 Dmitry Savransky (ds264@cornell.edu)

if ~exist('f','var')
    f = 1;
end

if ~exist('IsIt','var')
    IsIt = 1.5;
end

%set constants
IwIt = 0.06;
IwIs = IwIt/IsIt;
DIt = 0.5;

z0 = [0.2, 2, 0, 0]; % rad/s

[t,z] = ode45(@nutationDamper_eq,[0,400],z0);

%plot
f = figure(f);
clf
set(f,'Position',[100 100 1000 600]);
subplot(3,1,1)
plot(t,z(:,2),'k')
set(gca,'FontName','Times','FontSize',18)
ylabel('$\omega_2$ (rad)','Interpreter','LaTeX')
subplot(3,1,2)
set(gca,'LineStyleOrder','-|--','ColorOrder',[0,0,0;],...
    'NextPlot','replacechildren','FontName','Times','FontSize',18)
plot(t,z(:,[1,3]))
ylabel('(rad)')
legend({'$\omega_1$','$\omega_3$'},'Interpreter','LaTeX',...
    'Location','northeast','FontName','Times','FontSize',18)
subplot(3,1,3)
plot(t,z(:,4),'k')
set(gca,'FontName','Times','FontSize',18)
ylabel('$\Omega$ (rad)','Interpreter','LaTeX')
xlabel('Time (s)')

    function dz = nutationDamper_eq(t,z)
        %z = [omega_1, omega_2, omega_3, Omega]
        w1 = z(1); w2 = z(2); w3 = z(3); O = z(4);
        
        dz = [(IsIt - 1)*w2*w3 - IwIt*O*w2;
            IwIs*O*w1;
            (-(IsIt - 1)*w1*w2 + DIt*O)/(1 - IwIt);
            ((IsIt - 1)*w1*w2 - DIt/IwIt*O)/(1 - IwIt)];
    end
end