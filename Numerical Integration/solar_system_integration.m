%% 2) Solar system
% read in solar system data
% order: Mercury - Pluto, Sun
data = dlmread('solar_system_data.txt','\t',1,1);

% grab everything except pluto
mus = data([1:8,10],1); %gravtiational params (AU^3/day^2)
p0 = reshape(data([1:8,10],2:4).',3*length(mus),1); %initial positions
v0 = reshape(data([1:8,10],5:7).',3*length(mus),1); %initial velocities

%% 2a) integrate - with 10 day time step recording 5000 data points
tic;[T,Y,DY] = do_nbody_leapfrog(0:2000:365.25e5,p0,v0,10,mus);toc

%% a) Plot the orbits
cs = colormap(lines);
figure(1) %inner planets and sun
set(1,'Position', [252, 697, 1000 500])
clf
hold on
for j = [1:4,9], plot(0,0,'color',cs(j,:),'visible','off'); end
for j = [1:4,9]
    plot3(Y(:,(j-1)*3+1),Y(:,(j-1)*3+2),Y(:,(j-1)*3+3),'.','color',cs(j,:))
end
hold off
axis equal
set(gca,'FontName','Times','FontSize',16)
grid on
view(3)
axis equal
legend({'Mercury','Venus','Earth','Mars','Sol'})
xlabel('$\hat{\mathbf{e}}_1$ (AU)','Interpreter','Latex')
ylabel('$\hat{\mathbf{e}}_2$ (AU)','Interpreter','Latex')
zlabel('$\hat{\mathbf{e}}_3$ (AU)','Interpreter','Latex')

figure(2) %outer planets and sun
set(2,'Position', [252, 697, 1000 500])
clf
hold on
for j = 5:9, plot(0,0,'color',cs(j,:),'visible','off'); end
for j = 5:9
    plot3(Y(:,(j-1)*3+1),Y(:,(j-1)*3+2),Y(:,(j-1)*3+3),'.','color',cs(j,:))
end
hold off
axis equal
set(gca,'FontName','Times','FontSize',16)
grid on
view(3)
axis equal
legend({'Jupiter','Saturn','Neptune','Uranus','Sol'})
xlabel('$\hat{\mathbf{e}}_1$ (AU)','Interpreter','Latex')
ylabel('$\hat{\mathbf{e}}_2$ (AU)','Interpreter','Latex')
zlabel('$\hat{\mathbf{e}}_3$ (AU)','Interpreter','Latex')

%% b) Total energy
%energy = 1/2 v^2 + sum_j>i mu_j/|r_ij|
n = length(mus);
inds = reshape(1:3*n,3,n);
cols = repmat(1:n,1,n);
ks = cols(logical(ones(n) - eye(n)));
js = reshape(repmat(1:n,n-1,1),1,n*(n-1));
inds2 = reshape(1:n*(n-1),n-1,n).';

Etot = zeros(length(Y),1);

G = 6.67259e-11;%m^3/kg/s^2
mAU = 1.495978706910000e11;%AU in meters
G = G/mAU^3*86400^2; %AU^3/kg/day^2

for j = 1:length(Etot)
    x = Y(j,:).';
    dx = DY(j,:).';
    
    % r_k/j = r_k/o - r_j/o
    rkj = x(inds(:,ks)) - x(inds(:,js));
    
    Vs = ((mus(ks)/G.*mus(js)).')./sqrt(sum(rkj.^2));
    Ts = sum(reshape(dx,3,n).^2).*mus.'/G/2;
    
    Etot(j) = -sum(Vs)/2 + sum(Ts);
end

figure(3)
clf
plot(T,Etot)
set(gca,'FontName','Times','FontSize',16)
xlabel('Time (days)','Interpreter','Latex')
ylabel('Energy (kg AU$^2$/day$^2$','Interpreter','Latex')
%% c) Earth's osculating orbital elements
mu = mus(3) + mus(end); %total gravitational parameter
r = Y(:,(3-1)*3+1:(3-1)*3+3) - Y(:,end-2:end); %heliocentric r,v
v = DY(:,(3-1)*3+1:(3-1)*3+3) - DY(:,end-2:end);

[a,e,E,I,omega,Omega] = vec2orbElem(r.',v.',mu);


figure(4)
set(4,'Position',[1000,900,1000,1200])
clf
subplot(5,1,1)
plot(T,a)
set(gca,'FontName','Times','FontSize',16)
xlim([min(T),max(T)])
ylim([min(a),max(a)])
ylabel('a (AU)')
subplot(5,1,2)
plot(T,e)
set(gca,'FontName','Times','FontSize',16)
xlim([min(T),max(T)])
ylim([min(e),max(e)])
ylabel('e')
subplot(5,1,3)
plot(T,I)
set(gca,'FontName','Times','FontSize',16)
xlim([min(T),max(T)])
ylim([min(I),max(I)])
ylabel('I (rad)')
subplot(5,1,4)
plot(T,omega)
set(gca,'FontName','Times','FontSize',16)
xlim([min(T),max(T)])
ylim([-pi,pi])
ylabel('\omega (rad)')
subplot(5,1,5)
plot(T,Omega)
set(gca,'FontName','Times','FontSize',16)
xlim([min(T),max(T)])
ylim([-pi,pi])
ylabel('\Omega (rad)')
xlabel('Time (days)')

%%
%look for periodicities
figure(5)
periodogram(e)
title('')
set(gca,'FontName','Times','FontSize',16)
ylabel('Power/Frequency','Interpreter','LateX')
xlabel('Normalized Frequency (x$\pi$ rad/sample)','Interpreter','LateX')

figure(6)
periodogram(omega)
title('')
set(gca,'FontName','Times','FontSize',16)
ylabel('Power/Frequency','Interpreter','LateX')
xlabel('Normalized Frequency (x$\pi$ rad/sample)','Interpreter','LateX')