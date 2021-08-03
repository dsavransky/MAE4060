%Generate the Mars Porkchop Plot for the 2020 Transfer Opportunity

%Copyright (c) 2020 Dmitry Savransky (ds264@cornell.edu)

%% Set up general position/velocity functions
Afun = @(a,e,I,omega,Omega) [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega));...
    a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega));...
    a.*sin(I).*sin(omega)];

Bfun = @(a,e,I,omega,Omega) [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + ...
    sin(Omega).*cos(I).*cos(omega));...
    a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + ...
    cos(Omega).*cos(I).*cos(omega));...
    a.*sqrt(1-e.^2).*sin(I).*cos(omega)];

Mfun = @(n,t,t_p) mod(n*(t - t_p),2*pi);

rfun = @(A,B,E,e) A*(cos(E) - e) + B*sin(E);
vfun = @(A,B,E,e,n) (n/(1 - e*cos(E)))*(-A*sin(E) + B*cos(E));

%% Define orbital elements for Earth and Mars
a_e =  1.000373836656026E+00 ; %AU
a_m = 1.523627165631874E+00;
e_e = 1.712127710968187E-02;
e_m = 9.350045592043167E-02;
I_e = 2.777040607882003E-03*pi/180; %rad
I_m = 1.848059606867434E+00*pi/180; 
w_e = 3.043573249748720E+02*pi/180;
w_m = 2.866772838275621E+02*pi/180;
O_e = 1.596967974767415E+02*pi/180;
O_m = 4.950122742914472E+01*pi/180;
t_p_e = 2458853.731945450883; %JD
t_p_m = 2459064.988962875679; 
mu_e = 2.9591309705483544E-04; %AU^3/day^2
mu_m = 2.9591230378107436E-04; 

% Mean motions
n_e = sqrt(mu_e/a_e^3);
n_m = sqrt(mu_m/a_m^3);

A_e = Afun(a_e,e_e,I_e,w_e,O_e);
B_e = Bfun(a_e,e_e,I_e,w_e,O_e);
A_m = Afun(a_m,e_m,I_m,w_m,O_m);
B_m = Bfun(a_m,e_m,I_m,w_m,O_m);

EarthE = @(t) invKepler(Mfun(n_e,t,t_p_e),e_e);
MarsE = @(t) invKepler(Mfun(n_m,t,t_p_m),e_m);

Earthr = @(t) rfun(A_e,B_e,EarthE(t),e_e);
Earthv = @(t) vfun(A_e,B_e,EarthE(t),e_e,n_e);

Marsr = @(t) rfun(A_m,B_m,MarsE(t),e_m);
Marsv = @(t) vfun(A_m,B_m,MarsE(t),e_m,n_m);

mu = 2.9591220828559093e-4; %sun Gm in AU^3/day^2
kmAU = 149597870.700; %1 AU in km 
R_e = 6378.1366; %Earth radius (km)
mue = 3.986004418e5; %km^3/s^2 Earth Gm
%% check 1 - full orbits

t0 = juliandate(datetime('2020-06-01'));

dts = 0:1:365*3;
r_es = zeros(3,length(dts));
r_ms = zeros(3,length(dts));

for j = 1:length(dts)
    r_es(:,j) = Earthr(t0+dts(j));
    r_ms(:,j) = Marsr(t0+dts(j));
end

figure(1)
clf
plot3(r_es(1,:),r_es(2,:),r_es(3,:),r_ms(1,:),r_ms(2,:),r_ms(3,:))

%% check 2 - reference trajectory
t1 = juliandate(datetime('2020-08-03'));
t2 = juliandate(datetime('2021-02-11'));

%t1 = juliandate(datetime('2020-09-04'));
%t2 = juliandate(datetime('2020-12-09'));
r1 = Earthr(t1);
v1 = Earthv(t1);
r2 = Marsr(t2);
v2 = Marsv(t2);

[v0l,vl] = lambertUniversal(r1,r2,t2-t1,1,mu);
[v0s,vs] = lambertUniversal(r1,r2,t2-t1,-1,mu);

orbdt = linspace(0,t2-t1,100);
transorbl = zeros(3,length(orbdt));
transorbs = zeros(3,length(orbdt));
transorbs(:,1) = r1;
transorbl(:,1) = r1;
xl = [r1;v0l];
xs = [r1;v0s];
for k = 2:length(orbdt)
    try
        xl = keplerSTM(xl,orbdt(k)-orbdt(k-1),mu)*xl;
        xs = keplerSTM(xs,orbdt(k)-orbdt(k-1),mu)*xs;
        transorbs(:,k) = xs(1:3);
        transorbl(:,k) = xl(1:3);
    catch
        transorbs(:,:) = nan;
        transorbl(:,:) = nan;
        disp('oops')
        break
    end
end

figure(1)
clf
plot3(r_es(1,:),r_es(2,:),r_es(3,:),r_ms(1,:),r_ms(2,:),r_ms(3,:))
hold on
plot3(transorbl(1,:),transorbl(2,:),transorbl(3,:),transorbs(1,:),transorbs(2,:),transorbs(3,:))
plot3(r1(1),r1(2),r1(3),'.',r2(1),r2(2),r2(3),'.','MarkerSize',20)

C3l0 = norm((v0l - v1)*kmAU/86400)^2 %geocentric velocity km/s
C3s0 = norm((v0s - v1)*kmAU/86400)^2

vinfl0 = norm((vl - v2)*kmAU/86400) % mars rel vel km/s
vinfs0 = norm((vs - v2)*kmAU/86400)

%% Now lets populate the whole porkchop
dtf = juliandate(datetime('2020-10-01'));
tdep = t0:1:dtf;
tarr0 = juliandate(datetime('2020-12-01'));
tarr = tarr0:1:juliandate(datetime('2021-11-01'));
[d,r] = meshgrid(tdep,tarr);
dts = r-d;

C3l = zeros(length(tdep),length(tarr));
C3s = zeros(length(tdep),length(tarr));
vinfl = zeros(length(tdep),length(tarr));
vinfs = zeros(length(tdep),length(tarr));
for j = 1:length(tdep)
    disp(j)
    r1 = Earthr(tdep(j));
    v1 = Earthv(tdep(j));
    
    for k = 1:length(tarr)
        
        r2 = Marsr(tarr(k));
        v2 = Marsv(tarr(k));
        
        try
            [v0l,vl] = lambertUniversal(r1,r2,tarr(k)-tdep(j),1,mu);
            C3l(j,k) = (norm(v0l - v1)*kmAU/86400)^2;
            vinfl(j,k) = (norm(vl - v2)*kmAU/86400)^2;
        catch
            C3l(j,k) = nan;
            vinfl(j,k) = nan;
            disp('could not find long way')
        end
        try
            [v0s,vs] = lambertUniversal(r1,r2,tarr(k)-tdep(j),-1,mu);
            C3s(j,k) = (norm(v0s - v1)*kmAU/86400)^2;
            vinfs(j,k) = (norm(vs - v2)*kmAU/86400)^2;
        catch
            C3s(j,k) = nan;
            vinfs(j,k) = nan;
            disp('could not find short way')
        end
    end
end

%% C3 plot
nmonthsdep = calmonths(caldiff([datetime(tdep(1),'ConvertFrom','JD'),datetime(tdep(end),'ConvertFrom','JD')]));
xticks = [juliandate(datetime(tdep(1),'ConvertFrom','JD')+caldays(14))];
for j = 1:nmonthsdep-1
    xticks(j+1) = juliandate(datetime(tdep(1),'ConvertFrom','JD')+caldays(14)+calmonths(j));
end
xticklabels = {};
for j = 1:length(xticks)
    xticklabels{j} = char(datetime(xticks(j),'ConvertFrom','JD','Format','dd-MMM'));
end

nmonthsarr = calmonths(caldiff([datetime(tarr(1),'ConvertFrom','JD'),datetime(tarr(end),'ConvertFrom','JD')]));
yticks = [juliandate(datetime(tarr(1),'ConvertFrom','JD')+caldays(14))];
for j = 1:nmonthsarr-1
    yticks(j+1) = juliandate(datetime(tarr(1),'ConvertFrom','JD')+caldays(14)+calmonths(j));
end
yticklabels = {};
for j = 1:length(yticks)
    yticklabels{j} = char(datetime(yticks(j),'ConvertFrom','JD','Format','dd-MMM'));
end

lvls = [ceil(min(C3l(:))*100)/100,ceil(min(C3l(:))):2:20,25,30,50,500];
lvls2 = [ceil(min(C3s(:))*100)/100,ceil(min(C3s(:))):2:20,25,30,100];

figure(2)
clf
[C1,h1] = contour(tdep,tarr,C3l.',lvls,'Linewidth',2);
clabel(C1,h1,'Color','b','FontSize',16)
hold on
[C2,h2] = contour(tdep,tarr,C3s.',lvls2,'Linewidth',2);
clabel(C2,h2,'Color','r','FontSize',16)
caxis([round(min(C3l(:)),2),40])
contour(tdep,tarr,dts,'k--','ShowText','on')

set(gca,'FontName','Times','FontSize',16,'XTick',xticks,...
    'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)
xlabel('Earth Departure (2020)')
ylabel('Mars Arrival (2020-2021)')
title('Departure $C_3$ (km$^2$ s$^{-2}$)','Interpreter','LaTeX')

%% vinf plot

lvls = [ceil(min(vinfl(:))*100)/100,ceil(min(vinfl(:)))+1:2:11,15,25,50];
lvls2 = [ceil(min(vinfs(:))*100)/100,ceil(min(vinfs(:))):2:14,17,100];

figure(3)
clf
[C1,h1] = contour(tdep,tarr,vinfl.',lvls,'Linewidth',2);
clabel(C1,h1,'Color','b','FontSize',16)
hold on
[C2,h2] = contour(tdep,tarr,vinfs.',lvls2,'Linewidth',2);
clabel(C2,h2,'Color','r','FontSize',16)
caxis([round(min(C3l(:)),2),40])
contour(tdep,tarr,dts,'k--','ShowText','on')

set(gca,'FontName','Times','FontSize',16,'XTick',xticks,...
    'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)
xlabel('Earth Departure (2020)')
ylabel('Mars Arrival (2020-2021)')
title('Arrival $v_\infty$ (km s$^{-1}$)','Interpreter','LaTeX')

%% type 2 trajectories

[v,j] = min(C3s);
[vv,k] = min(v);
j = j(k);

t2 = tarr(k);
r2 = Marsr(t2);

t1s = tdep(j)-60:20:tdep(j)+40;
transorbs = zeros(3,100,length(t1s));
r1s = zeros(3,length(t1s));

for j = 1:length(t1s)
    t1 = t1s(j);
    r1 = Earthr(t1);
    r1s(:,j) = r1;
    [v1,v2] = lambertUniversal(r1,r2,t2-t1,-1,mu);
    
    orbdt = linspace(0,t2-t1,100);
    
    transorbs(:,1,j) = r1;
    x = [r1;v1];
    for k = 2:length(orbdt)
        x = keplerSTM(x,orbdt(k)-orbdt(k-1),mu)*x;
        transorbs(:,k,j) = x(1:3);
    end
end

%%
figure(4)
clf
plot3(r_es(1,:),r_es(2,:),r_es(3,:),'b',r_ms(1,:),r_ms(2,:),r_ms(3,:),'r',...
    r2(1),r2(2),r2(3),'r.','MarkerSize',20)
hold on
for j = 1:length(t1s)
    plot3(transorbs(1,:,j),transorbs(2,:,j),transorbs(3,:,j),'k')
    if j == 4
        plot3(transorbs(1,:,j),transorbs(2,:,j),transorbs(3,:,j),'m','Linewidth',4)
    end
    plot3(r1s(1,j),r1s(2,j),r1s(3,j),'b.','MarkerSize',20)
end
view(2)
axis equal
set(gca,'Visible','off')

%% type 1 trajectories

[v,j] = min(C3l);
[vv,k] = min(v);
j = j(k);

t2 = tarr(k);
r2 = Marsr(t2);

t1s = tdep(j)-20:20:tdep(j)+100;
transorbs = zeros(3,100,length(t1s));
r1s = zeros(3,length(t1s));

for j = 1:length(t1s)
    t1 = t1s(j);
    r1 = Earthr(t1);
    r1s(:,j) = r1;
    [v1,v2] = lambertUniversal(r1,r2,t2-t1,1,mu);
    
    orbdt = linspace(0,t2-t1,100);
    
    transorbs(:,1,j) = r1;
    x = [r1;v1];
    for k = 2:length(orbdt)
        x = keplerSTM(x,orbdt(k)-orbdt(k-1),mu)*x;
        transorbs(:,k,j) = x(1:3);
    end
end

%%
figure(5)
clf
plot3(r_es(1,:),r_es(2,:),r_es(3,:),'b',r_ms(1,:),r_ms(2,:),r_ms(3,:),'r',...
    r2(1),r2(2),r2(3),'r.','MarkerSize',20)
hold on
for j = 1:length(t1s)
    plot3(transorbs(1,:,j),transorbs(2,:,j),transorbs(3,:,j),'k')
    if j == 2
        plot3(transorbs(1,:,j),transorbs(2,:,j),transorbs(3,:,j),'m','Linewidth',4)
    end
    plot3(r1s(1,j),r1s(2,j),r1s(3,j),'b.','MarkerSize',20)
end
view(2)
axis equal
set(gca,'Visible','off')