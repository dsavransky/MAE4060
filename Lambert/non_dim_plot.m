%Lambert Problem Non-Dimensionalized Solution Plot, based on Kaplan (1976)

%Copyright (c) 2021 Dmitry Savransky (ds264@cornell.edu)

%% Time functions
alph = @(s,a) 2*asin(sqrt(s/2./a));
bet = @(s,a,c) 2*asin(sqrt((s-c)/2./a));
gam = @(s,a) 2*asinh(sqrt(s/2./-a));
del  = @(s,a,c) 2*asinh(sqrt((s-c)/2./-a));

E1t = @(s,a,c) a.^(3/2).*(alph(s,a) - sin(alph(s,a)) - (bet(s,a,c) - sin(bet(s,a,c))));
E4t = @(s,a,c) a.^(3/2).*(alph(s,a) - sin(alph(s,a)) + (bet(s,a,c) - sin(bet(s,a,c))));

E2t = @(s,a,c) 2*pi*a.^(3/2) - a.^(3/2).*(alph(s,a) - sin(alph(s,a)) + (bet(s,a,c) - sin(bet(s,a,c))));
E3t = @(s,a,c) 2*pi*a.^(3/2) - a.^(3/2).*(alph(s,a) - sin(alph(s,a)) - (bet(s,a,c) - sin(bet(s,a,c))));

P1t = @(s,c) 1/3*sqrt(2)*(s.^(3/2) - (s - c).^(3/2));
P2t = @(s,c) 1/3*sqrt(2)*(s.^(3/2) + (s - c).^(3/2));

H1t = @(s,a,c) (-a).^(3/2).*(sinh(gam(s,a)) - gam(s,a) - (sinh(del(s,a,c)) - del(s,a,c)));
H2t = @(s,a,c) (-a).^(3/2).*(sinh(gam(s,a)) - gam(s,a) + (sinh(del(s,a,c)) - del(s,a,c)));

%%
mu = 1; %all canonical units
amin = 1;
s = 2*amin;
K = [0,0.9];
Estar = linspace(-1,1,200);

c = s*(1-K);
a = -amin./Estar;

tkern = @(r,a) r./sqrt(2*r - r.^2./a);

ts = zeros(length(a),length(c));
for j = 1:length(a)
    for k = 1:length(c)
        ts(j,k) = integral(@(r) tkern(r,a(j)),s - c(k),s);
    end
end

%%

Estar = linspace(-1,-1e-6,200);
a = -amin./Estar;

%c = s*(1-K) so K = 0 -> c = s, K=1 -> c = 0, K=0.5 -> c = s/2
E10K = E1t(s,a,s);
E41K = E4t(s,a,0);
E11K = E1t(s,a,0);
E15K = E1t(s,a,s/2);
E45K = E4t(s,a,s/2);

E20K = E2t(s,a,s);
E21K = E2t(s,a,0);
E31K = E3t(s,a,0);

E25K = E2t(s,a,s/2);
E35K = E3t(s,a,s/2);

Estar2 = linspace(1e-6,0.5,100);
a = -amin./Estar2;
H10K = H1t(s,a,s);
H11K = H1t(s,a,0);
H21K = H2t(s,a,0);

H15K = H1t(s,a,s/2);
H25K = H2t(s,a,s/2);

%%
figure(1)
clf
hold on
plot(-1,-1,'k','Linewidth',2)
plot(-1,-1,'k--')
plot(-1,-1,'k:','Linewidth',2)
fill([E10K,flip(E41K)],[Estar,flip(Estar)],'r','FaceAlpha',0.5)
fill([E10K,flip(E11K)],[Estar,flip(Estar)],'b','FaceAlpha',0.5)
fill([E20K,flip(E21K)],[Estar,flip(Estar)],'r','FaceAlpha',0.5)
fill([E20K,flip(E31K)],[Estar,flip(Estar)],'b','FaceAlpha',0.5)

fill([H10K,flip(H11K)],[Estar2,flip(Estar2)],[255,140,0]/255,'FaceAlpha',0.5)
fill([H10K,flip(H21K)],[Estar2,flip(Estar2)],[85,107,47]/255,'FaceAlpha',0.5)

plot(H10K,Estar2,'k',E10K,Estar,'k',E20K,Estar,'k','Linewidth',2)

plot(H11K,Estar2,'k:',H21K,Estar2,'k:',E11K,Estar,'k:',E21K,Estar,'k:',E31K,Estar,'k:',E41K,Estar,'k:','Linewidth',2)

plot(E15K,Estar,'k--')
plot(E45K,Estar,'k--')
plot(E25K,Estar,'k--')
plot(E35K,Estar,'k--')
plot(H15K,Estar2,'k--')
plot(H25K,Estar2,'k--')

pimax = 4;
xlim([0,pimax*pi])

plot([0,pimax*pi],[0,0],'k-.')

xticklabs = {'0','\pi'};
for j = 2:pimax
    xticklabs{j+1} = [num2str(j),'\pi'];
end
set(gca,'XTick',(0:1:6)*pi,'XTickLabel',xticklabs)

set(gca,'FontName','Times','FontSize',16)
ylabel('$\mathcal{E}^\star$','FontSize',20,'Interpreter','Latex')
xlabel('$T^\star$','FontSize',20,'Interpreter','Latex')

text(0.1,0.25,'$H_1$','FontSize',20,'Interpreter','Latex')
text(1.75,0.25,'$H_2$','FontSize',20,'Interpreter','Latex')

text(0.4,-0.5,'$E_1$','FontSize',20,'Interpreter','Latex')
text(2.15,-0.5,'$E_4$','FontSize',20,'Interpreter','Latex')


text(5.15,-0.825,'$E_2$','FontSize',20,'Interpreter','Latex')
text(6,-0.92,'$E_3$','FontSize',20,'Interpreter','Latex')

legend({'$K = 0$','$K = 0.5$','$K = 1$'},'FontSize',16,'Interpreter','Latex')
