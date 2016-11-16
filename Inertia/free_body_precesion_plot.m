It = 1;
Is = 2;
I = [It,It,Is];
ws = 1;
w0 = 0.1;

[t,w] = torque_free_motion([0,10],[w0,0,ws],[It,It,Is]);

figure(1)
clf
plot(t,w)

wn = (Is - It)/It*ws;

hold on
t1 = linspace(t(1),t(end),50);
plot(t1,w0*cos(wn*t1),'o',t1,w0*sin(wn*t1),'o')