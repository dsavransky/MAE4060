% Visualization of the Kinetic Energy Ellipsoid. 

%Copyright (c) 2009 Dmitry Savransky (ds264@cornell.edu)

% set up figure
if ishandle(256), close(256); end
figure(256)

%define principal moments of inertia and KE
I = [1 2 3];
T = 1;

%make the ellipsoid
[xe,ye,ze]=ellipsoid(0,0,0,sqrt(2*T/I(1)),sqrt(2*T/I(2)),...
    sqrt(2*T/I(3)),100);
s1 = surface(xe,ye,ze,'FaceAlpha',0.8,'SpecularExponent',10,...
    'SpecularStrength',0.3);
hold on
%set up lighting
l(1) = light('Position',[0 1.5 0]);
%l(2) = light('Position',[1.5 0.5 -0.5]);
lighting phong
shading interp
cs = colormap('gray');
%cs = cs(30:50,:); unnecessary in newer MATLAB (after 2020a)
colormap(cs)
view(3)
axis equal
view(-200,28)
pc = [0.3,0.3,0.3];
plot3([0 2],[0 0],[0,0],'Linewidth',2,'Color',pc)
plot3([0 0],[0 2],[0,0],'Linewidth',2,'Color',pc)
plot3([0 0],[0 0],[0,1.5],'Linewidth',2,'Color',pc)
set(gca,'Visible','off')
text(1.9,0,0,'I_1','FontName','Times','FontSize',14)
text(0.05,2.18,0,'I_2','FontName','Times','FontSize',14)
text(0.02,0,1.55,'I_3','FontName','Times','FontSize',14)

%define angular momenta of polhodes - one at the critical point, two above
%and two below
cp = sqrt(2*T*I(2));
minL = sqrt(2*T*I(1));
maxL = sqrt(2*T*I(3));
Ls = [minL*1.001,linspace(minL*1.01,cp,4),...
    linspace(cp,maxL*0.99,4),maxL*0.9999];

for j = 1:length(Ls)
    L = Ls(j);
    if L < sqrt(2*T*I(2))
        psi0 = acos((1/I(1) + 1/I(2) - 4*T/L^2)/(1/I(1) - 1/I(2)))/2;
        psi = [psi0:pi/1000:pi - psi0; pi+psi0:pi/1000:2*pi - psi0];
    else
        psi = 0:pi/1000:2*pi;
    end
    
    temp = ((2*T/L^2) - 1/I(3))./...
        ((cos(psi).^2/I(2) + sin(psi).^2/I(1)) - 1/I(3));
    sth = sqrt(temp);
    cth = sqrt(1 - temp);
    %clean up MATLAB's numerical garbage:
    if any(imag(cth)~=1), cth = real(cth); end
    
    %find cartesian coordinates (because they're easy to plot)
    x = -L/I(1)*sth.*sin(psi);
    y = L/I(2)*sth.*cos(psi);
    z = -L/I(3)*cth;
    
    
    %plot!
    if L < sqrt(2*T*I(2))
        for k=1:2
            plot3([x(k,:),x(k,end:-1:1)],[y(k,:),y(k,end:-1:1)],...
                [z(k,:),-z(k,end:-1:1)],'k','Linewidth',2)
        end
    else
        plot3(x,y,z,'k','Linewidth',2)
        plot3(x,y,-z,'k','Linewidth',2)
    end
end