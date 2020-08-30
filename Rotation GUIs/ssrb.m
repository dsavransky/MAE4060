function varargout = ssrb(varargin)
% SSRB Spinning Symmetric Rigid Body
%      Animate and explore the motion of a symmetric spinning body with a
%      constant moment applied to one transverse body frame axis.
%      
%      The main axis displays the body, the current orientation of the
%      body symmetry (spin) axis and a trace of the path taken by the spin
%      axis in inertia space.  The axis can be roated at any time in 3
%      dimensions to better explore the visualization.
%
%      The second axis display a trace of the current spin rate (Omega).
%      The texbox shows the value of the initial state vector, where psi is
%      the rotation about the inertial vertical axis (e_3), theta is the
%      rotation about the body transverse axis (b_1), and phi is rotation
%      about the body symmetry axis (such that Omega-the spin rate-is the
%      time derivative of phi). The values are updated when the simulation
%      is stopped to their current values.  All values can be edited by the
%      user to select a new initial integration state.
%
%      The first slider controls the ratio of the the moment of inertia
%      about the non-symmetry body axes (I_1) to the mometn of inertia
%      about the symmetry axis (I_2).  Values of I1/I2 greater than 1 make
%      the body a minor axis spinner.
%
%      The second slider controls the ratio of the moment applied about the
%      transverse body axis (b_1) to the moment of inertia about the
%      symmetry axis.
%
%      All user editable components of the interface (the sliders, textbox,
%      and main axis rotation) can be toggled while the simulation is
%      running.
%
%      SSRB, by itself, creates a new SSRB or raises the existing
%      singleton.
%
%      H = SSRB returns the handle to a new SSRB or the handle to
%      the existing singleton.
%
%      SSRB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSRB.M with the given input arguments.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Written by Dmitry Savransky 02-May-2011
% Copyright (c) 2015 Dmitry Savransky (ds264@cornell.edu)
% Major rewrite 2017

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ssrb_OpeningFcn, ...
                   'gui_OutputFcn',  @ssrb_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end
% End initialization code - DO NOT EDIT


% --- Executes just before ssrb is made visible.
function ssrb_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ssrb (see VARARGIN)

% add file location to path
handles.path0 = path();
s = which('ssrb');
pathstr = fileparts(s);
path(pathstr,path)

% Choose default command line output for ssrb
handles.output = hObject;
    
%set up 3d rotation for main axis
handles.main3d = rotate3d(handles.mainAx);
handles.baseView = [40,25];
if ~isempty(varargin)
    tmp = varargin{1};
    if isnumeric(tmp) && numel(tmp) == 2
        handles.baseView = tmp;
    end
end
set(handles.main3d,'RotateStyle','box','Enable','on');

%rotation matrices
rotMatE3 = @(ang) [cos(ang) -sin(ang) 0;sin(ang) cos(ang) 0;0 0 1];
rotMatE2 = @(ang) [cos(ang) 0 sin(ang);0 1 0;-sin(ang) 0 cos(ang)];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) -sin(ang);0 sin(ang) cos(ang)];
handles.rotMats = {rotMatE1,rotMatE2,rotMatE3};

% Update handles structure & set default vals
guidata(hObject, handles);
ssrb_set_def_vals(hObject,handles)
handles = guidata(hObject);

% This sets up the initial plot - only do when we are invisible
% so window can get raised
if strcmp(get(hObject,'Visible'),'off')
    %fill in state label
    cla(handles.statelabelax);
    text(handles.statelabelax,...
        'String','$\theta,\dot\theta,\psi,\dot\psi,\phi,\Omega=$',...
        'Interpreter','Latex','HorizontalAlignment','left',...
        'VerticalAlignment','bottom','FontSize',18);
    
    %set up plot axis
    cla(handles.plotAx)
    xlabel(handles.plotAx,'Time (s)','FontSize',16)
    ylabel(handles.plotAx,'$\Omega$','Interpreter','Latex','FontSize',16)
    hold(handles.plotAx,'on')
    
    
    ssrb_makeEllipsoid(hObject, handles)
end
end

function ssrb_set_def_vals(hObject,handles)
%set default values

handles.I2 = 1/(10 - 0.5);     %kg m^2
handles.I1 = 10*handles.I2;    %kg m^2
handles.dims = [sqrt(handles.I2/2), sqrt(handles.I1 - handles.I2/2)];
handles.M1 = 5*handles.I2;      %N m
set(handles.MOIslider,'Value',log10(handles.I1/handles.I2));
set(handles.MOItext,'String',sprintf('I1/I2 = %2.3f',handles.I1/handles.I2))
set(handles.M1slider,'Value',handles.M1/handles.I2);
set(handles.M1text,'String',sprintf('M1/I2 = %2.3f',handles.M1/handles.I2))

%default initial state vector [th,thd,psi,psid,phi,Omega] (rad, rad/s)
handles.statevec0 = [pi/6,0,0,-0.1,0,20];

guidata(hObject, handles);

end


% --- Executes when user attempts to close mainWindow.
function figure1_CloseRequestFcn(hObject, ~, handles)

path(handles.path0);
delete(hObject);

end

% --- Outputs from this function are returned to the command line.
function varargout = ssrb_OutputFcn(~, ~, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

end

% --- Executes on slider movement - slider value is I1/I2
function MOIslider_Callback(hObject, ~, handles)

anim = false;
%if animation is going, stop it
if strcmp(get(handles.runAnim,'String'),'Stop')
    set(handles.runAnim,'String','Go')
    anim = true;
end

rat = 10^get(hObject,'Value');
handles.I2 = 1/(rat - 0.5);     %kg m^2
handles.I1 = rat*handles.I2;    %kg m^2
handles.dims = [sqrt(handles.I2/2), sqrt(handles.I1 - handles.I2/2)];
set(handles.MOItext,'String',sprintf('I1/I2 = %2.3f',handles.I1/handles.I2))

% Update handles structure and regenerate ellipse
guidata(hObject, handles);
ssrb_reset_Omegatrace(handles)
ssrb_makeEllipsoid(hObject, handles)
handles = guidata(hObject);

if anim
    set(handles.runAnim,'String','Stop')
    ssrb_intFun(hObject,handles);
end

end

% --- Executes during object creation, after setting all properties.
function MOIslider_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% --- Executes on slider movement.
function M1slider_Callback(hObject, ~, handles)
anim = false;
%if animation is going, stop it
if strcmp(get(handles.runAnim,'String'),'Stop')
    set(handles.runAnim,'String','Go')
    anim = true;
end

handles.M1 = get(hObject,'Value')*handles.I2;
set(handles.M1text,'String',sprintf('M1/I2 = %2.3f',handles.M1/handles.I2))
guidata(hObject, handles);
ssrb_reset_Omegatrace(handles)

if anim
    set(handles.runAnim,'String','Stop')
    ssrb_intFun(hObject,handles);
end

end

% --- Executes during object creation, after setting all properties.
function M1slider_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


function state_Callback(hObject, ~, handles)

[val,count] = sscanf(get(hObject,'String'),'%f, %f, %f, %f, %f, %f');
if count ~= 6
    ssrb_update_displaystate(handles);
    return
else
    %if animation is going, stop it
    anim = false;
    if strcmp(get(handles.runAnim,'String'),'Stop')
        set(handles.runAnim,'String','Go')
        anim = true;
    end
    
    handles.statevec = val.';
    th = handles.statevec(1);
    psi = handles.statevec(3);
    handles.b2cur = handles.rotMats{3}(psi)*handles.rotMats{1}(th)*[0,1,0].';
    handles.b1cur = handles.rotMats{3}(psi)*handles.rotMats{1}(th)*[1,0,0].';
    
    % Update handles structure and regenerate ellipse
    guidata(hObject, handles);
    ssrb_reset_Omegatrace(handles)
    ssrb_makeEllipsoid(hObject, handles)
    handles = guidata(hObject);
    
    if anim
        set(handles.runAnim,'String','Stop')
        ssrb_intFun(hObject,handles);
    end
    
end
end


% --- Executes during object creation, after setting all properties.
function state_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%generate the ellipsoid and traces
function ssrb_makeEllipsoid(hObject, handles)

[xe,ye,ze]=ellipsoid(0,0,0,handles.dims(1),handles.dims(2),handles.dims(1),100);

if ~isfield(handles,'s1') || ~ishandle(handles.s1)
    axes(handles.mainAx)
    handles.statevec = handles.statevec0;

    handles.s1 = surface(handles.mainAx,xe,ye,ze,'FaceAlpha',0.9);
    shading interp
    view(3)
    axis equal
    view(handles.baseView)
    handles.b2cur = [0,1,0].';
    handles.b1cur = [1,0,0].';
    
    th = handles.statevec(1);
    rotate(handles.s1,handles.b1cur,th*180/pi)
    handles.b2cur = handles.rotMats{1}(th)*handles.b2cur;
        
    hold on
    handles.e3 = quiver3(0,0,0,0,0,1.1,0,'Linewidth',2,'Color','k');
    handles.b2 = quiver3(0,0,0,handles.b2cur(1)*1.1,handles.b2cur(2)*1.1,...
        handles.b2cur(3)*1.1,0,'LineWidth',2);
    handles.b2trace = plot3(handles.b2cur(1)*1.1,handles.b2cur(2)*1.1,...
        handles.b2cur(3)*1.1,'-');
    hold off
    
    handles.Omegatrace = plot(handles.plotAx,0,handles.statevec(end));
else
    temp = [handles.b1cur,handles.b2cur,cross(handles.b1cur,handles.b2cur)]*...
        [xe(:).';ye(:).';ze(:).'];
    set(handles.s1,'XData',reshape(temp(1,:),size(xe)),...
        'YData',reshape(temp(2,:),size(ye)),...
        'ZData',reshape(temp(3,:),size(ze)));
    set(handles.b2,'UData',handles.b2cur(1)*1.1,...
        'VData',handles.b2cur(2)*1.1,...
        'WData',handles.b2cur(3)*1.1)
    if ishandle(handles.b2trace), delete(handles.b2trace); end
    hold on
    handles.b2trace = plot3(handles.b2cur(1)*1.1,handles.b2cur(2)*1.1,...
        handles.b2cur(3)*1.1,'-');
    hold off
end
axis([-1,1,-1,1,-1,1]*ceil(max(handles.dims)*11)/10)
grid on

ssrb_update_displaystate(handles);

% Update handles structure
guidata(hObject, handles);

end

%write the current statevec to the textbox
function ssrb_update_displaystate(handles)

tmp = handles.statevec;
for j=[1,3,5]
    tmp(j) = mod(tmp(j),2*pi);
end

set(handles.state,'String',...
    sprintf('%3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f',...
    tmp));

end

% --- Executes on button press in resetFig.
function resetFig_Callback(hObject, ~, handles)

if strcmp(get(handles.runAnim,'String'),'Stop')
    set(handles.runAnim,'String','Go')
end

if ishandle(handles.s1), delete(handles.s1); end
if ishandle(handles.e3), delete(handles.e3); end
if ishandle(handles.b2), delete(handles.b2); end
if ishandle(handles.b2trace), delete(handles.b2trace); end
if ishandle(handles.Omegatrace), delete(handles.Omegatrace); end

ssrb_set_def_vals(hObject,handles)
handles = guidata(hObject);

ssrb_makeEllipsoid(hObject, handles)

end

% --- Executes on button press in runAnim.
function runAnim_Callback(hObject, ~, handles)

if strcmp(get(hObject,'String'),'Go')
    set(hObject,'String','Stop')
    ssrb_intFun(hObject,handles);
else
    set(hObject,'String','Go')
    ssrb_update_displaystate(handles);
    ssrb_reset_Omegatrace(handles)
end

end

%rest the omega trace
function ssrb_reset_Omegatrace(handles)

set(handles.Omegatrace,'XData',0,'YData',handles.statevec(end))

end

%performs numerical integration of eoms
function ssrb_intFun(hObject,handles)

M1 = handles.M1;
I1 = handles.I1;
I2 = handles.I2;
t = 0:0.1:100-0.1;
%z = zeros(length(t),6);

[~,z] = ode45(@ssrb_eom,t,handles.statevec);
handles.t = t;
handles.z = z;
th = z(:,1); psi = z(:,3);
handles.b2hist = [-sin(psi).*cos(th),cos(psi).*cos(th), sin(th)];
handles.b1hist = [cos(psi),sin(psi), zeros(size(psi))];

% Update handles structure
guidata(hObject, handles);

%check if user stopped anim before integration terminated, otherwise
%animate
if strcmp(get(handles.runAnim,'String'),'Go')
    return;
else
    ssrb_animate(hObject,handles)
end

% [~,tmp] = ode45(@ssrb_eom,t(1:10),handles.statevec);
% z(1:10,:) = tmp;
% for j = 2:50
%     drawnow
%     handles = guidata(handles.figure1);
%     if strcmp(get(handles.runAnim,'String'),'Go')
%         return;
%     end
%     inds = (j-1)*10:j*10;
%     [~,tmp] = ode45(@ssrb_eom,t(inds),z(inds(1),:));
%     if length(tmp) == length(inds)
%         z(inds,:) = tmp;
%     else
%         set(handles.runAnim,'String','Go')
%         t = -1;
%         return
%     end
% end

    function dz = ssrb_eom(~,z)
        %z = [th,thd,psi,psid,phi,Omega]
        theta = z(1);
        thetadot = z(2);
        psidot = z(4);
        Omega = z(6);
        
        thetaddot = (-I1.*psidot.^2.*sin(2*theta)/2 + I2.*Omega.*psidot.*cos(theta) +...
            I2.*psidot.^2.*sin(2*theta)/2 - M1)./I1;
        psiddot = thetadot.*(2*I1.*psidot.*tan(theta) - I2.*Omega./cos(theta) - ...
            I2.*psidot.*tan(theta))./I1;
        Omegadot = thetadot.*(-I1.*psidot.*sin(theta).^2 - I1.*psidot + ...
            I2.*Omega.*sin(theta) + I2.*psidot.*sin(theta).^2)./(I1.*cos(theta));
        
        dz = [thetadot;
            thetaddot;
            psidot;
            psiddot;
            Omega;
            Omegadot];
    end
end

%animates ellipsoid
function ssrb_animate(hObject,handles)

t = handles.t;
z = handles.z;

hold on
for j = 2:length(t)
    drawnow
    if strcmp(get(handles.runAnim,'String'),'Go')
        hold off;
        return;
    end
    dth = z(j,1) - z(j-1,1);
    dpsi = z(j,3) - z(j-1,3);
    dt = t(j) - t(j-1);
    dphi = z(j,5)-z(j-1,5);
    
    rotate(handles.s1,[0,0,1],dpsi*180/pi) %psi rotation about e3
    rotate(handles.s1,handles.b1hist(j,:),dth*180/pi) %th rotation about b1
    rotate(handles.s1,handles.b2hist(j,:),dphi*180/pi) %phi rotation about b2
    
%     rotate(handles.s1,[0,0,1],dpsi*180/pi)
%     handles.symAx = handles.rotMats{3}(dpsi)*handles.symAx;
%     handles.b1 = handles.rotMats{3}(dpsi)*handles.b1;
%     rotate(handles.s1,handles.b1,dth*180/pi)
%     handles.symAx = handles.rotMats{1}(dth)*handles.symAx;
%     rotate(handles.s1,handles.symAx,dphi*180/pi)
%     plot3(handles.mainAx,handles.symAx(1)*handles.dims(2),...
%         handles.symAx(2)*handles.dims(2),...
%         handles.symAx(3)*handles.dims(2),'k.')
%     
    set(handles.b2,'UData',handles.b2hist(j,1)*1.1,...
        'VData',handles.b2hist(j,2)*1.1,...
        'WData',handles.b2hist(j,3)*1.1)
    
    set(handles.b2trace,'XData',[get(handles.b2trace,'XData'),handles.b2hist(j,1)*1.1],...
        'YData',[get(handles.b2trace,'YData'),handles.b2hist(j,2)*1.1],...
        'ZData',[get(handles.b2trace,'ZData'),handles.b2hist(j,3)*1.1])
    
    otime = get(handles.Omegatrace,'XData');
    if t(j) < otime(end)
        set(handles.Omegatrace,'XData',0,'YData',z(j,end));
    else
        set(handles.Omegatrace,'XData',[otime,t(j)],...
            'YData',[get(handles.Omegatrace,'YData'),z(j,end)])
    end
    
    handles.statevec = z(j,:);
    handles.b1cur = handles.b1hist(j,:).';
    handles.b2cur = handles.b2hist(j,:).';
    
    % Update handles structure
    guidata(hObject, handles);

    pause(dt/2)
end

%if we haven't cancelled yet, go again
ssrb_intFun(hObject,handles);

end
