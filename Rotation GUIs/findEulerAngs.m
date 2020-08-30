function varargout = findEulerAngs(varargin)
% FINDEULERANGS Interactive Euler Angles Demonstration
%      FINDEULERANGS, by itself, creates a new FINDEULERANGS or raises 
%      the existing singleton.
%
%      H = FINDEULERANGS returns the handle to a new FINDEULERANGS or the
%      handle to the existing singleton.  Data can be accessed via
%      GUIDATA(H).
%
%      H = FINDEULERANGS([az,el]) sets the base figure view (baseView) to
%      [az, el],measured in degrees. Default is (100,15).
%
%      FINDEULERANGS('CALLBACK',hObject,eventData,handles,...) calls the
%      local function named CALLBACK in FINDEULERANGS.M with the given
%      input arguments.
%
%      Use the mouse to change the orientation of the box.  A wireframe of
%      the original box position will appear.  The 'Euler Axis' button
%      animates the rotation about the Euler Axis between the original and
%      new box positions.  The 'Rotate' button calculates a set of Euler
%      angles for the rotation type currently selected ('Body' or 'Space')
%      and the axis set selected in the dropdown menu to orient the axes to
%      the new box position. 'Derotate' reverses this process and brings
%      the axes back to their previous position. 
%      Once the axes have been rotated, subsequent rotations are calculated
%      from this new orientation.  This means that Space rotations use the
%      coordinate system associated with the previous cube orientation
%      (i.e. the axes of the wireframe cube). 
%      'Axes DCM' is the direction cosine matrix of the current axes
%      orientation in the inertial (baseView) frame.  'Box DCM (inertial)'
%      is the DCM of the box's current orientation in the baseView frame.
%      'Box DCM (Axes)' is the matrix product of the previous two DCMs.
%      The Euler axis is calculated in the baseView frame. 
%
%      Note: findEulerAngs adds the folder in which it is placed to the
%      path while it is running to ensure that the m file is always
%      accessible.  The path is restored to its original form on a clean
%      exit.
%      
% See also: GUIDE, GUIDATA, GUIHANDLES

% Written by Dmitry Savransky 21-Apr-2011 dsavrans@princeton.edu
% 5-May-2011 Added Space rotations
% 01-Jan-2012 Offloaded backedn to calcEulerAngs

% Copyright (c) 2015 Dmitry Savransky (ds264@cornell.edu)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @findEulerAngs_OpeningFcn, ...
                   'gui_OutputFcn',  @findEulerAngs_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes just before findEulerAngs is made visible.
function findEulerAngs_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to findEulerAngs (see VARARGIN)

% add file location to path
handles.path0 = path();
s = which('findEulerAngs');
pathstr = fileparts(s);
path(pathstr,path)

% Choose default command line output for findEulerAngs
handles.output = hObject;
handles.main3d = rotate3d(handles.mainAx);
handles.baseView = [100,15];
if ~isempty(varargin)
    tmp = varargin{1};
    if isnumeric(tmp) && numel(tmp) == 2
        handles.baseView = tmp;
    end
end
set(handles.main3d,'RotateStyle','box','Enable','on',...
    'ActionPostCallback',@postRotateFun);

%these aren't backwards, we're just defining the aCi DCM
rotMatE3 = @(ang) [cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
rotMatE2 = @(ang) [cos(ang) 0 -sin(ang);0 1 0;sin(ang) 0 cos(ang)];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];
handles.rotMats = {rotMatE1,rotMatE2,rotMatE3};

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using findEulerAngs.
if strcmp(get(hObject,'Visible'),'off')
    findEulerAngs_make3daxes(hObject, handles)
end

% --- Executes when user attempts to close mainWindow.
function mainWindow_CloseRequestFcn(hObject, ~, handles)

path(handles.path0);
delete(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = findEulerAngs_OutputFcn(~, ~, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in rotStyle.
function rotStyle_Callback(~, ~, ~)
%contents = cellstr(get(hObject,'String')) returns rotStyle contents as cell array
%contents{get(hObject,'Value')} returns selected item from rotStyle

function rotStyle_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset.
function reset_Callback(hObject, ~, handles)

set(handles.doRotation,'String','Rotate')
set(handles.rotMsg,'String','');
findEulerAngs_make3daxes(hObject, handles);

% --- Executes on button press in doRotation.
function doRotation_Callback(hObject, ~, handles)

%if newRot means we derotate, otherwise, we rotate
if handles.newRot
    if sum(handles.prevRot) == 0
        return;
    end
    for j=3:-1:1
        if handles.prevRotType
            %previous rotation was body
            %rotate function takes angles in degrees
            findEulerAngs_rot3daxes(handles,handles.coordSys(handles.prevRot(1,j),:),...
                -handles.prevRot(2,j)*180/pi,1);
            %update coordinate system
            handles.coordSys = handles.rotMats{handles.prevRot(1,j)}(-handles.prevRot(2,j))*handles.coordSys;
        else
            %previous rotaton was space
            findEulerAngs_rot3daxes(handles,handles.coordSys0(handles.prevRot(1,j),:),...
                -handles.prevRot(2,j)*180/pi,1);  
        end
    end
    if ~handles.prevRotType, handles.coordSys = handles.coordSys0; end
    handles.newRot = false;
    set(handles.doRotation,'String','Rotate')
    set(handles.rotMsg,'String','');
else
    %handles.dcm is wrt the orig inertial axes
    %we want to calculate angles with respect to the current axes,
    %in the frame of the current axes (aCi * (bCi iCa))
    dcm = handles.coordSys*handles.dcm*handles.coordSys.';
    
    %which rotStyle
    bodyRot = get(handles.rotTypeBody,'Value');
    rotStyle = get(handles.rotStyle,'Value');
    contents = get(handles.rotStyle,'String');
    if bodyRot, rstring = 'Body'; else rstring = 'Space'; end
    
    switch rotStyle
        case 1 %1-2-3
            rotSet = [1,2,3];
        case 2 %2-3-1
            rotSet = [2,3,1];
        case 3 %3-1-2
            rotSet = [3,1,2];
        case 4 %1-3-2
            rotSet = [1,3,2];
        case 5 %2-1-3
            rotSet = [2,1,3];
        case 6 %3-2-1
            rotSet = [3,2,1];
        case 7 %1-2-1
            rotSet = [1,2,1];
        case 8 %1-3-1
            rotSet = [1,3,1];
        case 9 %2-1-2
            rotSet = [2,1,2];
        case 10 %2-3-2
            rotSet = [2,3,2];
        case 11 %3-1-3
            rotSet = [3,1,3];
        case 12 %3-2-3
            rotSet = [3,2,3];
    end
    [th1,th2,th3] = calcEulerAngs(dcm,rotSet,~bodyRot);
     
    angs = [th1,th2,th3];
    handles.coordSys0 = handles.coordSys;
    %inds = 1:3;
    set(handles.currRotAx,'Visible','on')
    if bodyRot
        for j=1:3
            %set(handles.as(inds(inds ~= rotSet(j))),'FaceAlpha',0.2,'EdgeAlpha',0.2)
            %set(handles.as(rotSet(j)),'FaceAlpha',1,'EdgeAlpha',1)
            
            ca = get(handles.as(rotSet(j)),{'Xdata','Ydata','Zdata'});
            ca = [mean(ca{1},2),mean(ca{2},2),mean(ca{3},2)];
            ca = ca(2,:) - ca(1,:);
            set(handles.currRotAx,'Xdata',[-1,1]*ca(1)*10,...
                'Ydata',[-1,1]*ca(2)*10,'Zdata',[-1,1]*ca(3)*10);
            
            %rotate function takes angles in degrees
            findEulerAngs_rot3daxes(handles,handles.coordSys(rotSet(j),:),...
                angs(j)*180/pi);
            %update coordinate system
            handles.coordSys = handles.rotMats{rotSet(j)}(angs(j))*handles.coordSys;
        end
    else
        for j=1:3
            ca = handles.coordSys(rotSet(j),:);
            set(handles.currRotAx,'Xdata',[-1,1]*ca(1)*10,...
                'Ydata',[-1,1]*ca(2)*10,'Zdata',[-1,1]*ca(3)*10);
            findEulerAngs_rot3daxes(handles,handles.coordSys(rotSet(j),:),...
                angs(j)*180/pi);
        end
        handles.coordSys = handles.coordSys*handles.dcm.';
    end
    set(handles.currRotAx,'Visible','off');
    
    %set(handles.as,'FaceAlpha',1,'EdgeAlpha',1)
    handles.newRot = true;
    handles.prevRot = [rotSet;angs];
    handles.prevRotType = bodyRot;
    set(handles.doRotation,'String','Derotate')
    set(handles.rotMsg,'String',[rstring,' ',contents{rotStyle},...
        ' Rotation angles: ',num2str(angs*180/pi,'%3.3f  '),' degrees.']);
end

% Update handles and display
set(handles.DCMA,'String',findEulerAngs_numPrint(handles.coordSys))
set(handles.DCMB,'String',findEulerAngs_numPrint(handles.coordSys*handles.dcm))

guidata(hObject, handles);


% --- Executes on button press in EulerAxis.
function EulerAxis_Callback(~, ~, handles)

set(handles.rotVec,'Xdata',[-1,1]*handles.rotAx(1)*5,...
                   'Ydata',[-1,1]*handles.rotAx(2)*5,...
                   'Zdata',[-1,1]*handles.rotAx(3)*5,'Visible','on');
nsteps = round(abs(handles.th)/3.1304); %180 deg in 2.5 s at 23 fps
for j=1:nsteps
    rotate(handles.bx,handles.rotAx,-handles.th/nsteps,[0,0,0])
    rotate(handles.circs,handles.rotAx,-handles.th/nsteps,[0,0,0])
    pause(1/23);
end
pause(1);
for j=1:nsteps
    rotate(handles.bx,handles.rotAx,handles.th/nsteps,[0,0,0])
    rotate(handles.circs,handles.rotAx,handles.th/nsteps,[0,0,0])
    pause(1/23);
end
set(handles.rotVec,'Visible','off')


%create or reset a dextral set of coordinate axes
function findEulerAngs_make3daxes(hObject, handles)

%cleanup any existing objects
delete(get(handles.mainAx,'Children'));

%create an axis
cylinderScale = 0.05;
[xrod,yrod,zrod] = cylinder(cylinderScale,12);

hold on

%create and place the axes
as = zeros(1,3);
as(3) = surface(xrod,yrod,zrod,'FaceColor','b');
as(2) = surface(xrod,yrod,zrod,'FaceColor','g');
as(1) = surface(xrod,yrod,zrod,'FaceColor','r');
rotate(as(2),[1 0 0],-90,[0,0,0])
rotate(as(1),[0 1 0],90,[0,0,0])

%define the cube
cubeScale = 0.6;
X =[-1    -1    -1    -1    -1     1;
     1    -1     1     1     1     1;
     1    -1     1     1     1     1;
    -1    -1    -1    -1    -1     1]*cubeScale;

Y =[-1    -1    -1    -1     1    -1;
    -1     1    -1    -1     1     1;
    -1     1     1     1     1     1;
    -1    -1     1     1     1    -1]*cubeScale;

Z =[-1    -1     1    -1    -1    -1;
    -1    -1     1    -1    -1    -1;
     1     1     1    -1     1     1;
     1     1     1    -1     1     1]*cubeScale;

C = [0.5 1 0 0 0.5 1];
bx = fill3(X,Y,Z,C,'FaceAlpha',0.3);
obx = zeros(1,6);

for j = 1:6
    tmp = get(bx(j),{'Xdata','Ydata','Zdata'});
    tmp = cat(2,tmp{:});
    tmp2 = [tmp;tmp(1,:)];
    obx(j) = plot3(tmp2(:,1),tmp2(:,2),tmp2(:,3),'k--');
end

%Euler axis
rotVec = plot3([0,1],[0,1],[0,1],'k','LineWidth',2);
set(rotVec,'Visible','off')

%current rotation axis
currRotAx = plot3([0,1],[0,1],[0,1],'k--','LineWidth',2);
set(currRotAx,'Visible','off')

%previous axes
oas = zeros(1,3);
col = ['r','g','b'];
for j=1:3
    a = get(as(j),{'Xdata','Ydata','Zdata'});
    oas(j) = plot3(mean(a{1},2),mean(a{2},2),mean(a{3},2),...
        col(j),'LineWidth',3);
end

%circles on cube sides
th = linspace(0,2*pi,15).';
circs = zeros(1,3);
circScale = cylinderScale*1.1;
circs(1) = plot3(circScale*cos(th), circScale*sin(th), zeros(size(th))+cubeScale,'k','LineWidth',3);
circs(2) = plot3(zeros(size(th))+cubeScale, circScale*sin(th), circScale*cos(th),'k','LineWidth',3);
circs(3) = plot3(circScale*sin(th), zeros(size(th))+cubeScale, circScale*cos(th),'k','LineWidth',3);
hold off

%clean up display
axis(handles.mainAx);
axis equal;
axis([-1 1 -1 1 -1 1]*1.25)
view(handles.baseView)
grid on

% Update handles structure and display
handles.as = as;
handles.oas = oas;
handles.bx = bx;
handles.obx = obx;
handles.rotVec = rotVec;
handles.currRotAx = currRotAx;
handles.circs = circs;
handles.newRot = true;
handles.coordSys = eye(3);
handles.dcm = eye(3);
handles.rotAx = zeros(1,3);
handles.prevRot = zeros(2,3);
handles.prevRotType = 1; %1 for Body, 0 for space
handles.th = 0;
set(handles.DCMI,'String',num2str(handles.dcm,'%1.3f   '))
set(handles.DCMA,'String',num2str(handles.coordSys,'%1.3f   '))
set(handles.DCMB,'String',num2str(handles.coordSys*handles.dcm,'%1.3f   '))
set(handles.lastRot,'String',sprintf('%3.3f about inertial z\n\n%3.3f about inertial y',[0,0]))
set(handles.EAxis,'String',[num2str(handles.rotAx,'%1.3f   '),sprintf('\n\n%3.3f degrees',handles.th)])
set(handles.rotMsg,'String','')
guidata(hObject, handles);


function postRotateFun(obj,event_obj)
handles = guidata(obj);

%if this is a new rotation
if handles.newRot
    handles.newRot = false;
    handles.dcm = eye(3); %reset dcm
    set(handles.rotMsg,'String','');
    set(handles.doRotation,'String','Rotate')
    
    %set old box to new box
    for j = 1:6
        tmp = get(handles.bx(j),{'Xdata','Ydata','Zdata'});
        tmp = cat(2,tmp{:});
        tmp2 = [tmp;tmp(1,:)];
        set(handles.obx(j),'XData',tmp2(:,1),'YData',tmp2(:,2),...
            'ZData',tmp2(:,3));
    end
    %set old axes to new axes
    for j=1:3
        a = get(handles.as(j),{'Xdata','Ydata','Zdata'});
        a = [mean(a{1},2),mean(a{2},2),mean(a{3},2)];
        a = [0;1]*(a(2,:) - a(1,:));
        set(handles.oas(j),'Xdata',a(:,1),'Ydata',a(:,2),'Zdata',a(:,3));
        
%         set(handles.oas(j),'XData',mean(a{1},2),'YData',mean(a{2},2),...
%             'ZData',mean(a{3},2));
    end
end

%calculate current rotation
newView = get(event_obj.Axes,'View');
dth = handles.baseView - newView;
set(handles.lastRot,'String',sprintf('%3.3f about inertial z\n\n%3.3f about inertial y',dth))
view(handles.baseView) %reset view
dth = dth*pi/180;

%DCM about absolute inertial view (coordSys = I)
dcm = [cos(dth(1))*cos(dth(2)),-sin(dth(1)),-cos(dth(1))*sin(dth(2));
       sin(dth(1))*cos(dth(2)),cos(dth(1)),-sin(dth(1))*sin(dth(2));
       sin(dth(2)), 0 , cos(dth(2))];
for j = 1:6
    newbx = dcm*[get(handles.bx(j),'XData'),...
        get(handles.bx(j),'YData'),...
        get(handles.bx(j),'ZData')].';
    set(handles.bx(j),'XData',newbx(1,:),'YData',newbx(2,:),...
        'ZData',newbx(3,:));
end
for j = 1:3
    newCirc = dcm*[get(handles.circs(j),'XData');...
        get(handles.circs(j),'YData');...
        get(handles.circs(j),'ZData')];
    set(handles.circs(j),'XData',newCirc(1,:),'YData',newCirc(2,:),...
        'ZData',newCirc(3,:));
end

% Euler parameters (again from Kane & Levinson)
dcm = dcm*handles.dcm;
e4 = sqrt(1 + sum(diag(dcm)))/2;
e1 = (dcm(3,2) - dcm(2,3))/(4*e4);
e2 = (dcm(1,3) - dcm(3,1))/(4*e4);
e3 = (dcm(2,1) - dcm(1,2))/(4*e4);

% Update handles structure and display
handles.dcm = dcm;
handles.rotAx = [e1,e2,e3];
handles.th = 2*acos(e4)*180/pi;
set(handles.DCMI,'String',findEulerAngs_numPrint(handles.dcm))
set(handles.DCMB,'String',findEulerAngs_numPrint(handles.coordSys*handles.dcm))
set(handles.EAxis,'String',[findEulerAngs_numPrint(handles.rotAx),sprintf('\n\n%3.3f degrees',handles.th)])
guidata(obj, handles);

%animate rotation
function findEulerAngs_rot3daxes(handles,l,ang,nsteps)

if ~exist('nsteps','var')
    nsteps = round(abs(ang)/3.1304); %180 deg in 2.5 s at 23 fps
end
for j=1:nsteps
    for i=1:3
        rotate(handles.as(i),l,ang/nsteps,[0,0,0])
    end
    if nsteps > 1, pause(1/23);end
end

%format number strings
function out =  findEulerAngs_numPrint(in)
out = num2str(round(in*1000)/1000+1-1,'%1.3f   ');
