function varargout = planet_shape(varargin)
% PLANET_SHAPE GUI for exploring spherical harmonics and resulting Earth
% shapes.
%      PLANET_SHAPE, by itself, creates a new PLANET_SHAPE or raises the existing
%      singleton*.
%
%      H = PLANET_SHAPE returns the handle to a new PLANET_SHAPE or the handle to
%      the existing singleton*.
%
%      PLANET_SHAPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLANET_SHAPE.M with the given input arguments.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2016 Dmitry Savransky (ds264@cornell.edu)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @planet_shape_OpeningFcn, ...
                   'gui_OutputFcn',  @planet_shape_OutputFcn, ...
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


% --- Executes just before planet_shape is made visible.
function planet_shape_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to planet_shape (see VARARGIN)

% Choose default command line output for planet_shape
handles.output = hObject;

%create internal data structure
topo = load('topo.mat','topo');
handles.data.phi = linspace(0,pi,50);
handles.data.theta = linspace(0,2*pi,50);
[Theta,Phi] = meshgrid(handles.data.theta,handles.data.phi);
handles.data.Theta = Theta;
handles.data.Phi = Phi;
%handles.data.isPlaying = false;
handles.data.amplitude = 0.1;
handles.data.l = 0;
handles.data.m = 0;

% Update handles structure
guidata(hObject, handles);
set(handles.lval,'String',num2str(handles.data.l));
set(handles.mval,'String',num2str(handles.data.m));
set(handles.amplitude,'Value',handles.data.amplitude);

%generate initial spherical harmonic
[x,z,y] = planet_shape_calc_spherical_harmonic(hObject);
handles.data.x0 = x;
handles.data.y0 = y;
handles.data.z0 = z;

%generate surfaces
cla(handles.ax1,'reset')
cla(handles.ax2,'reset')
props.EdgeColor = 'none';
props.FaceLighting = 'phong';

handles.data.harmonic = surface(handles.ax1,handles.data.x0,...
    handles.data.y0,handles.data.z0,props);
axis(handles.ax1,'tight','equal','off')
view(handles.ax1,3)
camzoom(handles.ax1,1.25)

props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.Cdata = topo.topo;

handles.data.planet_shape = surface(handles.ax2,handles.data.x0,...
    handles.data.y0,handles.data.z0,props);
axis(handles.ax2,'tight','equal','off')
view(handles.ax2,3)
camzoom(handles.ax2,1.25)

% %make both axes rotatable 
handles.ax3d = rotate3d;
set(handles.ax3d,'Enable','on');

% Update handles structure
guidata(hObject, handles);


% function axClickFun(obj,event_obj)
% 
% %rotate3d(get(obj,'Parent'),'ON')
% disp(obj)
% 

function [x,y,z] = planet_shape_calc_spherical_harmonic(hObject)

handles = guidata(hObject);

% Calculate the harmonic
Pmls = legendre(handles.data.l,cos(handles.data.phi));
Pml = Pmls(handles.data.m+1,:).';

U = Pml*(cos(handles.data.m*handles.data.theta)+...
    sin(handles.data.m*handles.data.theta));
r = U/max(abs(U(:)));
x = r.*sin(handles.data.Phi).*cos(handles.data.Theta);
y = r.*sin(handles.data.Phi).*sin(handles.data.Theta);
z = r.*cos(handles.data.Phi);


function planet_shape_update_plots(hObject)

handles = guidata(hObject);
[x,z,y] = planet_shape_calc_spherical_harmonic(hObject);
set(handles.data.harmonic,'XData',x,'YData',y,'ZData',z)
set(handles.data.planet_shape,...
    'XData',handles.data.x0 - x*handles.data.amplitude,...
    'YData',handles.data.y0 - y*handles.data.amplitude,...
    'ZData',handles.data.z0 - z*handles.data.amplitude)


% --- Outputs from this function are returned to the command line.
function varargout = planet_shape_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hObject;


function lval_Callback(hObject, eventdata, handles)
% hObject    handle to lval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject,'String'));
if isnan(value)
    warning('planet_shape:NaNInput','Input must be numeric.');
    set(handles.lval,'String',num2str(handles.data.l));
    return
end
if (value < 0) 
    warning('planet_shape:OOBInput','l must be >= 0.');
    set(handles.lval,'String',num2str(handles.data.l));
    return
end
if (value < handles.data.m) 
    warning('planet_shape:OOBInput','m must be <= l.');
    set(handles.lval,'String',num2str(handles.data.m));
    value = handles.data.m;
end

handles.data.l = value;

% Update handles structure
guidata(hObject, handles);

planet_shape_update_plots(hObject);


% --- Executes during object creation, after setting all properties.
function lval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mval_Callback(hObject, eventdata, handles)
% hObject    handle to mval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject,'String'));
if isnan(value)
    warning('planet_shape:NaNInput','Input must be numeric.');
    set(handles.mval,'String',num2str(handles.data.m));
    return
end
if (value < 0) || (value > handles.data.l)
    warning('planet_shape:OOBInput','m must be >= 0 and < l.');
    set(handles.mval,'String',num2str(handles.data.m));
    return
end
handles.data.m = value;

% Update handles structure
guidata(hObject, handles);

planet_shape_update_plots(hObject);


% --- Executes during object creation, after setting all properties.
function mval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.amplitude = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

planet_shape_update_plots(hObject);


% --- Executes during object creation, after setting all properties.
function amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
