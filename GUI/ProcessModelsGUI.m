function varargout = ProcessModelsGUI(varargin)
% PROCESSMODELSGUI MATLAB code for ProcessModelsGUI.fig
%      PROCESSMODELSGUI, by itself, creates a new PROCESSMODELSGUI or raises the existing
%      singleton*.
%
%      H = PROCESSMODELSGUI returns the handle to a new PROCESSMODELSGUI or the handle to
%      the existing singleton*.
%
%      PROCESSMODELSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSMODELSGUI.M with the given input arguments.
%
%      PROCESSMODELSGUI('Property','Value',...) creates a new PROCESSMODELSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcessModelsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcessModelsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcessModelsGUI

% Last Modified by GUIDE v2.5 09-Oct-2013 15:38:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessModelsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessModelsGUI_OutputFcn, ...
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


% --- Executes just before ProcessModelsGUI is made visible.
function ProcessModelsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcessModelsGUI (see VARARGIN)

% Choose default command line output for ProcessModelsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProcessModelsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProcessModelsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
Xname = get(handles.Xname,'String');
Mname = get(handles.Mname,'String');
Yname = get(handles.Yname,'String');
switch popup_sel_index
    case 1
    text(0.3,0.3,Xname)
    text(0.5,0.6,Mname)
    text(0.7,0.3,Yname)
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mname_Callback(hObject, eventdata, handles)
% hObject    handle to Mname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mname as text
%        str2double(get(hObject,'String')) returns contents of Mname as a double
Mname = get(hObject,'String');

% --- Executes during object creation, after setting all properties.
function Mname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xname_Callback(hObject, eventdata, handles)
% hObject    handle to Xname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xname as text
%        str2double(get(hObject,'String')) returns contents of Xname as a double
Xname = get(hObject,'String');

% --- Executes during object creation, after setting all properties.
function Xname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Yname_Callback(hObject, eventdata, handles)
% hObject    handle to Yname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Yname as text
%        str2double(get(hObject,'String')) returns contents of Yname as a double
Yname = get(hObject,'String');


% --- Executes during object creation, after setting all properties.
function Yname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Yname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.


handles.Xname = 'X';
handles.Mname = 'M';
handles.Yname = 'Y';


set(handles.Xname, 'String', handles.Xname);
set(handles.Mname, 'String', handles.Mname);
set(handles.Yname, 'String', handles.Yname);
