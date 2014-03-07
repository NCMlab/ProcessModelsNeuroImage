function varargout = Work_030514(varargin)
% WORK_030514 MATLAB code for Work_030514.fig
%      WORK_030514, by itself, creates a new WORK_030514 or raises the existing
%      singleton*.
%
%      H = WORK_030514 returns the handle to a new WORK_030514 or the handle to
%      the existing singleton*.
%
%      WORK_030514('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORK_030514.M with the given input arguments.
%
%      WORK_030514('Property','Value',...) creates a new WORK_030514 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Work_030514_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Work_030514_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Work_030514

% Last Modified by GUIDE v2.5 06-Mar-2014 20:22:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Work_030514_OpeningFcn, ...
                   'gui_OutputFcn',  @Work_030514_OutputFcn, ...
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


% --- Executes just before Work_030514 is made visible.
function Work_030514_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Work_030514 (see VARARGIN)

% Choose default command line output for Work_030514
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Work_030514 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Work_030514_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in InputData.
function InputData_Callback(hObject, eventdata, handles)
% hObject    handle to InputData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns InputData contents as cell array
%        contents{get(hObject,'Value')} returns selected item from InputData


% --- Executes during object creation, after setting all properties.
function InputData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DataSelector.
function DataSelector_Callback(hObject, eventdata, handles)
% hObject    handle to DataSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DataSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DataSelector
fprintf(1,'Hello\n');
Values = get(hObject,'String');
SelectionValue = get(hObject,'Value');
Selection = Values{SelectionValue};
switch Selection
    case  'Load data from text file'
        fprintf(1,'Loading data from text file\n');
    case  'Load images'
        fprintf(1,'Loading images\n');
        P = spm_select(Inf,'images');
        Name = inputdlg('Please enter a name for this variable');
    case  'Load data from workspace'
        fprintf(1,'Loading workspace variables\n');
end
InputDataList = get(handles.InputData,'String');

Nvar = length(InputDataList);
if Nvar>0
    InputDataList{Nvar+1} = Name{1};
    ColNames = get(handles.Direct,'ColumnName');
    RowNames = get(handles.Direct,'RowName');
    ColNames{Nvar+1} = Name{1};
    RowNames{Nvar+1} = Name{1};
    set(handles.Direct,'ColumnName',ColNames);
    set(handles.Direct,'RowName',RowNames);
    set(handles.Direct,'Data',num2cell(logical(zeros(Nvar+1,Nvar+1))));
    set(handles.Direct,'ColumnEditable',logical(ones(1,Nvar+1)));
    set(handles.Interactions,'ColumnName',ColNames);
    set(handles.Interactions,'RowName',RowNames);
    set(handles.Interactions,'Data',num2cell(logical(zeros(Nvar+1,Nvar+1))));
    set(handles.Interactions,'ColumnEditable',logical(ones(1,Nvar+1)));
    set(handles.Paths,'ColumnName',ColNames);
    set(handles.Paths,'RowName',RowNames);
    
    PathSteps = {};
    for i = 1:Nvar+1
        temp = {'0'};
        for j = 1:Nvar+1
            temp{j+1} = num2str(j);
        end
        PathSteps{i} = temp;
    end
    PathSteps
    set(handles.Paths,'ColumnFormat',PathSteps);
    set(handles.Paths,'ColumnEditable',num2cell(logical(ones(Nvar+1,1))));
  %  set(handles.Paths,'Data',PathSteps);
    ColFormat = get(handles.Direct,'ColumnFormat');
    ColFormat{Nvar+1} = 'logical';
    set(handles.Direct,'ColumnFormat',ColFormat);
    set(handles.Interactions,'ColumnFormat',ColFormat);
    % Create possible pathways 
   

   
else
    InputDataList = Name;
    set(handles.Direct,'Data',num2cell(logical(0)));
    set(handles.Direct,'ColumnName',Name);
    set(handles.Direct,'RowName',Name);
    set(handles.Direct,'ColumnFormat',{'logical'})
    set(handles.Interactions,'Data',num2cell(logical(0)));
    set(handles.Interactions,'ColumnName',Name);
    set(handles.Interactions,'RowName',Name);
    set(handles.Interactions,'ColumnFormat',{'logical'})
   % set(handles.Paths,'Data',num2cell(logical(0)));
    set(handles.Paths,'ColumnName',Name);
    set(handles.Paths,'RowName',Name);
    PathSteps = {};
    for i = 1:Nvar+1
        temp = {'0'};
        for j = 1:Nvar+1
            temp{j+1} = num2str(j);
        end
        PathSteps{i} = temp;
    end
    PathSteps
    
    set(handles.Paths,'ColumnFormat',PathSteps)
    set(handles.Paths,'ColumnEditable',logical(1))
    set(handles.Direct,'ColumnEditable',logical(1));
end
set(handles.InputData,'String',InputDataList)


% --- Executes during object creation, after setting all properties.
function DataSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
