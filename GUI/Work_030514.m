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

% Last Modified by GUIDE v2.5 19-Mar-2014 21:48:19

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
% User defined variables
handles.PathData = [];
% Store these new values into the GUI data so they can be retrieved
% elsewhere in the program
guidata(hObject,handles);
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
    set(handles.Paths,'Data',zeros(Nvar+1,Nvar+1))
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
%    set(handles.Paths,'ColumnEditable',num2cell(logical(ones(Nvar+1,1))));
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
    set(handles.Paths,'Data',[0])
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


% --- Executes on selection change in PathSelector.
function PathSelector_Callback(hObject, eventdata, handles)
% hObject    handle to PathSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PathSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PathSelector
SelectedPath = get(handles.PathSelector,'value')
set(handles.Paths,'Data',handles.PathData(:,:,SelectedPath));
%handles.Paths(:,:,SelectedPath) = get(hObject,'Data');

% --- Executes during object creation, after setting all properties.
function PathSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PathSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in Paths.
function Paths_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Paths (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
fprintf(1,'Hello\n');
size(handles.PathData,3)
% Which path is being edited?
SelectedPath = get(handles.PathSelector,'value')
% get the updated data from this path and add it to the handles
handles.PathData(:,:,SelectedPath) = get(hObject,'Data');
guidata(hObject,handles);

% --- Executes on button press in AddPath.
function AddPath_Callback(hObject, eventdata, handles)
% hObject    handle to AddPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf(1,'Hello\n')

Nvar = length(get(handles.InputData,'String'));
NumberOfPaths = str2num(get(handles.NumberOfPaths,'String'))+1;

handles.PathData(:,:,NumberOfPaths+1) = zeros(Nvar,Nvar);
guidata(hObject,handles);

set(handles.NumberOfPaths,'String',num2str(NumberOfPaths))
% Add the new path number to the path selector piece
PathSelectorString = {};
for j = 1:NumberOfPaths
    PathSelectorString{j} = num2str(j);
end
set(handles.Paths,'Data',handles.PathData(:,:,NumberOfPaths+1));
set(handles.PathSelector,'String',PathSelectorString)
set(handles.PathSelector,'value',NumberOfPaths)





% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over AddPath.
function AddPath_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to AddPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectMask.
function SelectMask_Callback(hObject, eventdata, handles)
% hObject    handle to SelectMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
P = spm_select(1,'image','Select Mask image');
set(handles.MaskImage,'String',P)



function MaskImage_Callback(hObject, eventdata, handles)
% hObject    handle to MaskImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaskImage as text
%        str2double(get(hObject,'String')) returns contents of MaskImage as a double


% --- Executes during object creation, after setting all properties.
function MaskImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BaseFolder_Callback(hObject, eventdata, handles)
% hObject    handle to BaseFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaseFolder as text
%        str2double(get(hObject,'String')) returns contents of BaseFolder as a double


% --- Executes during object creation, after setting all properties.
function BaseFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaseFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectBaseFolder.
function SelectBaseFolder_Callback(hObject, eventdata, handles)
% hObject    handle to SelectBaseFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
P = spm_select(1,'dir');
set(handles.BaseFolder,'String',P);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over SelectMask.
function SelectMask_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to SelectMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MaskImage.
function MaskImage_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MaskImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
