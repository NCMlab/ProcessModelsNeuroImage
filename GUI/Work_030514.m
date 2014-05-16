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

% Last Modified by GUIDE v2.5 28-Apr-2014 15:15:14

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

Values = get(hObject,'String');
SelectionValue = get(hObject,'Value');
Selection = Values{SelectionValue};
% Determine if the data structure has been initialized
if ~isfield(handles,'data')
        handles.data = {};
        guidata(hObject,handles);
end
CurrentNSub = str2num(get(handles.NSub,'String'));
switch Selection
    case  'Load data from text file'
        % Check to see what the current number of subjects is set to
        
        % select a text file of data having one or more variables
        P = spm_select(1,'txt','Select text file');
        data = textread(P);
        [m n] = size(data);
        
        if m < n
            % rotate so that rows are "subjects"
            data = data';
        end
        % find out how many subjects are in the file
        [NSub ~] = size(data);
        if isempty(CurrentNSub) || CurrentNSub == NSub
            % update the GUI front panel
            set(handles.NSub,'String',num2str(NSub));
            Name = inputdlg('Please enter a name for this variable');
            Ndata = length(handles.data);
            % Add this data to the data structure held by the GUI
            handles.data{Ndata + 1} = data;
            fprintf(1,'Loading data from text file\n');
        else
            % Input data has different number of subjects as existing data
            errordlg('Input data does not have the correct number of subjects');
        end
    case  'Load images'
        fprintf(1,'Loading images\n');
        P = spm_select(Inf,'image');
        % Read the header information 
        V = spm_vol(P);
        NSub = length(V);
        if isempty(CurrentNSub) || CurrentNSub == NSub
            % read the data
            data = spm_read_vols(V);
            Name = inputdlg('Please enter a name for this variable');
            Ndata = length(handles.data);
            % Add this data to the data structure held by the GUI
            handles.data{Ndata + 1} = data;
            
            clear data
        else
            errordlg('Input data does not have the correct number of subjects');
        end
           
    case  'Load data from workspace'
        fprintf(1,'Loading workspace variables\n');
end
guidata(hObject,handles);
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


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Is a mask file loaded?
if isempty(get(handles.Nvoxels,'String'))
    errordlg('Select a mask image');
end
% are there any paths specified?
if isempty(handles.PathData)
    errordlg('Enter at least one path');
else
    % are they all non empty
    for i = 1:size(handles.PathData,3)-1
        if isempty(find(handles.PathData(:,:,i)))
            errordlg('Enter values into the path');
        end
    end
end
% Is an output folder specified?
if isempty(get(handles.BaseFolder,'String'))
    errordlg('Specify a basefolder for the model.');
end
% Is an output name specified?
if isempty(get(handles.ModelName,'String'))
    errordlg('Specify a name for the model.');
end
% Check the data for neuroimaging data and apply the mask
for i = 1:size(handles.data,2);
    [m n o p ] = size(handles.data{i});
    if p > 1
        % found 4-D data
        data = zeros(handles.NVox,p);
        for j = 1:p
            tempData = handles.data{i}(:,:,:,j);
            data(:,j) = tempData(handles.Indices);
        end
        handles.data{i} = data;
    end
end
    
Model = {};
Model.data = handles.data;
Model.Direct = cell2mat(get(handles.Direct,'Data'));
Model.Inter = cell2mat(get(handles.Interactions,'Data'));
Model.Paths = handles.PathData;
Model.names = get(handles.InputData,'String');
Model.BaseFolder = get(handles.BaseFolder,'String');
Model.MaskImage = get(handles.MaskImage,'String');
Model.ModelName = get(handles.ModelName,'String');
% Check to see of the model folder exists
OutPath = fullfile(get(handles.BaseFolder,'String'),get(handles.ModelName,'String'));
if exist(OutPath)
    warndlg('Folder exists, overwriting contents');
    rmdir(OutPath)
end
    mkdir(OutPath)
    Str = sprintf('save %s Model',fullfile(OutPath,'Model'));
    eval(Str);

fprintf(1,'Save Pressed\n');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ModelName_Callback(hObject, eventdata, handles)
% hObject    handle to ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ModelName as text
%        str2double(get(hObject,'String')) returns contents of ModelName as a double


% --- Executes during object creation, after setting all properties.
function ModelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModelName (see GCBO)
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
Vm = spm_vol(P);
mask = spm_read_vols(Vm);
U = unique(mask);
if length(U) ~= 2
    errordlg('Selected mask is not binary.');
else
    % save the included voxel indices
    handles.Indices = find(mask);
    % save the number of voxels in the mask
    handles.NVox = length(handles.Indices);
    % save the header structure of the mask
    Vm.descrip = '';
    Vm.fname = '';
    handles.V = Vm;
    
    set(handles.Nvoxels,'String',num2str(handles.NVox));
    guidata(hObject,handles);
    set(handles.MaskImage,'String',P)
    
end



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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Save.
function Save_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
