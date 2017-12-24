function varargout = evalBPC_GUI(varargin)
% EVALBPC_GUI MATLAB code for evalBPC_GUI.fig
%      EVALBPC_GUI, by itself, creates a new EVALBPC_GUI or raises the existing
%      singleton*.
%
%      H = EVALBPC_GUI returns the handle to a new EVALBPC_GUI or the handle to
%      the existing singleton*.
%
%      EVALBPC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVALBPC_GUI.M with the given input arguments.
%
%      EVALBPC_GUI('Property','Value',...) creates a new EVALBPC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evalBPC_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evalBPC_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help evalBPC_GUI

% Last Modified by GUIDE v2.5 21-Sep-2017 05:56:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @evalBPC_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @evalBPC_GUI_OutputFcn, ...
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


% --- Executes just before evalBPC_GUI is made visible.
function evalBPC_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evalBPC_GUI (see VARARGIN)

% Choose default command line output for evalBPC_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes evalBPC_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = evalBPC_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in folderset.
function folderset_Callback(hObject, eventdata, handles)
% hObject    handle to folderset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir();
set(handles.folderdisp,'string',folder);
folderpath = strcat(folder,'/');
fnames = dir(strcat(folderpath,'*.csv'));
samplelist = cell(size(fnames));
for i = 1:length(fnames)
    samplelist{i} = strcat('(',num2str(i),') ',fnames(i).name(1:end-4));
end
set(handles.samplelist,'string',samplelist);


% --- Executes on selection change in samplelist.
function samplelist_Callback(hObject, eventdata, handles)
% hObject    handle to samplelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns samplelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from samplelist


% --- Executes during object creation, after setting all properties.
function samplelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samplelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function samporder_Callback(hObject, eventdata, handles)
% hObject    handle to samporder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of samporder as text
%        str2double(get(hObject,'String')) returns contents of samporder as a double


% --- Executes during object creation, after setting all properties.
function samporder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samporder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
samporder = get(handles.samporder,'string');
if ~isempty(samporder)&~strcmp(samporder,'y')&~strcmp(samporder,'Sample order')&~strcmp(samporder,'auto')
    orderlist = str2num(samporder);
    [results,namelist] = evalBPC(strcat(get(handles.folderdisp,'string'),'/'),orderlist);
elseif strcmp(samporder,'auto')
    [results,namelist] = evalBPC(strcat(get(handles.folderdisp,'string'),'/'));
else
    fnames = dir(strcat(get(handles.folderdisp,'string'),'/*.csv'));
    [results,namelist] = evalBPC(strcat(get(handles.folderdisp,'string'),'/'),[1:size(fnames,1)]);
end
set(handles.BPCval,'Data',results(:,:,1));
set(handles.BPCval,'ColumnName',namelist);
set(handles.BPCval,'RowName',namelist);
set(handles.BPCunc,'Data',results(:,:,2));
set(handles.BPCunc,'ColumnName',namelist);
set(handles.BPCunc,'RowName',namelist);
