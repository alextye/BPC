function varargout = makePME_GUI(varargin)
% MAKEPME_GUI MATLAB code for makePME_GUI.fig
%      MAKEPME_GUI, by itself, creates a new MAKEPME_GUI or raises the existing
%      singleton*.
%
%      H = MAKEPME_GUI returns the handle to a new MAKEPME_GUI or the handle to
%      the existing singleton*.
%
%      MAKEPME_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKEPME_GUI.M with the given input arguments.
%
%      MAKEPME_GUI('Property','Value',...) creates a new MAKEPME_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before makePME_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to makePME_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help makePME_GUI

% Last Modified by GUIDE v2.5 09-Nov-2017 14:27:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @makePME_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @makePME_GUI_OutputFcn, ...
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


% --- Executes just before makePME_GUI is made visible.
function makePME_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to makePME_GUI (see VARARGIN)

% Choose default command line output for makePME_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes makePME_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = makePME_GUI_OutputFcn(hObject, eventdata, handles) 
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



function loage_Callback(hObject, eventdata, handles)
% hObject    handle to loage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loage as text
%        str2double(get(hObject,'String')) returns contents of loage as a double


% --- Executes during object creation, after setting all properties.
function loage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upage_Callback(hObject, eventdata, handles)
% hObject    handle to upage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upage as text
%        str2double(get(hObject,'String')) returns contents of upage as a double


% --- Executes during object creation, after setting all properties.
function upage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upage (see GCBO)
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
cores = get(handles.cores,'string');
if ~isempty(cores)
    corespec = str2num(cores);
else
    corespec = 0;
end
makePME(strcat(get(handles.folderdisp,'string'),'/'),str2double(get(handles.loage,'string')),str2double(get(handles.upage,'string')),corespec)


function cores_Callback(hObject, eventdata, handles)
% hObject    handle to cores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cores as text
%        str2double(get(hObject,'String')) returns contents of cores as a double


% --- Executes during object creation, after setting all properties.
function cores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
