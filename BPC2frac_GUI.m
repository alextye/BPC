function varargout = BPC2frac_GUI(varargin)
% BPC2FRAC_GUI MATLAB code for BPC2frac_GUI.fig
%      BPC2FRAC_GUI, by itself, creates a new BPC2FRAC_GUI or raises the existing
%      singleton*.
%
%      H = BPC2FRAC_GUI returns the handle to a new BPC2FRAC_GUI or the handle to
%      the existing singleton*.
%
%      BPC2FRAC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BPC2FRAC_GUI.M with the given input arguments.
%
%      BPC2FRAC_GUI('Property','Value',...) creates a new BPC2FRAC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BPC2frac_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BPC2frac_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BPC2frac_GUI

% Last Modified by GUIDE v2.5 21-Sep-2017 08:49:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BPC2frac_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BPC2frac_GUI_OutputFcn, ...
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


% --- Executes just before BPC2frac_GUI is made visible.
function BPC2frac_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BPC2frac_GUI (see VARARGIN)

% Choose default command line output for BPC2frac_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BPC2frac_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BPC2frac_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function BPCval_Callback(hObject, eventdata, handles)
% hObject    handle to BPCval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPCval as text
%        str2double(get(hObject,'String')) returns contents of BPCval as a double


% --- Executes during object creation, after setting all properties.
function BPCval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPCval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BPCunc_Callback(hObject, eventdata, handles)
% hObject    handle to BPCunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPCunc as text
%        str2double(get(hObject,'String')) returns contents of BPCunc as a double


% --- Executes during object creation, after setting all properties.
function BPCunc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPCunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n1_Callback(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n1 as text
%        str2double(get(hObject,'String')) returns contents of n1 as a double


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n2_Callback(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n2 as text
%        str2double(get(hObject,'String')) returns contents of n2 as a double


% --- Executes during object creation, after setting all properties.
function n2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f1_Callback(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1 as text
%        str2double(get(hObject,'String')) returns contents of f1 as a double


% --- Executes during object creation, after setting all properties.
function f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
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
BPCval = str2num(get(handles.BPCval,'string'));
BPCunc = str2num(get(handles.BPCunc,'string'));
n1 = str2num(get(handles.n1,'string'));
n2 = str2num(get(handles.n2,'string'));
f1 = get(handles.f1,'string');
if ~strcmp(f1,'')&size(f1,2)>0
    BPC2frac(BPCval,BPCunc,n1,n2,0,str2num(f1))
else
    BPC2frac(BPCval,BPCunc,n1,n2)
end
