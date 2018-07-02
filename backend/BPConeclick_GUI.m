function varargout = BPConeclick_GUI(varargin)
% BPCONECLICK_GUI MATLAB code for BPConeclick_GUI.fig
%      BPCONECLICK_GUI, by itself, creates a new BPCONECLICK_GUI or raises the existing
%      singleton*.
%
%      H = BPCONECLICK_GUI returns the handle to a new BPCONECLICK_GUI or the handle to
%      the existing singleton*.
%
%      BPCONECLICK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BPCONECLICK_GUI.M with the given input arguments.
%
%      BPCONECLICK_GUI('Property','Value',...) creates a new BPCONECLICK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BPConeclick_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BPConeclick_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BPConeclick_GUI

% Last Modified by GUIDE v2.5 14-Jun-2018 12:27:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BPConeclick_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BPConeclick_GUI_OutputFcn, ...
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


% --- Executes just before BPConeclick_GUI is made visible.
function BPConeclick_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BPConeclick_GUI (see VARARGIN)

% Choose default command line output for BPConeclick_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BPConeclick_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BPConeclick_GUI_OutputFcn(hObject, eventdata, handles) 
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
cand_folder = uigetdir();
if cand_folder~=0
    folder = cand_folder();
    set(handles.folderdisp,'string',folder);
    folderpath = strcat(folder,'/');
    fnames = dir(strcat(folderpath,'*.csv'));
    samplelist = cell(size(fnames));
    for i = 1:length(fnames)
        samplelist{i} = strcat('(',num2str(i),') ',fnames(i).name(1:end-4));
    end
    set(handles.samplelist,'string',samplelist);
end

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

cores = get(handles.cores,'string');
if ~isempty(cores)
    corespec = str2num(cores);
else
    corespec = 0;
end

auflag = get(handles.aucheck,'value');
BPCuncswitch = get(handles.BPCuncswitch,'value');
quickmode = get(handles.quickmode,'value');
overwrite = get(handles.overwrite,'value');
maxmodels = get(handles.maxmodels,'string');

%if user desires to use quick mode, prevent BPC uncertainties from being
%calculated, as they are unreliable.
if quickmode & BPCuncswitch
    BPCuncswitch = 0;
    msgbox('BPC uncertainties will not be calculated when quick mode is in use.');
end

if ~isempty(maxmodels)
    maxmodelspec = str2num(maxmodels);
else
    maxmodelspec = 0;
end

%hard minimum on number of models retained.
if maxmodelspec > 0 & maxmodelspec < 1000
    maxmodelspec = 1000;
end

%write run settings out to a log file
out_settings = fopen(strcat(get(handles.folderdisp,'string'),'/','BPCrunsettingslog.txt'),'a');
fprintf(out_settings, strcat(datestr(now),'\n'));
fprintf(out_settings, 'directory = %s/\n',get(handles.folderdisp,'string'));
fprintf(out_settings, 'min age = %s Ma\n',get(handles.loage,'string'));
fprintf(out_settings, 'max age = %s Ma\n',get(handles.upage,'string'));
fprintf(out_settings, 'Account for analytical uncertainties = %s\n',num2str(auflag));
fprintf(out_settings, 'Calculate BPC uncertainties = %s\n',num2str(BPCuncswitch));
fprintf(out_settings, 'Quick mode = %s\n',num2str(quickmode));
fprintf(out_settings, 'Overwrite existing PME and BPC unc files = %s\n',num2str(overwrite));
fprintf(out_settings, 'Max models (0 indicates all) = %s\n',num2str(maxmodelspec));
fprintf(out_settings, 'Number of cores (0 indicates not specified) = %s\n',num2str(corespec));

%if either makePME or BPCunc is canceled by user, the subsequent analysis
%and plotting will not occur.

tPMEstart = tic;
%canceledPME flags whether the user canceled the operation during execution
canceledPME = makePME(strcat(get(handles.folderdisp,'string'),'/'),str2double(get(handles.loage,'string')),str2double(get(handles.upage,'string')),corespec,auflag,quickmode,overwrite,maxmodelspec);
tmakePME = toc(tPMEstart);

%do not run the second part of the operation if the user canceled the first
%part
if ~canceledPME & BPCuncswitch
    tBPCuncstart = tic;
    canceledunc = BPCunc(strcat(get(handles.folderdisp,'string'),'/'),str2double(get(handles.loage,'string')),str2double(get(handles.upage,'string')),corespec,overwrite);
    tBPCunc = toc(tBPCuncstart);
else
    canceledunc = 0;
end

fprintf(out_settings, 'PME inferrence runtime = %s seconds\n',num2str(tmakePME));
fprintf(out_settings, 'PME inferrence cancelled during run = %s\n',num2str(canceledPME));
if exist('tBPCunc')
    fprintf(out_settings, 'Calculate BPC uncertainty runtime = %s seconds\n',num2str(tBPCunc));
    fprintf(out_settings, 'Calculate BPC uncertainty cancelled during run = %s\n',num2str(canceledunc));
end

fprintf(out_settings,'***\n\n\n');

if ~canceledPME & ~canceledunc
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

    %output BPC data to tables
    BPCtableout = cell(length(namelist)+1);
    BPCtableout{1,1} = 'Sample';
    BPCtableout([2:end],1) = namelist;
    BPCtableout(1,[2:end]) = namelist';
    BPCtableout([2:end],[2:end]) = mat2cell(results(:,:,1),ones(1,length(namelist)),ones(1,length(namelist)));
    set(handles.BPCval,'Data',BPCtableout);
    BPCunctableout = cell(length(namelist)+1);
    BPCunctableout{1,1} = 'Sample';
    BPCunctableout([2:end],1) = namelist;
    BPCunctableout(1,[2:end]) = namelist';
    BPCunctableout([2:end],[2:end]) = mat2cell(results(:,:,2),ones(1,length(namelist)),ones(1,length(namelist)));
    set(handles.BPCunc,'Data',BPCunctableout);
end

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


% --- Executes on button press in aucheck.
function aucheck_Callback(hObject, eventdata, handles)
% hObject    handle to aucheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aucheck


% --- Executes on button press in BPCuncswitch.
function BPCuncswitch_Callback(hObject, eventdata, handles)
% hObject    handle to BPCuncswitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPCuncswitch


% --- Executes on button press in quickmode.
function quickmode_Callback(hObject, eventdata, handles)
% hObject    handle to quickmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of quickmode


% --- Executes on button press in overwrite.
function overwrite_Callback(hObject, eventdata, handles)
% hObject    handle to overwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overwrite



function maxmodels_Callback(hObject, eventdata, handles)
% hObject    handle to maxmodels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxmodels as text
%        str2double(get(hObject,'String')) returns contents of maxmodels as a double


% --- Executes during object creation, after setting all properties.
function maxmodels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxmodels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
