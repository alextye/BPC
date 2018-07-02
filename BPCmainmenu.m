function varargout = BPCmainmenu(varargin)
% BPCMAINMENU MATLAB code for BPCmainmenu.fig
%      BPCMAINMENU, by itself, creates a new BPCMAINMENU or raises the existing
%      singleton*.
%
%      H = BPCMAINMENU returns the handle to a new BPCMAINMENU or the handle to
%      the existing singleton*.
%
%      BPCMAINMENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BPCMAINMENU.M with the given input arguments.
%
%      BPCMAINMENU('Property','Value',...) creates a new BPCMAINMENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BPCmainmenu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BPCmainmenu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BPCmainmenu

% Last Modified by GUIDE v2.5 15-Feb-2018 15:12:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BPCmainmenu_OpeningFcn, ...
                   'gui_OutputFcn',  @BPCmainmenu_OutputFcn, ...
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


% --- Executes just before BPCmainmenu is made visible.
function BPCmainmenu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BPCmainmenu (see VARARGIN)

% Choose default command line output for BPCmainmenu
handles.output = hObject;

%dummy = spmak([0 0 0 0.5 1 1 1], [1 1 1 1 1]);

license_names = {'Optimization Toolbox\n', 'Statistics and Machine Learning Toolbox\n', 'Curve Fitting Toolbox\n', 'Parallel Computing Toolbox\n', 'Global Optimization Toolbox\n'};
installed = zeros(1,5);
installed(1) = ~isempty(ver('optim'));
installed(2) = ~isempty(ver('stats'));
installed(3) = ~isempty(ver('curvefit'));
installed(4) = ~isempty(ver('distcomp'));
installed(5) = ~isempty(ver('globaloptim'));
mess = sprintf(strcat('The following toolboxes are necessary and not installed:\n\n',license_names{find(installed==0)}));
if sum(installed) < 5
    msgbox(mess);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BPCmainmenu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BPCmainmenu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BPConeclick.
function BPConeclick_Callback(hObject, eventdata, handles)
% hObject    handle to BPConeclick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('backend');
BPConeclick_GUI

% --- Executes on button press in PMEplot.
function PMEplot_Callback(hObject, eventdata, handles)
% hObject    handle to PMEplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('backend');
PMEplot_GUI

% --- Executes on button press in PME2CSV.
function PME2CSV_Callback(hObject, eventdata, handles)
% hObject    handle to PME2CSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('backend');
PME2CSV_GUI

% --- Executes on button press in BPC2frac.
function BPC2frac_Callback(hObject, eventdata, handles)
% hObject    handle to BPC2frac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('backend');
BPC2frac_GUI
