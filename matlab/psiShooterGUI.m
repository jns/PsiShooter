function varargout = psiShooterGUI(varargin)
% PSISHOOTERGUI M-file for psiShooterGUI.fig
%      PSISHOOTERGUI, by itself, creates a new PSISHOOTERGUI or raises the existing
%      singleton*.
%
%      H = PSISHOOTERGUI returns the handle to a new PSISHOOTERGUI or the handle to
%      the existing singleton*.
%
%      PSISHOOTERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSISHOOTERGUI.M with the given input arguments.
%
%      PSISHOOTERGUI('Property','Value',...) creates a new PSISHOOTERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before psiShooterGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to psiShooterGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help psiShooterGUI

% Last Modified by GUIDE v2.5 11-May-2009 20:19:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @psiShooterGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @psiShooterGUI_OutputFcn, ...
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

%%

% --- Executes just before psiShooterGUI is made visible.
function psiShooterGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to psiShooterGUI (see VARARGIN)

%Global data set
global potential psi;

% Choose default command line output for psiShooterGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes psiShooterGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = psiShooterGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%


% --- Executes on selection change in loadMenu.
function loadMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns loadMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from loadMenu

loadMenu_index = get(handles.loadMenu, 'Value');
switch loadMenu_index
    case 1
        %temp = get(handles.loadInput, 'String');
        %if isempty(temp)
        set(handles.loadInput, 'String', 'Input File Path');
        %end
    case 2
        set(handles.loadInput, 'String', '');
    case 3
        set(handles.loadInput, 'String', 'Input Function');
    case 4
        set(handles.loadInput, 'String', '');
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function loadMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject,'String',{'Load From File','Parametric Selection',...
    'Load From Function', 'Load Test Potential'});

%%


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.simPlot);%set up to plot on the main plot window
global potential; %access global potential variable

loadMenu_index = get(handles.loadMenu, 'Value');
switch loadMenu_index
    case 1
        currSysMessText = get(handles.systemMessages, 'String');
        filePath = get(handles.loadInput,'String');
        
        %get the file path if none is specified in the text box
        if or(isempty(filePath),strcmp(filePath,'Input File Path'))
            filePath = uigetfile({'*.*'});
            set(handles.loadInput,'String',filePath);
        end
        
        %reformat default string into cell data structure (so we can have
        %dissimilar length strings)
        if ~iscell(currSysMessText)
            temp = currSysMessText;
            clear('currSysMessText');
            currSysMessText{1} = temp;
        end
        
        %Remind the user what they just did (for debug)
        currSysMessText =[{['Load From File: ' filePath]};currSysMessText];
        set(handles.systemMessages, 'String', currSysMessText);
        
        %Get the potential
        [potential,messages] = loadPotentialData(filePath);
        
        if ~isempty(messages)
            currSysMessText =[messages;currSysMessText];
            set(handles.systemMessages, 'String', currSysMessText);
        end
        
        cla;%clear the plot window
        if isempty(potential.x)
            return
        elseif or(length(potential.x)==1,length(potential.y)==1)
            if length(potential.x) == 1
                visualize1D(potential.y,potential.data,'red');
            else
                visualize1D(potential.x,potential.data,'red');
            end
        elseif and(length(potential.x)>1,length(potential.y)>1)
            visualize2D(potential.x,potential.y,potential.data,0)
        end
        
    case 2
        cla;
    case 3
        cla;
    case 4
        cla;
    case 5
        cla;
end

guidata(hObject,handles)

%%


function loadInput_Callback(hObject, eventdata, handles)
% hObject    handle to loadInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadInput as text
%        str2double(get(hObject,'String')) returns contents of loadInput as a double


% --- Executes during object creation, after setting all properties.
function loadInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%


% --- Executes on button press in simulateButton.
function simulateButton_Callback(hObject, eventdata, handles)
% hObject    handle to simulateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%


% --- Executes on selection change in simulateMenu.
function simulateMenu_Callback(hObject, eventdata, handles)
% hObject    handle to simulateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns simulateMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from simulateMenu


% --- Executes during object creation, after setting all properties.
function simulateMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simulateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject,'String',{'1-D Shooting Simulator','2-D Shooting Simulator',...
    'Digitized Potential w/ Transfer Matrix', 'Run Test Simulation'});

%%


function systemMessages_Callback(hObject, eventdata, handles)
% hObject    handle to systemMessages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of systemMessages as text
%        str2double(get(hObject,'String')) returns contents of
%        systemMessages as a double


% --- Executes during object creation, after setting all properties.
function systemMessages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to systemMessages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','gray');
end

%%


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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double



% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%


function visualize1D(varargin)
%function visualize1D(varargin)
%ex: visualize1D([1 2 3 4 5],[1 2 3 2 1],'red');
%input data structure:
% varargin{1+3*n} = X-axis array
% barargin{2+3*n} = Y-axis array
% style{3+3*n} = display style
%where n is the number of items being plotted
%style will pass matlab plot options directly to the plot command.

if mod(nargin,3) ~= 0 && nargin ~= 0;
    fprintf(['Bad input arguments. There should be a multiple \n' ...
        'of three input arguments. There are ' num2str(nargin) '.\n']);
    return;
end

for n = 1:3:nargin
    X(:,ceil(n/3)) = varargin{n};
    Y(:,ceil(n/3)) = varargin{n+1};
    color{ceil(n/3)} = varargin{n+2};
end

%axes(handles.simPlot);
for n = 1:1:nargin/3
    plot(X(:,n),Y(:,n),color{n})
    hold on;
end
hold off
'end';% a good place to stick a breakpoint for playing with the plots.

%%


function visualize2D(varargin)
%function visualize2D(varargin)
%visualize2D([1:100],[1:100],peaks(100),0,[1:20],[1:20],peaks(20)+10,5);
%input data structure:
% varargin{1+3*n} = x basis vector ie: [x1,x2,x3]
% varargin{2+3*n} = y basis vector ie: {y1,y2,y3]
% varargin{3+3*n} = 2D data array ie: [d_11,d_12,d_13;d_21,d_22,d_23]
% level{4+3*n} = plot offset from zero energy (for the wavefunctions)
%where n is the number of items being plotted
%style will pass matlab plot options directly to the plot command.

if mod(nargin,4) ~= 0 && nargin ~= 0;
    fprintf(['Bad input arguments. There should be a multiple \n' ...
        'of four input arguments. There are ' num2str(nargin) '.\n']);
    return;
end

%plot the wavefunctions and potential together.

%axes(handles.simPlot);%Can't select this one of multiple inputs unless we
%pass the handle structure into this function.

for n = 1:4:nargin
    if n == 1 %plot styles set up for the potential (more transparent)
        h(1)=surf(varargin{n},varargin{n+1},varargin{n+2});
        set(h(1),'facealpha',0.35);
        set(h(1),'edgealpha',0.00);
        hold on;
    else %plot styles set up for wavefunctions
        %h(n) is the handle for the graphic objects I am creating. They can
        %be accessed by get or set(h(n)) with whatever property you want to
        %mess with.
        h(n)=surf(varargin{n},varargin{n+1},varargin{n+2}+varargin{n+3});
        set(h(n),'facealpha',0.7);
        set(h(n),'edgealpha',0.7);
    end
end
hold off;
'end';% a good place to stick a breakpoint for playing with the plots.
%make sure to turn 'rotate3d on;'


%%

function [potential,messages] = loadPotentialData(path)
%function potential = loadPotentialData(path)
%potential.x = [x1,x2,x3,x4,...]
%potential.y = [y1,y2,y3,y4,...]
%potential.data = [d_11,d_12,d_13,d_14;d_21,d_22,d_23,d24;...]
%Works with binary or delimited ascii/unicode data. Save the data file as a
%.txt file if you want to read it in as ascii or unicode. Otherwise, it
%will default to binary.

try
    fidP = fopen(path);
catch
    potential = [];
    messages = {'FILE FAILED TO LOAD'};
    return
end
if fidP == -1
    potential = [];
    messages = {'BAD FILE PATH'};
    return
end

try
    messages = [];
    if strcmp(path(end-3:end),'.txt')
        %ascii or unicode
        xNum = fscanf(fidP,'%d',1);
        yNum = fscanf(fidP,'%d',1);
        potential.x = fscanf(fidP,'%d',xNum);
        potential.y = fscanf(fidP,'%d',yNum);
        potential.data = fscanf(fidP,'%d',[xNum,yNum]);
        if size(potential.data,1)*size(potential.data,2) ~= xNum*yNum
            messages =[{'POTENTIAL DATA MISLOADED (size != xNum*yNum)'};...
                {'Is your file formatted correctly?'};...
                {'It might not like how your decimal numbers are formatted'}];
        end
    else
        %binary
        xNum = fread(fidP,1,'flat64');
        yNum = fread(fidP,1,'flat64');
        potential.x = fread(fidP,xNum,'float64');
        potential.y = fread(fidP,yNum,'float64');
        potential.data = fread(fidP,[xNum,yNum],'float64');
        if size(potential.data,1)*size(potential.data,2) ~= xNum*yNum
            messages =[{'POTENTIAL DATA MISLOADED (size != xNum*yNum)'};...
                {'Is your file formatted correctly?'}];
        end
    end
    
    fclose(fidP);
catch
    potential = [];
    messages = {'POTENTIAL DATA FAILED TO LOAD'};
    fclose(fidP);
    return
end

