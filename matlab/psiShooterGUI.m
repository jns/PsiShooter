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

% Last Modified by GUIDE v2.5 21-May-2009 20:49:58

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
global data psi energies vPath binPath;
%clear the global variables in case this program is being run twice in the
%same instance of matlab.
data = [];
psi = [];
vPath = [];
binPath = [];
energies = [];

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
        getPotential();
        set(handles.loadInput, 'String', 'SELECT YOUR NEW FILE');
    case 3
        set(handles.loadInput, 'String', 'Ex:"sin(2*pi/1e-6*x) 0:1e-9:1e-6"');
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
global data vPath; %access global potential variable

loadMenu_index = get(handles.loadMenu, 'Value');
currSysMessText = get(handles.systemMessages, 'String');
switch loadMenu_index
    case 1
        filePath = get(handles.loadInput,'String');
        
        %get the file path if none is specified in the text box
        if or(isempty(filePath),strcmp(filePath,'Input File Path'))
            [name,path] = uigetfile({'*.*'});
            if isequal(name,0) || isequal(path,0)
                return;
            end
            filePath = [path name];
            set(handles.loadInput,'String',name);
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
        vPath = filePath;
        [data,messages] = loadData(filePath);
        potential = data(1);
        
        if ~isempty(messages)
            currSysMessText =[messages;currSysMessText];
            set(handles.systemMessages, 'String', currSysMessText);
        end
        
        cla;%clear the plot window
        if isempty(potential.x)
            return
        elseif or(length(potential.x) == 1,length(potential.y)==1)
            if length(potential.x) == 1
                visualize1D(potential.y,potential.data,'red',0);
            else
                visualize1D(potential.x,potential.data,'red',0);
            end
        elseif and(length(potential.x)>1,length(potential.y)>1)
            visualize2D(potential.x,potential.y,potential.data,0)
        end
        %Enable the plot control toggle buttons
        set(handles.rotateEnable,'Enable','on');
        set(handles.panEnable,'Enable','on');
        set(handles.zoomEnable,'Enable','on');
        
    case 2
        set(handles.loadInput,'String','LOAD THE FILE YOU CREATED');
    case 3
        input = get(handles.loadInput,'String');
        %Example: "sin(x) 0:1e-9:1e-6"
        seperationIndex = find(input==' ');
        x = str2num(input(seperationIndex:end));
        potentialFunction = input(1:seperationIndex);
        %Scale into cm
        potential.x = x*1e2;
        potential.y = 0;
        potential.data = eval(potentialFunction);
        data = potential(1);
        data.data = data.data*1.60217646e-12;
        
        currSysMessText =[{['Load Function: ' potentialFunction ...
            ' from x=' num2str(min(potential.x)) ' to ' ...
            num2str(max(potential.x))]};currSysMessText];
        set(handles.systemMessages, 'String', currSysMessText);
        
        cla;%clear the plot window
        if isempty(potential.x)
            return
        elseif or(length(potential.x) == 1,length(potential.y)==1)
            if length(potential.x) == 1
                visualize1D(potential.y,potential.data,'red',0);
            else
                visualize1D(potential.x,potential.data,'red',0);
            end
        elseif and(length(potential.x)>1,length(potential.y)>1)
            visualize2D(potential.x,potential.y,potential.data,0)
        end
        %Enable the plot control toggle buttons
        set(handles.rotateEnable,'Enable','on');
        set(handles.panEnable,'Enable','on');
        set(handles.zoomEnable,'Enable','on');
        
        vPath = 'potentialFromFunction';
        
        %write the function defined potential to a file scaled into cm and ergs.
        writeFile(potential.x,potential.y,1.60217646e-12*potential.data,vPath);
    case 4
        currSysMessText =[{'Load Default Potential'};currSysMessText];
        set(handles.systemMessages, 'String', currSysMessText);
        data.x = [1,2,3,4,5,6,7,8,9];
        data.y = [1,2,3,4,5,6,7,8,9];
        data.data = [5,5,5,5,5,5,5,5,5;5,4,4,4,4,4,4,4,5;...
            5,4,3,2,2,2,3,4,5;5,4,2,2,1,2,2,4,5;5,4,2,1,1,1,2,4,5;...
            5,4,2,2,1,2,2,4,5;5,4,3,2,2,2,3,4,5;5,4,4,4,4,4,4,4,5;...
            5,5,5,5,5,5,5,5,5];
        visualize2D(data.x,data.y,data.data,0)
        %Enable the plot control toggle buttons
        set(handles.rotateEnable,'Enable','on');
        set(handles.panEnable,'Enable','on');
        set(handles.zoomEnable,'Enable','on');
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

%Precede the shell command with ! to execute.
global data psi energies vPath binPath;

simulateMenu_index = get(handles.simulateMenu, 'Value');
currSysMessText = get(handles.systemMessages, 'String');

if isempty(vPath)
    msgbox('No Potential Loaded. Please load the potential file first.');
    return;
end

unixOS = 1;

%get(handles.)

if unixOS
    [status,result] = unix('ls psi_shooter');
    if length(result) == 12
        %Then it found the binary
        binPath = './psi_shooter';
    else
        [status,result] = unix('ls ../src/psi_shooter');
        if length(result) == 19
            %then it found the binary in the source directory.
            binPath = '../src/psi_shooter ';
        else
            msgbox('Give me the path to the Psi Shooter binary!')
            [name,path] = uigetfile('','Select the Psi Shooter Binary');
            if isequal(filename,0) || isequal(pathname,0)
                return;
            end
            binPath = [path name];
        end
        [status,messages] = unix([binPath ' ' vPath]);
    end
    if ~isempty(messages)
        if length(messages) > 1000
            currSysMessText = [{'---PSISHOOTER BINARY MESSAGES---'};...
                {'Messages too long. Contents clipped.'};...
                {messages(end-256:end)}; ...
            {'---PSISHOOTER BINARY MESSAGES---'};...
            currSysMessText];
        else
        currSysMessText = [{'---PSISHOOTER BINARY MESSAGES---'};...
            {messages}; ...
            {'---PSISHOOTER BINARY MESSAGES---'};...
            currSysMessText];
        end
    end
else
    msgbox(['Your operating system is not yet fully supported. ' ...
        'Please edit psiShooterGUI.m to add support']);
end
set(handles.systemMessages, 'String',currSysMessText);

%messages(end-42:end-40) should be the place where the binary the number of
%solutions it found. If it found 
if str2num(messages(end-43:end-40)) == 0
    currSysMessText = [{'NO SOLUTIONS FOUND'};currSysMessText];
    set(handles.systemMessages, 'String',currSysMessText);
    return
end
%Assuming that these are in the present working directory
[psi,messages] = loadData('BS.dat');
currSysMessText =[messages;currSysMessText];

[energies,messages] = loadEnergies('E.txt');
currSysMessText =[messages;currSysMessText];
if isempty(energies)
    energies(1:length(psi)) = [1:length(psi)];
end
if length(energies) ~= length(psi)
    currSysMessText =[{'Different number of bound states'}; ...
        {'found than solutions found in file'};currSysMessText];
end
set(handles.systemMessages, 'String', currSysMessText);

try
    cla;%clear the plot window
    if isempty(psi(1).x)
        return
    elseif or(length(psi(1).x)==1,length(psi(1).y)==1)
        if length(psi(1).x) == 1 
            %checking to see if the length of the potential data is larget
            %than the length of the solution data is a cheap hack to see if
            %the data structure  is 2d or not. There are better ways of
            %doing this, but it is late, and I'd rather write a long
            %comment than correct what I just wrote.
            if isempty(data) || length(data(1).data) > length(psi(1).data)
                solInput = [{psi(1).y},...
                    {zeros(1,length(psi(1).y))},{'black'},{0}];
            else
                solInput = [{data(1).y},...
                    {data(1).data/1.60217646e-12},{'black'},{0}];
            end
            for index = 1:length(psi)
                if index == 1%Set the colors for each fucntion.
                    %This should be done with a function, but here's a
                    %quick soltion.
                    color = 'red';
                elseif index == 2
                    color = 'green';
                elseif index == 3
                    color = 'blue';
                else
                    color = 'randomize';
                end
                solInput = [solInput,{psi(index).y}, ...
                    {psi(index).data},{color},{energies(index)}];
            end
        else
            if isempty(data) || length(data(1).data) > length(psi(1).data)
                solInput = [{psi(1).x},...
                    {zeros(1,length(psi(1).x))},{'black'},{0}];
            else
                solInput = [{data(1).x},...
                    {data(1).data/1.60217646e-12},{'black'},{0}];
            end
            for index = 1:length(psi)
                if index == 1%Set the colors for each fucntion.
                    %This should be done with a function, but here's a
                    %quick soltion.
                    color = 'red';
                elseif index == 2
                    color = 'green';
                elseif index == 3
                    color = 'blue';
                else
                    color = 'randomize';
                end
                solInput = [solInput,{psi(index).x}, ...
                    {psi(index).data},{color},{energies(index)}];
            end
        end
        visualize1D(solInput);
        set(handles.rotateEnable,'Enable','off');
    elseif and(length(psi.x)>1,length(psi.y)>1)
        if isempty(data)
            solInput = [{psi(1).x},{psi(1).y},...
                {zeros(length(psi(1).x),length(psi(1).y))},{0}];
        else
            solInput = [{data(1).x},{data(1).y},...
                {data(1).data},{0}];
        end
        for index = 1:length(psi)
            solInput = [solInput,{psi(index).x},{psi(index).y}, ...
                {psi(index).data},{energies(index)}];
        end
        visualize2D(solInput)
        set(handles.rotateEnable,'Enable','on');
    end
    %Enable the plot control toggle buttons
    set(handles.panEnable,'Enable','on');
    set(handles.zoomEnable,'Enable','on');
catch
    
    currSysMessText =['Solutions failed to plot';currSysMessText];

    set(handles.systemMessages,'String',currSysMessText);
end


%%


% --- Executes on selection change in simulateMenu.
function simulateMenu_Callback(hObject, eventdata, handles)
% hObject    handle to simulateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns simulateMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        simulateMenu

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

'Arbitrary Breakpoint';


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


function filterButton_Callback(hObject, eventdata, handles)
% hObject    handle to filterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterButton as text
%        str2double(get(hObject,'String')) returns contents of filterButton as a double

global psi data energies;
filterString = get(handles.filterInput,'String');
parseCommas = find(filterString == ',');
parseDashes = find(filterString == '-');
parseSpaces = find(filterString == ' ');
if ~isempty(parseCommas)
    plotIndices = str2num(filterString(:,parseCommas(1)-1));
    for parseIndex = 1:length(parseCommas)
        if length(parseCommas) > parseIndex
            plotIndices(parseIndex+1) = str2num(filterString(...
                parseCommas(parseIndex):parseCommas(parseIndex+1)));
        else
            plotIndices(parseIndex+1) = str2num(filterString(...
                parseCommas(parseIndex):end));
        end
    end
elseif ~isempty(parseDashes)
    startValue = str2num(filterString(1:parseDashes(1)-1));
    endValue = str2num(filterString(parseDashes(1)+1:end));
    plotIndices = startValue:endValue;
elseif ~isempty(parseSpaces)
    plotIndices = str2num(filterString(:,parseSpaces(1)-1));
    for parseIndex = 1:length(parseSpaces)
        if length(parseSpaces) > parseIndex
            plotIndices(parseIndex+1) = str2num(filterString(...
                parseSpaces(parseIndex):parseSpaces(parseIndex+1)));
        else
            plotIndices(parseIndex+1) = str2num(filterString(...
                parseSpaces(parseIndex):end));
        end
    end
else
    plotIndices = 1:length(psi);
end

cla;%clear the plot window
if isempty(psi(1).x)
    return
elseif or(length(psi(1).x)==1,length(psi(1).y)==1)
    if length(psi(1).x) == 1
        %checking to see if the length of the potential data is larget
        %than the length of the solution data is a cheap hack to see if
        %the data structure  is 2d or not. There are better ways of
        %doing this, but it is late, and I'd rather write a long
        %comment than correct what I just wrote.
        if isempty(data) || length(data(1).data) > length(psi(1).data)
            solInput = [{psi(1).y},...
                {zeros(1,length(psi(1).y))},{'black'},{0}];
        else
            solInput = [{data(1).y},...
                {length(data(1).data)},{'black'},{0}];
        end

        for index = plotIndices
            if index == 1%Set the colors for each fucntion.
                %This should be done with a function, but here's a
                %quick soltion.
                color = 'red';
            elseif index == 2
                color = 'green';
            elseif index == 3
                color = 'blue';
            else
                color = 'randomize';
            end
            solInput = [solInput,{psi(index).y}, ...
                {psi(index).data},{color},{energies(index)}];
        end
    else
        if isempty(data) || length(data(1).data) > length(psi(1).data)
            solInput = [{psi(1).x},...
                {zeros(1,length(psi(1).x))},{'black'},{0}];
        else
            solInput = [{data(1).x},...
                {data(1).data/1.60217646e-12},{'black'},{0}];
        end
        for index = plotIndices
            if index == 1%Set the colors for each fucntion.
                %This should be done with a function, but here's a
                %quick soltion.
                color = 'red';
            elseif index == 2
                color = 'green';
            elseif index == 3
                color = 'blue';
            else
                color = 'randomize';
            end
            solInput = [solInput,{psi(index).x}, ...
                {psi(index).data},{color},{energies(index)}];
        end
    end
    visualize1D(solInput);
    set(handles.rotateEnable,'Enable','off');
elseif and(length(psi.x)>1,length(psi.y)>1)
    if isempty(data)
        solInput = [{psi(1).x},{psi(1).y},...
            {zeros(length(psi(1).x),length(psi(1).y))},{0}];
    else
        solInput = [{data(1).x},{data(1).y},...
            {data(1).data},{0}];
    end
    for index = plotIndices
        solInput = [solInput,{psi(index).x},{psi(index).y}, ...
            {psi(index).data},{energies(index)}];
    end
    visualize2D(solInput)
    set(handles.rotateEnable,'Enable','on');
end
%Enable the plot control toggle buttons
set(handles.panEnable,'Enable','on');
set(handles.zoomEnable,'Enable','on');

% --- Executes during object creation, after setting all properties.
function filterButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function filterInput_Callback(hObject, eventdata, handles)
% hObject    handle to filterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterInput as text
%        str2double(get(hObject,'String')) returns contents of filterInput as a double


% --- Executes during object creation, after setting all properties.
function filterInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%

function loadSolutionButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadSolutionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadSolutionButton as text
%        str2double(get(hObject,'String')) returns contents of loadSolutionButton as a double

global psi data energies;

currSysMessText = get(handles.systemMessages, 'String');

if ~iscell(currSysMessText)
    temp = currSysMessText;
    clear('currSysMessText');
    currSysMessText{1} = temp;
end

[name,path] = uigetfile({'*.*'});
filePath = [path name];
currSysMessText =[{['Loading Solutions in ' name]};currSysMessText];

[psi,messages] = loadData(filePath);
currSysMessText =[messages;currSysMessText];

[energies,messages] = loadEnergies([path 'E.txt']);
currSysMessText =[messages;currSysMessText];
if isempty(energies)
    energies(1:length(psi)) = [1:length(psi)];
end
if length(energies) ~= length(psi)
    currSysMessText =[{'Different number of bound states'}; ...
        {'found than solutions found in file'};currSysMessText];
end

set(handles.systemMessages, 'String', currSysMessText);

%Put    EVERYTHING below in the try-catch statements into a seperate
%function so that it can be called from the simulate  button too.
try
    cla;%clear the plot window
    if isempty(psi(1).x)
        return
    elseif or(length(psi(1).x)==1,length(psi(1).y)==1)
        if length(psi(1).x) == 1 
            %checking to see if the length of the potential data is larget
            %than the length of the solution data is a cheap hack to see if
            %the data structure  is 2d or not. There are better ways of
            %doing this, but it is late, and I'd rather write a long
            %comment than correct what I just wrote.
            if isempty(data) || length(data(1).data) > length(psi(1).data)
                solInput = [{psi(1).y},...
                    {zeros(1,length(psi(1).y))},{'black'},{0}];
            else
                solInput = [{data(1).y},...
                    {length(data(1).data)},{'black'},{0}];
            end
            for index = 1:length(psi)
                if index == 1%Set the colors for each fucntion.
                    %This should be done with a function, but here's a
                    %quick soltion.
                    color = 'red';
                elseif index == 2
                    color = 'green';
                elseif index == 3
                    color = 'blue';
                else
                    color = 'randomize';
                end
                solInput = [solInput,{psi(index).y}, ...
                    {psi(index).data},{color},{energies(index)}];
            end
        else
            if isempty(data) || length(data(1).data) > length(psi(1).data)
                solInput = [{psi(1).x},...
                    {zeros(1,length(psi(1).x))},{'black'},{0}];
            else
                solInput = [{data(1).x},...
                    {data(1).data/1.60217646e-12},{'black'},{0}];
            end
            for index = 1:length(psi)
                if index == 1%Set the colors for each fucntion.
                    %This should be done with a function, but here's a
                    %quick soltion.
                    color = 'red';
                elseif index == 2
                    color = 'green';
                elseif index == 3
                    color = 'blue';
                else
                    color = 'randomize';
                end
                solInput = [solInput,{psi(index).x}, ...
                    {psi(index).data},{color},{energies(index)}];
            end
        end
        visualize1D(solInput);
        set(handles.rotateEnable,'Enable','off');
    elseif and(length(psi.x)>1,length(psi.y)>1)
        if isempty(data)
            solInput = [{psi(1).x},{psi(1).y},...
                {zeros(length(psi(1).x),length(psi(1).y))},{0}];
        else
            solInput = [{data(1).x},{data(1).y},...
                {data(1).data},{0}];
        end
        for index = 1:length(psi)
            solInput = [solInput,{psi(index).x},{psi(index).y}, ...
                {psi(index).data},{energies(index)}];
        end
        visualize2D(solInput)
        set(handles.rotateEnable,'Enable','on');
    end
    %Enable the plot control toggle buttons
    set(handles.panEnable,'Enable','on');
    set(handles.zoomEnable,'Enable','on');
catch
    
    currSysMessText =['Solutions failed to plot';currSysMessText];

    set(handles.systemMessages,'String',currSysMessText);
end



% --- Executes during object creation, after setting all properties.
function loadSolutionButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadSolutionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%


function zoomEnable_Callback(hObject, eventdata, handles)
% hObject    handle to zoomEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoomEnable as text
%        str2double(get(hObject,'String')) returns contents of zoomEnable as a double
zoomEnableStatus = get(handles.zoomEnable,'Value');
if ~zoomEnableStatus
    zoom on;
    %set(handles.rotateEnable,'Value',1);
else
    zoom off;
    %set(handles.rotateEnable,'Value',0);
end

% --- Executes during object creation, after setting all properties.
function zoomEnable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%

function panEnable_Callback(hObject, eventdata, handles)
% hObject    handle to panEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of panEnable as text
%        str2double(get(hObject,'String')) returns contents of panEnable as a double
panEnableStatus = get(handles.panEnable,'Value');
if ~panEnableStatus
    pan on;
    %set(handles.rotateEnable,'Value',1);
else
    pan off;
    %set(handles.rotateEnable,'Value',0);
end

% --- Executes during object creation, after setting all properties.
function panEnable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%

function rotateEnable_Callback(hObject, eventdata, handles)
% hObject    handle to rotateEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotateEnable as text
%        str2double(get(hObject,'String')) returns contents of rotateEnable as a double
rotateEnableStatus = get(handles.rotateEnable,'Value');
if ~rotateEnableStatus
    rotate3d on;
    %set(handles.rotateEnable,'Value',1);
else
    rotate3d off;
    %set(handles.rotateEnable,'Value',0);
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function rotateEnable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotateEnable (see GCBO)
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
%ex: visualize1D([1 2 3 4 5],[1 2 3 2 1],'red',offset);
%input data structure:
% varargin{1+4*n} = X-axis array
% barargin{2+4*n} = Y-axis array
% style{3+4*n} = display style
% offset{4+4*n} = plot offset from 0.
%where n is the number of items being plotted
%style will pass matlab plot options directly to the plot command.

if nargin == 1 && length(varargin{1}) > 1 
    %fix the input issue. should probably put this back into properly
    %matlab data structure format at some later point, but it is just too
    %much of a pain at the moment.
    vararginC = varargin{1};
    narginC = length(vararginC);
    
else
    narginC = nargin;
    vararginC = varargin;
end

if mod(narginC,4) ~= 0 && narginC ~= 0;
    fprintf(['Bad input arguments. There should be a multiple \n' ...
        'of four input arguments. There are ' num2str(narginC) '.\n']);
    return;
end

for n = 1:4:narginC
    X(:,ceil(n/4)) = vararginC{n};
    Y(:,ceil(n/4)) = vararginC{n+1};
    if n > 1
        %Normalize the 1D plot values the Y values
        if ~((max(Y(:,1))-min(Y(:,1))) == 0)
            Y(:,ceil(n/4)) = (max(Y(:,1))-min(Y(:,1))) * (2e-7) * Y(:,ceil(n/4)).^2 /...
                (sum(Y(:,ceil(n/4)).^2) *(X(2)-X(1)));
        else
            Y(:,ceil(n/4)) = (2e-7) * Y(:,ceil(n/4)).^2 /...
                (sum(Y(:,ceil(n/4)).^2) *(X(2)-X(1)));
        end
    end
    color{ceil(n/4)} = vararginC{n+2};
    offset{ceil(n/4)} = vararginC{n+3};
end

legendText = [];
for n = 1:1:narginC/4
    %X axis multiplied by 1e7 to convert from cm to nm in display.
    if strcmp(color{n},'randomize')
        plot(X(:,n)*1e7,Y(:,n)+offset{n},'Color',[rand(1),rand(1),rand(1)]);
        legendText = [legendText;{['#' num2str(n) ' ' ...
            num2str(offset{n}) ' eV']}];
    elseif strcmp(color{n},'black')
        plot(X(:,n)*1e7,Y(:,n)+offset{n},color{n})
    else
        plot(X(:,n)*1e7,Y(:,n)+offset{n},color{n})
        legendText = [legendText;{['#' num2str(n) ' ' ...
            num2str(offset{n}) ' eV']}];
    end
    if n == 1
        legendText = {'Potential'};
    end
    hold on;
end
legend(legendText);
xlabel('Position (nm)');
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

if nargin == 1 && length(varargin{1}) > 1 
    %fix the input issue. should probably put this back into properly
    %matlab data structure format at some later point, but it is just too
    %much of a pain at the moment.
    vararginC = varargin{1};
    narginC = length(vararginC);
    
else
    narginC = nargin;
    vararginC = varargin;
end

if mod(narginC,4) ~= 0 && narginC ~= 0;
    fprintf(['Bad input arguments. There should be a multiple \n' ...
        'of four input arguments. There are ' num2str(narginC) '.\n']);
    return;
end

%plot the wavefunctions and potential together.

%axes(handles.simPlot);%Can't select this one of multiple inputs unless we
%pass the handle structure into this function.

for n = 1:4:narginC
    if n == 1 %plot styles set up for the potential (more transparent)
        h(1)=surf(vararginC{n}*1e7,vararginC{n+1}*1e7,vararginC{n+2}/1.60217646e-12);
        set(h(1),'facealpha',0.35);
        set(h(1),'edgealpha',0.01);
        hold on;
    else %plot styles set up for wavefunctions
        %h(n) is the handle for the graphic objects I am creating. They can
        %be accessed by get or set(h(n)) with whatever property you want to
        %mess with.
        h(n)=surf(vararginC{n}*1e7,vararginC{n+1}*1e7,vararginC{n+2}+vararginC{n+3});
        set(h(n),'facealpha',0.7);
        set(h(n),'edgealpha',0.1);
    end
end
hold off;
'end';% a good place to stick a breakpoint for playing with the plots.
%make sure to turn 'rotate3d on;'

%%

function [energies,messages] = loadEnergies(path)
%function energies = loadEnergies(path)
try
    fidP = fopen(path,'r','ieee-le.l64');
    header = fgetl(fidP);
    energies = fscanf(fidP,'%f');
    messages = {'Solution loaded Successfully'};
catch
    energies = [];
    messages = {'Energy File Failed to load'};
    return;
end

fclose(fidP);

%%

function [data,messages] = loadData(path)
%function data = loadData(path)
%data.x = [x1,x2,x3,x4,...]
%data.y = [y1,y2,y3,y4,...]
%data.data = [d_11,d_12,d_13,d_14;d_21,d_22,d_23,d24;...]
%Works with binary or delimited ascii/unicode data. Save the data file as a
%.txt file if you want to read it in as ascii or unicode. Otherwise, it
%will default to binary.

try
    %fidP = fopen(path,'r','ieee-le.l64');
    fidP = fopen(path,'r');
catch
    data = [];
    messages = {'DATA FILE FAILED TO LOAD'};
    return
end
if fidP == -1
    data = [];
    messages = {'BAD FILE PATH'};
    return
end

try
    index = 0;
    fseek(fidP,0,'eof');
    fileEnd = ftell(fidP);
    fseek(fidP,0,'bof');
    while ftell(fidP) < fileEnd - 8
        index = index + 1;
        messages = [];
        if strcmp(path(end-3:end),'.txt')
            %ascii or unicode
            xNum(index) = fscanf(fidP,'%e',1);
            yNum(index) = fscanf(fidP,'%e',1);
            data(index).x = fscanf(fidP,'%e',xNum(index));
            data(index).y = fscanf(fidP,'%e',yNum(index));
            if yNum(index) == 0
                data(index).data = fscanf(fidP,'%e',[xNum(index),1]);
                if size(data.data,1)*size(data.data,2) ~= xNum(index)*yNum(index)
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'};...
                        {'It might not like how your decimal numbers are formatted'}];
                end
            else
                data(index).data = fscanf(fidP,'%e',[xNum(index),yNum(index)]);
                if size(data.data,1)*size(data.data,2) ~= xNum(index)*yNum(index)
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'};...
                        {'It might not like how your decimal numbers are formatted'}];
                end
            end
            if isempty(xNum)
                messages ={'DATA MISLOADED. xNum is empty...'};
            end
        else
            %binary
            xNum(index) = fread(fidP,1,'float64');
            if(xNum(1) < 1) %If it looks like the file is wrong endian,
                %close the file and re-open it with little endian.
                fclose(fidP);
                fidP = fopen(path,'r','ieee-be.l64');
                xNum(index) = fread(fidP,1,'float64');
            end
            yNum(index) = fread(fidP,1,'float64');
            data(index).x = fread(fidP,xNum(index),'float64');
            data(index).y = fread(fidP,yNum(index),'float64');
            if or(yNum(index) == 0, yNum(index) == 1)
                data(index).data = fread(fidP,[xNum(index),1],'float64');
                if size(data(index).data,1)*1 ~= xNum(index)*1
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'}];
                end
            else
                data(index).data = fread(fidP,[xNum(index),yNum(index)],'float64');
                if size(data(index).data,1)*size(data(index).data,2) ~= xNum(index)*yNum(index)
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'}];
                end
            end
            if isempty(xNum)
                messages =[{'DATA MISLOADED'};...
                    {'Is your file opening a big endian (PowerPC) file'};...
                    {'file on a little endian (x86) machine?'}];
            end
        end
    end
    
    fclose(fidP);
catch
    data = [];
    messages = {'DATA FAILED TO LOAD'};
    fclose(fidP);
    return
end

%%
%grabbed from getPotential.m
function writeFile(X,Y,Data,path)
%writeFile(X,Y,Data,fileName)

if Y == 0
    yLength = 1;
else
    yLength=length(Y);
end

try
    %fidP = fopen(path,'w','ieee-le.l64');
    fidP = fopen(path,'w');
    if strcmp(path(end-3:end),'.txt')
        %ascii or unicode
        fprintf(fidP,'%e\n',length(X));
        fprintf(fidP,'%e\n',yLength);
        fprintf(fidP,'%e ',X);
        %fprintf(fidP,'%s','\n');
        fprintf(fidP,'%e ',Y);
        %fprintf(fidP,'%s','\n');
        fprintf(fidP,'%e ',Data);
    else
        %binary
        fwrite(fidP,length(X),'float64');
        fwrite(fidP,yLength,'float64');
        fwrite(fidP,X,'float64');
        fwrite(fidP,Y,'float64');
        fwrite(fidP,Data,'float64');
    end
    
    fclose(fidP);
catch
    fclose(fidP);
    return
end


% --------------------------------------------------------------------
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuExit_Callback(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fclose('all');
clear('all');
close('all');

% --------------------------------------------------------------------
function menuCredits_Callback(hObject, eventdata, handles)
% hObject    handle to menuCredits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox([{'PSI SHOOTER TEAM'};{''};{'Jere Harrison'};{'Joshua Shapiro'};{'Cyrus Haselby'}])
