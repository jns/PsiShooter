function varargout = getPotential(varargin)
% GETPOTENTIAL M-file for getPotential.fig
%      GETPOTENTIAL, by itself, creates a new GETPOTENTIAL or raises the existing
%      singleton*.
%
%      H = GETPOTENTIAL returns the handle to a new GETPOTENTIAL or the handle to
%      the existing singleton*.
%
%      GETPOTENTIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETPOTENTIAL.M with the given input arguments.
%
%      GETPOTENTIAL('Property','Value',...) creates a new GETPOTENTIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getPotential_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getPotential_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getPotential

% Last Modified by GUIDE v2.5 19-May-2009 23:31:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getPotential_OpeningFcn, ...
                   'gui_OutputFcn',  @getPotential_OutputFcn, ...
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

% --- Executes just before getPotential is made visible.
function getPotential_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getPotential (see VARARGIN)

%Global variable list
global layerData layerList;
%Need to clear the global variable in case it is left around in memory from
%the last time it was used 
clear('layerData');
layerData = [];
layerData(1).name = 'Free Space';
layerData(1).V = 0;
layerData(1).thickness = '0';%NEED TO CONVERT THIS STRING TO A NUMBER LATER!
%this table is full of values from all over the internet. 
%I would suggest that you update it with experimental data or more
%authroritative date before you put too much faith in it.
layerList(1).name = 'C (diamond)';  layerList(1).eg=5.48;
layerList(1).mEff(1) = 1.4;         layerList(1).mEff(2) = 0.36; %longitudinal and transferse
layerList(2).name = 'Si';           layerList(2).eg = 1.11;
layerList(3).name = 'Ge';           layerList(3).eg = 0.67;
layerList(4).name = 'AlN';          layerList(4).eg = 6.3;
layerList(5).name = 'AlP';          layerList(5).eg = 2.45;
layerList(6).name = 'AlAs';         layerList(6).eg = 2.16;
layerList(7).name = 'GaN';          layerList(7).eg = 3.4;
layerList(8).name = 'GaP';          layerList(8).eg = 2.26;
layerList(9).name = 'GaAs';         layerList(9).eg = 0.36;
layerList(10).name ='SiO2';         layerList(10).eg = 9.00;
layerList(11).name ='Si3N4';        layerList(11).eg = 5.30;
layerList(12).name ='H';            layerList(12).eg = 0;

%Until this table is properly filled out, taking the potentials for the
%array to be the bandgap.
layerlist(:).V = layerList(:).eg;

% Choose default command line output for getPotential
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes getPotential wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = getPotential_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%

% --- Executes on selection change in symmetrySelect.
function symmetrySelect_Callback(hObject, eventdata, handles)
% hObject    handle to symmetrySelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns symmetrySelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from symmetrySelect


% --- Executes during object creation, after setting all properties.
function symmetrySelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to symmetrySelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',{'Circular','Octagonal','Hexagonal','Square', ...
    'Triangular','1D'});

%%

% --- Executes on selection change in layerTypeSelect.
function layerTypeSelect_Callback(hObject, eventdata, handles)
% hObject    handle to layerTypeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns layerTypeSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from layerTypeSelect


% --- Executes during object creation, after setting all properties.
function layerTypeSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layerTypeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global layerList;
for index = 1:length(layerList)
    dropBoxItems{index} = layerList(index).name;
end

set(hObject,'String',dropBoxItems);

%%

function layerThicknessInput_Callback(hObject, eventdata, handles)
% hObject    handle to layerThicknessInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of layerThicknessInput as text
%        str2double(get(hObject,'String')) returns contents of layerThicknessInput as a double


% --- Executes during object creation, after setting all properties.
function layerThicknessInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layerThicknessInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%

% --- Executes on button press in addLayerButton.
function addLayerButton_Callback(hObject, eventdata, handles)
% hObject    handle to addLayerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global layerList layerData; %access global potential variable

layerMenu_index = get(handles.layerTypeSelect, 'Value');
index = length(layerData);

layerData(index+1).name = layerData(index).name;
layerData(index+1).V = layerData(index).V;
layerData(index+1).thickness = layerData(index).thickness;
layerData(index).name = layerList(layerMenu_index).name;
layerData(index).V = layerList(layerMenu_index).V;
layerData(index).thickness = get(handles.layerThicknessInput,'String');

for displayIndex = 1:length(layerData)
    layerDisplayText(displayIndex) = {[layerData(displayIndex).thickness ...
        'A  ' layerData(displayIndex).name]};
end
set(handles.layerDisplay, 'String',layerDisplayText);

%%

% --- Executes on button press in deleteLayerButton.
function deleteLayerButton_Callback(hObject, eventdata, handles)
% hObject    handle to deleteLayerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

layerMenu_index = get(handles.layerDisplay, 'Value');
global layerData;
if layerMenu_index ~= length(layerData)
    layerData(layerMenu_index) = [];
end

for displayIndex = 1:length(layerData)
    layerDisplayText(displayIndex) = {[layerData(displayIndex).thickness ...
        'A  ' layerData(displayIndex).name]};
end
set(handles.layerDisplay, 'String',layerDisplayText);

%%

% --- Executes on selection change in layerDisplay.
function layerDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to layerDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns layerDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from layerDisplay



% --- Executes during object creation, after setting all properties.
function layerDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layerDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%

% --- Executes on button press in genPotentialButton.
function genPotentialButton_Callback(hObject, eventdata, handles)
% hObject    handle to genPotentialButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Gather the required data
meshResolution = str2num(get(handles.meshResolutionInput,'String'));
geometrySelect_index = get(handles.symmetrySelect,'Value');
global layerData;
fileName = get(handles.fileNameInput,'String');
layerDataTemp = layerData;%Preserve the global data in string form
for index = 1:length(layerData)
    layerDataTemp(index).thickness = str2num(layerData(index).thickness);
end

meshGen = constructPotential(layerDataTemp,meshResolution, ...
    geometrySelect_index);

X = meshGen.X;
if geometrySelect_index ~= 6
    Y = meshGen.Y;
else
    Y = 0;
end
Data = meshGen.Data;
writeFile(X,Y,Data,fileName);

%%

function fileNameInput_Callback(hObject, eventdata, handles)
% hObject    handle to fileNameInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileNameInput as text
%        str2double(get(hObject,'String')) returns contents of fileNameInput as a double


% --- Executes during object creation, after setting all properties.
function fileNameInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileNameInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%

function meshResolutionInput_Callback(hObject, eventdata, handles)
% hObject    handle to meshResolutionInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meshResolutionInput as text
%        str2double(get(hObject,'String')) returns contents of meshResolutionInput as a double


% --- Executes during object creation, after setting all properties.
function meshResolutionInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meshResolutionInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%

function meshGen = constructPotential(layerData,meshResolution, ...
    geometrySelect_index)
%meshGen = constructPotential(layerData,meshResolution, ...
%geometrySelect_index)
%Generates the mesh. Very inefficiently written, but gets the job done
%implemented for linear potentials and circularly symmetric potentials now.

%this should all be in angstroms. Might need to change that later
geometrySizeExtra = 10;
geometrySize = geometrySizeExtra;
for index = 1:length(layerData)
    geometrySize = geometrySize + layerData(index).thickness;
end

if geometrySelect_index == 6
    meshGen.X = 0:meshResolution:geometrySize;
else
    meshGen.X = -geometrySize:meshResolution:geometrySize;
    meshGen.Y = -geometrySize:meshResolution:geometrySize;
end

%Yeah, I know this is an inefficient way to do this, I'm just throwing it
%together.
thickness = zeros(1,length(layerData));
thickness(1) = layerData(1).thickness;
for index = 2:length(layerData)
    thickness(index) = thickness(index-1) + layerData(index).thickness;
end

if geometrySelect_index == 6%Just 1-D
    meshGen.Data = zeros(1,length(meshGen.X));
    for indX = 1:length(meshGen.X)
        for index = 1:length(layerData)
            if indX*meshResolution < thickness(index)
                meshGen.Data(indX) = layerData(index).V;
                break;
            end
        end
    end
%     meshGen.Data = layerData(1).V;%load the potential function
%     for index = 1:length(layerData)
%         %rounds boundaries up to the nearest whole number index. Will
%         %overestimate the width of the boundaries.
%         meshGen.Data(end:end+ceil(meshResolution* ...
%             layerData(index).thickness)) = layerData(index).V;
%     end
elseif geometrySelect_index == 1 %circular symmetry
    meshGen.Data = zeros(length(meshGen.X),length(meshGen.Y));
    %Go through each and every point and
    %figure out what the potential should be set to.
    for indX = 1:length(meshGen.X)
        for indY = 1:length(meshGen.Y)
            for index = 1:length(layerData)
                if sqrt(meshGen.X(indX)^2+meshGen.Y(indY)^2) ...
                        < thickness(index)
                    %This will keep changing it until it is the right one.
                    meshGen.Data(indX,indY) = layerData(index).V;
                    break;
                end
            end
        end
    end
end

%%

function writeFile(X,Y,Data,path)
%writeFile(X,Y,Data,fileName)

if Y == 0
    yLength = 0;
else
    yLength=length(Y);
end

try
    fidP = fopen(path,'w','ieee-le.l64');
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