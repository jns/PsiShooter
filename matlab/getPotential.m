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
%clear('layerData');
layerData = [];
layerData(1).name = 'Free Space';
layerData(1).V = 0;
layerData(1).thickness = '0';%NEED TO CONVERT THIS STRING TO A NUMBER LATER!
layerData(1).epsRel = 1;
layerData(1).mEff = 1;
%this table is full of values from all over the internet. 
%I would suggest that you update it with experimental data or more
%authroritative date before you put too much faith in it.

%The second values of things are either documented or the correction based
%on the alloy percentage.
layerList(1).name = 'C (diamond)';  layerList(1).eg=5.48;
layerList(1).mEff(1) = 1.4;         layerList(1).mEff(2) = 0.36; %longitudinal and transferse
layerList(1).epsRel = 5.7;
layerList(2).name = 'Si';           layerList(2).eg = 1.11;
layerList(2).epsRel = 11.9;         layerList(2).mEff = 0.2;
layerList(3).name = 'Ge';           layerList(3).eg = 0.67;
layerList(3).epsRel = 16;           layerList(3).mEff = 0.041;
layerList(4).name = 'AlN';          layerList(4).eg = 6.3;
layerList(4).epsRel = 8.5;          layerList(4).mEff = [];
layerList(5).name = '*AlP';          layerList(5).eg = 2.45;
layerList(5).epsRel = [];           layerList(5).mEff = [];
layerList(6).name = '*AlAs';         layerList(6).eg = 2.16;
layerList(6).epsRel = [];           layerList(6).mEff = [];
layerList(7).name = '*GaN';          layerList(7).eg = 3.4;
layerList(7).epsRel = [];           layerList(7).mEff = [];
layerList(8).name = '*GaP';          layerList(8).eg = 2.26;
layerList(8).epsRel = [];           layerList(8).mEff = [];
layerList(9).name = 'GaAs';         layerList(9).eg = 1.424;
layerList(9).epsRel = 12.5;         layerList(9).mEff = 0.067;
layerList(10).name = 'AlGaAs';      layerList(10).eg = 2.16;
layerList(10).eg(2)=0.1*(0.36-2.16);layerList(10).epsRel(1) = 12.9;
layerList(10).epsRel(2) = 2.84;     layerList(10).mEff(1) = 0.063;
layerList(10).mEff(2) = 0.083;
layerList(11).name ='SiO2';         layerList(11).eg = 9.00;
layerList(11).epsRel = 3.9;          layerList(11).mEff = [];
layerList(12).name ='Si3N4';        layerList(12).eg = 5.30;
layerList(12).epsRel = 7.5;          layerList(12).mEff = [];
layerList(13).name ='H';            layerList(13).eg = 0;
layerList(13).epsRel = 1;           layerList(13).mEff = 1;
layerList(14).name = 'Al_15_Ga_85_As';      layerList(10).eg = 1.611;
layerList(14).eg(2)=0.1*(0.36-2.16);layerList(10).epsRel(1) = 12.9;
layerList(14).epsRel(2) = 2.84;     layerList(10).mEff(1) = 0.063;
layerList(14).mEff(2) = 0.083;

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
layerData(index+1).epsRel = layerData(index).epsRel;
layerData(index+1).mEff = layerData(index).mEff;
layerData(index).name = layerList(layerMenu_index).name;
layerData(index).V = layerList(layerMenu_index).eg(1);%PULLING POTENTIAL STRAIGHT FROM BANDGAP.
layerData(index).thickness = get(handles.layerThicknessInput,'String');
layerData(index).epsRel = layerList(layerMenu_index).epsRel(1);
layerData(index).mEff = layerList(layerMenu_index).mEff(1);

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
bias = str2num(get(handles.biasInput,'String'));
geometrySelect_index = get(handles.symmetrySelect,'Value');

global layerData;
fileName = get(handles.fileNameInput,'String');
layerDataTemp = layerData;%Preserve the global data in string form
for index = 1:length(layerData)
    layerDataTemp(index).thickness = str2num(layerData(index).thickness);
end

meshGen = constructPotential(layerDataTemp,meshResolution, ...
    geometrySelect_index,bias);

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
    geometrySelect_index,bias)
%meshGen = constructPotential(layerData,meshResolution, ...
%geometrySelect_index,electrical_bias)
%Generates the mesh. Very inefficiently written, but gets the job done
%implemented for linear potentials and circularly symmetric potentials now.

%this should all be in angstroms. Converted to meters at the end.
geometrySizeExtra = 0;
geometrySize = geometrySizeExtra;
for index = 1:length(layerData)
    geometrySize = geometrySize + layerData(index).thickness;
end

if geometrySelect_index == 6
    meshGen.X = 0:meshResolution:geometrySize;
    meshGen.Y = [];
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

%>> div D = charge density = 0, so D is constant in regions of 0 charge
%density (ie, the dielectric we are modeling. I know that it won't be
%true for non charge-neutral areas, but it should be OK as long as we
%dont worry about build in potentials formed during device operation.
%
%>> Vbias = int(E dl)
%>> Vbias = int(D / eps1 dl)from 0 to d1+int(D / eps2 dl)from d1 to d2+...
%>> Vbias = D*d1/eps1+D*(d2-d1)/eps2+... = D*(d1/eps1+(d2-d1/eps2)+...)
%>> D = Vbias / (d1/eps1+(d2-d1/eps2)+...)
%>> E = D/Eps, so E1 = D/Eps1, E2 = D / Eps2, etc.
thOverEpsTot = 0;
for layerIndex = 1:length(layerData)
    thOverEpsTot = thOverEpsTot + layerData(layerIndex).thickness / ...
        layerData(layerIndex).epsRel(1);%Add in epsRel(2) to deal with alloys
end
D = bias / thOverEpsTot;
for layerIndex = 1:length(layerData)-1%Exclude the free space stuck in automatically. It will skew the last point.
    layerData(layerIndex).E = D / layerData(layerIndex).epsRel(1);
end
layerData(end).E = 0;
xOffset(1) = 0;%for calculation of meshed bias
for  index = 2:length(layerData)
    xOffset(index) = xOffset(index-1)+layerData(index-1).thickness;
end

if geometrySelect_index == 6%Just 1-D    
    %Look up the approprite layer data for the current mesh point.
    meshGen.Data = zeros(1,length(meshGen.X));
    meshBias = zeros(1,length(meshGen.X));
    for indX = 1:length(meshGen.X)
        for index = 1:length(layerData)
            if meshGen.X(indX) <= thickness(index)
                meshGen.Data(indX) = layerData(index).V(1);
                if indX > 1
                    meshBias(indX) = meshBias(indX-1)+meshResolution * layerData(index).E;
                else
                    meshBias(indX) = 0;
                end
                break;
            end
        end
    end
    meshGen.Data = meshGen.Data + meshBias;
    
elseif geometrySelect_index == 1 %circular symmetry
    meshGen.Data = zeros(length(meshGen.X),length(meshGen.Y));
    meshBias = (-1)*bias*ones(length(meshGen.X),length(meshGen.Y));
    %Go through each and every point and
    %figure out what the potential should be set to.
    for indX = 1:length(meshGen.X)
        for indY = 1:length(meshGen.Y)
            for index = 1:length(layerData)
                if sqrt(meshGen.X(indX)^2+meshGen.Y(indY)^2) ...
                        <= thickness(index)
                    if indX == 46 && indY == 94
                        'BREAKPOINT';
                    end
                    %This will keep changing it until it is the right one.
                    meshGen.Data(indX,indY) = layerData(index).V;
                    if index == 1
                        meshBias(indX,indY) = (-1)*(sqrt((meshGen.X(indX))^2 ...
                            +(meshGen.Y(indY))^2)-0)* ...
                            layerData(index).E;
                    else%if index < length(layerData) - 1 && index > 1
                        biasOffset = 0;
                        for biasIndex = 1:index-1
                            biasOffset = biasOffset+layerData(biasIndex).thickness* ...
                                layerData(biasIndex).E;
                        end
                        meshBias(indX,indY) = (-1)*(sqrt((meshGen.X(indX))^2 ...
                            +(meshGen.Y(indY))^2)-thickness(index-1))* ...
                            layerData(index).E - biasOffset;
                    end
                    break;
                end
            end
        end
    end
    meshBias = meshBias - min(min(meshBias));
    meshGen.Data = meshGen.Data+meshBias;
end


%Convert mesh from Angstroms to centimeters.
meshGen.X = meshGen.X*1e-8;
if ~isempty(meshGen.Y)
    meshGen.Y = meshGen.Y*1e-8;
end
%Convert mesh from eV to ergs
meshGen.Data = meshGen.Data*1.60217646e-12;
% if length(meshGen.X) == 0 || length(meshGen.Y) == 0
%     meshGen.Data = meshGen.Data - min(meshGen.Data);
% else
%     meshGen.Data = meshGen.Data - min(min(meshGen.Data));
% end

%%

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