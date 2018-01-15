function varargout = MCC(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MCC_OpeningFcn, ...
                   'gui_OutputFcn',  @MCC_OutputFcn, ...
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


% --- Executes just before MCC is made visible.
function MCC_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for MCC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


clc;

pan(handles.a1,'on');

set(handles.ac1,'Value',1);
global selectedAC;
selectedAC = get(handles.uib1,'SelectedObject');


global ErrorFlag1 ErrorFlag2 ErrorFlag3;
ErrorFlag1 = true;
ErrorFlag2 = true;
ErrorFlag3 = true;


% Initialize AutoSaveFlag
global AutoSaveFlag;
AutoSaveFlag = get(handles.c1,'Value');

% Initialize LoadFrameFlag
global LoadFrameFlag;
LoadFrameFlag = false;

global output;
output = struct;

% Save   Output file every "saveOutputStructTimeFreq"   sec
% Backup Output file every "backupOutputStructTimeFreq" sec

global saveOutputStructTimeFreq backupOutputStructTimeFreq;
saveOutputStructTimeFreq = 30;      % Sec
backupOutputStructTimeFreq = 5*60;  % Sec

% --- Outputs from this function are returned to the command line.
function varargout = MCC_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on selection change in b1.
function b1_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b1);

% --- Executes during object creation, after setting all properties.
function b1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in b2.
function b2_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b2);
    

% --- Executes during object creation, after setting all properties.
function b2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in b3.
function b3_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b3);


% --- Executes during object creation, after setting all properties.
function b3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in b5.
function b5_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b5);


% --- Executes during object creation, after setting all properties.
function b5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in b6.
function b6_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b6);

% --- Executes during object creation, after setting all properties.
function b6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in b7.
function b7_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b7);

% --- Executes during object creation, after setting all properties.
function b7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in b8.
function b8_Callback(hObject, eventdata, handles)
CheckListSelection(handles.b8);

% --- Executes during object creation, after setting all properties.
function b8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- CheckListSelection
function CheckListSelection(b)
% Reset the list selection to number 3 (unannotated) if 1 and 2 are chosen
if get(b,'Value') == 1 || get(b,'Value') == 2
   set(b,'Value',3);
end

function e1_Callback(hObject, eventdata, handles)
str = get(handles.e1,'String');

[str, errorStatus] = CheckVariableName(handles.e1, str);
if errorStatus
    return;
end
set(handles.e1,'String',str);

data = get(handles.b1,'String');
AddStructFieldFun(str,handles.b1,handles.e1,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function e2_Callback(hObject, eventdata, handles)
str = get(handles.e2,'String');


[str, errorStatus] = CheckVariableName(handles.e2, str);
if errorStatus
    return;
end
set(handles.e2,'String',str);

data = get(handles.b2,'String');
AddStructFieldFun(str,handles.b2,handles.e2,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e3_Callback(hObject, eventdata, handles)
str = get(handles.e3,'String');

[str, errorStatus] = CheckVariableName(handles.e3, str);
if errorStatus
    return;
end
set(handles.e3,'String',str);

data = get(handles.b3,'String');
AddStructFieldFun(str,handles.b3,handles.e3,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function e5_Callback(hObject, eventdata, handles)
str = get(handles.e5,'String');

[str, errorStatus] = CheckVariableName(handles.e5, str);
if errorStatus
    return;
end
set(handles.e5,'String',str);

data = get(handles.b5,'String');
AddStructFieldFun(str,handles.b5,handles.e5,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e6_Callback(hObject, eventdata, handles)
str = get(handles.e6,'String');

[str, errorStatus] = CheckVariableName(handles.e6, str);
if errorStatus
    return;
end
set(handles.e6,'String',str);

data = get(handles.b6,'String');
AddStructFieldFun(str,handles.b6,handles.e6,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function e7_Callback(hObject, eventdata, handles)
str = get(handles.e7,'String');

[str, errorStatus] = CheckVariableName(handles.e7, str);
if errorStatus
    return;
end
set(handles.e7,'String',str);

data = get(handles.b7,'String');
AddStructFieldFun(str,handles.b7,handles.e7,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e8_Callback(hObject, eventdata, handles)
str = get(handles.e8,'String');

[str, errorStatus] = CheckVariableName(handles.e8, str);
if errorStatus
    return;
end
set(handles.e8,'String',str);

data = get(handles.b8,'String');
AddStructFieldFun(str,handles.b8,handles.e8,data{1},handles);

% --- Executes during object creation, after setting all properties.
function e8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ac1.
function ac1_Callback(hObject, eventdata, handles)
val = get(handles.ac1,'Value');
if val == 1 
    set(handles.ac1,'ForegroundColor',[0 0 1]);
    pan(handles.a1,'on');
    
    set(handles.ac2,'ForegroundColor',[1 0 0]);
    set(handles.ac2,'Value',0);
    set(handles.ac3,'ForegroundColor',[1 0 0]);
    set(handles.ac3,'Value',0);
else
    set(handles.ac1,'ForegroundColor',[1 0 0]);
    pan(handles.a1,'off');
end

global selectedAC;
selectedAC = get(handles.uib1,'SelectedObject');

% --- Executes on button press in ac2.
function ac2_Callback(hObject, eventdata, handles)
val = get(handles.ac2,'Value');
if val == 1 
    set(handles.ac2,'ForegroundColor',[0 0 1]);
    rotate3d(handles.a1,'on');
    
    set(handles.ac1,'ForegroundColor',[1 0 0]);
    set(handles.ac1,'Value',0);
    set(handles.ac3,'ForegroundColor',[1 0 0]);
    set(handles.ac3,'Value',0);
else
    set(handles.ac2,'ForegroundColor',[1 0 0]);
    rotate3d(handles.a1,'off');
end

global selectedAC;
selectedAC = get(handles.uib1,'SelectedObject');

% --- Executes on button press in ac3.
function ac3_Callback(hObject, eventdata, handles)
val = get(handles.ac3,'Value');
if val == 1 
    set(handles.ac3,'ForegroundColor',[0 0 1]);
    zoom(handles.a1,'on');
    
    set(handles.ac1,'ForegroundColor',[1 0 0]);
    set(handles.ac1,'Value',0);
    set(handles.ac2,'ForegroundColor',[1 0 0]);
    set(handles.ac2,'Value',0);
else
    set(handles.ac3,'ForegroundColor',[1 0 0]);
    zoom(handles.a1,'off');
end

global selectedAC;
selectedAC = get(handles.uib1,'SelectedObject');

% --- Executes on button press in ac4.
function ac4_Callback(hObject, eventdata, handles)
global selectedAC;

% Reset figure view and x,y,z range
set(handles.a1,'Color','k');
grid(handles.a1,'on');
set(handles.a1,'Xcolor',[1 1 1]);
set(handles.a1,'Ycolor',[1 1 1]);
set(handles.a1,'Zcolor',[1 1 1]);

zlim(handles.a1,[-150 200]);
xlim(handles.a1,[-400 400]);
ylim(handles.a1,[-400 400]);
set(handles.a1,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
view(handles.a1,[-48, 12]);

set(handles.(selectedAC.Tag),'Value',1);
set(handles.uib1,'SelectedObject',selectedAC);


% --- Executes on button press in p2.
function p2_Callback(hObject, eventdata, handles)
global ErrorFlag1 output outputFieldNames numCols;
[FileName,FileDir] = uigetfile('*.xlsx','Select the Excel file');

if isnumeric(FileName) || isnumeric(FileDir)
    return;
end

FilePath = strcat(FileDir,FileName);

[FileDir,FileNameNoExt,FileNameExt] = fileparts(FilePath);

if length(FileNameExt) == 5
    if FileNameExt == '.xlsx'
        set(handles.t1,'String',strcat('Excel File: ',FileName));
        set(handles.t1,'ForegroundColor',[0,0,0]);
        ErrorFlag1 = false;
    else
        ErrorFlag1 = true;
    end
else
    ErrorFlag1 = true;
end

if ErrorFlag1
    set(handles.t1,'String','Error: load *.xlsx file only');
    set(handles.t1,'ForegroundColor',[1,0,0]);
    return;
end

ReadExcelFlag = false;

numCols = 7;

if ~ErrorFlag1
    [num,txt] = xlsread(FilePath);
    
    SizTXT = size(txt);
    
    if SizTXT(2) ~=numCols
        set(handles.t1,'String','Error: wrong number of columns');
        set(handles.t1,'ForegroundColor',[1,0,0]);
        return;
    end
    
    numRows = SizTXT(1);
    ReadExcelFlag = true;
    
end   

if ReadExcelFlag
    
    ErrorFlag1 = false;
    
    handlesMAT = [handles.b2 handles.b1 handles.b3 handles.b5 handles.b6 handles.b7 handles.b8];
    
    for j = 1:numCols
        
        % Clear column j
        set(handlesMAT(j),'String','');
        % Read column j
        col = txt(:,j);
        for i = 1:numRows
            if isempty(col{i})
                break; % goto next column, goto XX
            else
                data = get(handlesMAT(j),'String');
                data{end+1} = col{i};
                set(handlesMAT(j),'String',data);
                if i == 1
                    data = get(handlesMAT(j),'String');
                    data{end+1} = '-------';
                    set(handlesMAT(j),'String',data);
                    
                    % Build output struct variable
                    output.(col{i}) = [];
                else
                    % Build output struct variable
                    output.(col{1}).(col{i}) = [];
                end
            end
        end
        
        % XX: break comes here
        
        set(handlesMAT(j),'Value',3);
        
    end
    outputFieldNames = fieldnames(output);
   
else
    ErrorFlag1 = true;
end


% --- Executes on button press in p3.
function p3_Callback(hObject, eventdata, handles)
global mocapstruct ErrorFlag2;

[FileName,FileDir] = uigetfile('*.mat','Select the motion capture data file');

if isnumeric(FileName) || isnumeric(FileDir)
    return;
end

set(handles.f1,'Pointer','watch');

FilePath = strcat(FileDir,FileName);

[FileDir,FileNameNoExt,FileNameExt] = fileparts(FilePath);

if length(FileNameExt) == 4
    if FileNameExt == '.mat'
        set(handles.t2,'String',strcat('Data File: ',FileName));
        set(handles.t2,'ForegroundColor',[0,0,0]);
        ErrorFlag2 = false;
    else
        ErrorFlag2 = true;
    end
else
    ErrorFlag2 = true;
end

if ErrorFlag2
    set(handles.t2,'String','Error: load *.mat file only');
    set(handles.t2,'ForegroundColor',[1,0,0]);
    return;
end

if ~ErrorFlag2
    set(handles.t2,'String','Loading *.mat file ...');
    set(handles.t2,'ForegroundColor',[0,0,1]);
    pause(1);
    DataFile = load(FilePath);
    waitfor(DataFile);
    set(handles.t2,'String',strcat('Data File: ',FileName));
    set(handles.t2,'ForegroundColor',[0,0,0]);
    DataFile = struct2cell(DataFile);
    mocapstruct = DataFile{1};
end


% --- Executes on button press in P4.
function P4_Callback(hObject, eventdata, handles)
global mocapstruct ErrorFlag3;

var = uigetvariables({'Choose motion capture struct variable from MATLAB workspace'},'InputTypes',{'struct'});
if ~isempty(var) && isstruct(var{1})
    mocapstruct = var{1};
    set(handles.t3,'String',strcat('Struct Name: ',varname(mocapstruct)));
    set(handles.t3,'ForegroundColor',[0,0,0]);
    ErrorFlag3 = false;
else
    set(handles.t3,'String','Error: load struct var again');
    set(handles.t3,'ForegroundColor',[1,0,0]);
    ErrorFlag3 = true;
    return;
end

% --- isStructEmpty function to check if a nested struct is empty or not
function res = isStructEmpty(output)

rootLayerFieldNames = fields(output);
isOutPutEmpty = true;

for i = 1: numel(fieldnames(output))
    
    rootFieldName = rootLayerFieldNames{i};
    
    fieldOutput = output.(rootFieldName);
    
    fieldFieldNames = fields(fieldOutput);
    fieldNumElements = numel(fieldFieldNames);
   
    
    for j = 1:fieldNumElements
        fieldFieldSubName = fieldFieldNames{j};
        fielfFiledSubElement = fieldOutput.(fieldFieldSubName);
        
        if ~isempty(fielfFiledSubElement)
            
            isOutPutEmpty = false;
            res = isOutPutEmpty;
            return;
        end
    end
    
end

res = isOutPutEmpty;


% --- Executes on button press in p5.
function p5_Callback(hObject, eventdata, handles)
global mocapstruct ErrorFlag1 ErrorFlag2 ErrorFlag3; 
global CurrentFrame TotalFrameNo;

% Initial motion capture data load check
if ~((~ErrorFlag1 && ~ErrorFlag2) || (~ErrorFlag1 && ~ErrorFlag3))
        errordlg('Error in input struct. Load the motion caption struct first', ...
             'Motion Capture Data Check');
    return;
end

% ----- Output file operations

% Create output struct
global output;

global outputStructName outputStructPath outputStructFullPath backupDirectoryName backupDirectoryFullPath;
[outputStructName,outputStructPath] = uiputfile('*.mat','Save Output Struct As');

if ~ischar(outputStructName) && ~ischar(outputStructPath)
    if outputStructName == 0 || outputStructPath == 0
        errordlg('Error in output struct name','File Name Check');
        set(handles.p5,'Enable','On');
        return;
    end
end

% Check for valid file name 
if regexp(outputStructName, '[/\*:?"<>|]', 'once')
    errordlg('Error in output struct name','File Name Check');
    return;
end

% Check file extension

if isempty(strfind(outputStructName,'.mat'))
    errordlg('Error in output struct file name','Output File Name Check');
    return;
end

outputStructFullPath = strcat(outputStructPath,outputStructName);
disp(outputStructFullPath);

if exist(outputStructFullPath,'file')

    outputTemp = load(outputStructFullPath);
    outputTemp = outputTemp.output;
    if ~isStructEmpty(outputTemp)
        while true
            ansAppend = questdlg('The output struct exists and it is NOT empty. Do you want to append data to this file?','Motion Caputre Struct Data Operation','No');
            
            if ~isempty(ansAppend) && ~strcmp(ansAppend, 'Cancel')
                break
            end
        end
        
        switch ansAppend
            case 'Yes'
                output = outputTemp;
            case 'No'
                movefile(outputStructFullPath,strcat(outputStructPath,'OLD_',strrep(outputStructName,'.mat',''),'_',datestr(now,'mmmm-dd-yyyy-HH-MM-SS'),'.mat'));
        end
    end
end

% Save the output struct file
outputStructQuickSave(handles);

% Reset the time (Used in SaveOutPutStruct function)
global timeCheck;
timeCheck = tic;
global BackupTime;
BackupTime = 0;

% Create Backup Directoy
backupDirectoryName = strcat('Backup_',strrep(outputStructName,'.mat',''));
if ~exist(strcat(outputStructPath,backupDirectoryName),'dir')
    mkdir(outputStructPath,backupDirectoryName);
end

backupDirectoryFullPath = strcat(outputStructPath,backupDirectoryName);

% Backup the output struct file
backupOutputStruct(handles);

% Check data struct load status and make other elements visible

if (~ErrorFlag1 && ~ErrorFlag2) || (~ErrorFlag1 && ~ErrorFlag3)
    set(handles.uip3,'Visible','On');
    
    % Plot 1st Frame
    CurrentFrame = 1;
    AnimateMarkerMovie(mocapstruct, CurrentFrame, hObject, eventdata, handles);
    
    % Get Total Frame Number
    TotalFrameNo = length(mocapstruct.markers_preproc.HeadF);
    
    % Write TotalFrameNo in GUI
    set(handles.t6,'String',num2str(TotalFrameNo));
    
    % Let's diable Load Excel adn Load Data (File & Workspace) Buttons
    % set(handles.p2,'Enable','Off'); % Decision left up to user
    set(handles.p3,'Enable','Off');
    set(handles.P4,'Enable','Off');
    set(handles.p5,'Enable','Off');
else
    errordlg('Error in input struct. Load the motion caption struct first', ...
             'Motion Capture Data Check');
    return;
end


% --- AddStructFieldFun function
function AddStructFieldFun(str,b,e,field,handles)
global output;
if  ~isempty(str)
    data = get(b,'String');
    data{end+1} = str;
    set(b,'String',data);
    set(e,'String','');
    
    % Update output struct
    output.(field).(str) = [];

    % Save the output struct file
    outputStructQuickSave(handles);
end

% --- varname function
function out = varname(~)
  out = inputname(1);
  
  
% --- SetSliderSetting function
function SetSliderSetting(hObject, eventdata, handles)
ArrowStep = str2double(get(handles.es3,'String'));
AreaStep  = str2double(get(handles.es4,'String'));
minVal = get(handles.s1,'Min');
maxVal = get(handles.s1,'Max');
CurrentRange = get(handles.s1,'SliderStep');
CurrentRange(1) = ArrowStep/(maxVal-minVal);
CurrentRange(2) = AreaStep/(maxVal-minVal);
set(handles.s1,'SliderStep',CurrentRange);


% --- CheckVariableName function
function [strOut, errorStatus] = CheckVariableName(e, str)

errorStatus = false;
strOut = [];

% Remove any non-letter, non-digit, and white space characters
str(~isstrprop(str,'alphanum')) = [];
str(isstrprop(str,'wspace')) = [];

% Remove any trailing non-letter character
if ~isempty(str)
    while ~isletter(str(1))
        str(1) = [];
        % Break if empty
        if isempty(str)
            break;
        end
    end
end

% Return if string is empty
if isempty(str)
    set(e,'String',[]);
    errorStatus = true;
    errordlg('Error in struct field name, Enter again !','Input Variable Check');
    return;
else
    errorStatus = false;
    strOut = str;
end

% --- Executes on slider movement.
function s1_Callback(hObject, eventdata, handles)

val = get(handles.s1,'Value');
set(handles.s1,'Value',round(val));

% Play the selected frame
global mocapstruct CurrentFrame AutoSaveFlag LoadFrameFlag FrameArray;

CurrentFrame = round(val);

if LoadFrameFlag
    CurrentFrame = FrameArray(round(val));
end   

if AutoSaveFlag
    SetFrameProperties(handles);
end
    
AnimateMarkerMovie(mocapstruct, CurrentFrame, hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function s1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function es1_Callback(hObject, eventdata, handles)
Val = str2double(get(handles.es1,'String'));

maxVal = str2double(get(handles.es2,'String'));

if Val  >=  maxVal
    Val = max( round(maxVal/2) , 1);
    set(handles.es1,'String',num2str(Val));
end

set(handles.s1,'Min',Val);
set(handles.s1,'Value',Val);
SetSliderSetting(hObject, eventdata, handles);

s1_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function es1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function es2_Callback(hObject, eventdata, handles)
global TotalFrameNo;
Val = str2double(get(handles.es2,'String'));

minVal = str2double(get(handles.es1,'String'));

if Val  <=  minVal
    Val = min( round(minVal*2) , TotalFrameNo );
    set(handles.es2,'String',num2str(Val));
end

if Val  <=  get(handles.s1,'Value')
    set(handles.s1,'Value',Val);
    s1_Callback(hObject, eventdata, handles);
end

set(handles.s1,'Max',Val);
SetSliderSetting(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function es2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function es3_Callback(hObject, eventdata, handles)
Val = str2double(get(handles.es3,'String'));
AreaStep = str2double(get(handles.es4,'String'));

if Val >= AreaStep
    Val = AreaStep;
    set(handles.es3,'String',num2str(Val));
end

CurrentRange = get(handles.s1,'SliderStep');
minVal = get(handles.s1,'Min');
maxVal = get(handles.s1,'Max');
CurrentRange(1) = Val/(maxVal-minVal);
set(handles.s1,'SliderStep',CurrentRange);

% --- Executes during object creation, after setting all properties.
function es3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function es4_Callback(hObject, eventdata, handles)
Val = str2double(get(handles.es4,'String'));
ArrowStep = str2double(get(handles.es3,'String'));

if Val <= ArrowStep
    Val = ArrowStep;
    set(handles.es4,'String',num2str(Val));
end

CurrentRange = get(handles.s1,'SliderStep');
minVal = get(handles.s1,'Min');
maxVal = get(handles.s1,'Max');
CurrentRange(2) = Val/(maxVal-minVal);
set(handles.s1,'SliderStep',CurrentRange);

% --- Executes during object creation, after setting all properties.
function es4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tg1.
function tg1_Callback(hObject, eventdata, handles)
global CurrentFrame AutoSaveFlag;

% Auto Animate
SliderVal = get(handles.s1,'Value');
maxVal = get(handles.s1,'Max');
minVal = get(handles.s1,'Min');
CurrentRange = get(handles.s1,'SliderStep');

if SliderVal == maxVal
    set(handles.tg1,'Value', 0);
    return;
end

set(handles.tg1,'String', 'Pause');
set(handles.tg1,'ForegroundColor', [1 0 0]);

ArrowStep = CurrentRange(1) * (maxVal-minVal);

for CurrentFrame = SliderVal:ArrowStep:maxVal
    set(handles.s1,'Value',CurrentFrame);
    s1_Callback(hObject, eventdata, handles);
    
    pause(1/ArrowStep); % Pause frequency
    
    if gco(handles.f1) ~= handles.tg1 || get(handles.tg1,'Value') == 0
        set(handles.tg1,'Value', 0);
        set(handles.tg1,'String', '>');
        set(handles.tg1,'ForegroundColor', [0 0 1]);
        
        % Let's Cancel the Auto Save
        set(handles.c1,'Value',0);
        set(handles.c1,'ForegroundColor',[0 0 0]);
        AutoSaveFlag = get(handles.c1,'Value');
        
        return;
    end
end

% End of animation, reset tg1
set(handles.tg1,'Value', 0);
set(handles.tg1,'String', '>');
set(handles.tg1,'ForegroundColor', [0 0 1]);

% Let's Cancel the Auto Save
set(handles.c1,'Value',0);
set(handles.c1,'ForegroundColor',[0 0 0]);
AutoSaveFlag = get(handles.c1,'Value');


% --- Executes on button press in tg2.
function tg2_Callback(hObject, eventdata, handles)
global CurrentFrame AutoSaveFlag;

% Auto Animate
SliderVal = get(handles.s1,'Value');
maxVal = get(handles.s1,'Max');
minVal = get(handles.s1,'Min');
CurrentRange = get(handles.s1,'SliderStep');

if SliderVal == maxVal
    set(handles.tg2,'Value', 0);
    return;
end

set(handles.tg2,'String', 'Pause');
set(handles.tg2,'ForegroundColor', [1 0 0]);

AreaStep = CurrentRange(2) * (maxVal-minVal);

for CurrentFrame = SliderVal:AreaStep:maxVal
    set(handles.s1,'Value',CurrentFrame);
    s1_Callback(hObject, eventdata, handles);
    
    pause(1/AreaStep); % Pause frequency
    
    if gco(handles.f1) ~= handles.tg2 || get(handles.tg2,'Value') == 0
        set(handles.tg2,'Value', 0);
        set(handles.tg2,'String', '>>');
        set(handles.tg2,'ForegroundColor', [0 0 1]);
        
        % Let's Cancel the Auto Save
        set(handles.c1,'Value',0);
        set(handles.c1,'ForegroundColor',[0 0 0]);
        AutoSaveFlag = get(handles.c1,'Value');
        
        return;
    end
end

% End of animation, reset tg2
set(handles.tg2,'Value', 0);
set(handles.tg2,'String', '>>');
set(handles.tg2,'ForegroundColor', [0 0 1]);

% Let's Cancel the Auto Save
set(handles.c1,'Value',0);
set(handles.c1,'ForegroundColor',[0 0 0]);
AutoSaveFlag = get(handles.c1,'Value');

% --- Executes on button press in tg3.
function tg3_Callback(hObject, eventdata, handles)
global CurrentFrame AutoSaveFlag;

% Auto Animate
SliderVal = get(handles.s1,'Value');
maxVal = get(handles.s1,'Max');
minVal = get(handles.s1,'Min');
CurrentRange = get(handles.s1,'SliderStep');

if SliderVal == minVal
    set(handles.tg3,'Value', 0);
    return;
end

set(handles.tg3,'String', 'Pause');
set(handles.tg3,'ForegroundColor', [1 0 0]);

ArrowStep = CurrentRange(1) * (maxVal-minVal);

for CurrentFrame = SliderVal:-ArrowStep:minVal
    set(handles.s1,'Value',CurrentFrame);
    s1_Callback(hObject, eventdata, handles);
    
    pause(1/ArrowStep); % Pause frequency
    
    if gco(handles.f1) ~= handles.tg3 || get(handles.tg3,'Value') == 0
        set(handles.tg3,'Value', 0);
        set(handles.tg3,'String', '<');
        set(handles.tg3,'ForegroundColor', [0 0 1]);
        
        % Let's Cancel the Auto Save
        set(handles.c1,'Value',0);
        set(handles.c1,'ForegroundColor',[0 0 0]);
        AutoSaveFlag = get(handles.c1,'Value');
        
        return;
    end
end

% End of animation, reset tg3
set(handles.tg3,'Value', 0);
set(handles.tg3,'String', '<');
set(handles.tg3,'ForegroundColor', [0 0 1]);

% Let's Cancel the Auto Save
set(handles.c1,'Value',0);
set(handles.c1,'ForegroundColor',[0 0 0]);
AutoSaveFlag = get(handles.c1,'Value');


% --- Executes on button press in tg4.
function tg4_Callback(hObject, eventdata, handles)
global CurrentFrame AutoSaveFlag;

% Auto Animate
SliderVal = get(handles.s1,'Value');
maxVal = get(handles.s1,'Max');
minVal = get(handles.s1,'Min');
CurrentRange = get(handles.s1,'SliderStep');

if SliderVal == minVal
    set(handles.tg4,'Value', 0);
    return;
end

set(handles.tg4,'String', 'Pause');
set(handles.tg4,'ForegroundColor', [1 0 0]);

AreaStep = CurrentRange(2) * (maxVal-minVal);

for CurrentFrame = SliderVal:-AreaStep:minVal
    set(handles.s1,'Value',CurrentFrame);
    s1_Callback(hObject, eventdata, handles);
    
    pause(1/AreaStep); % Pause frequency
    
    if gco(handles.f1) ~= handles.tg4 || get(handles.tg4,'Value') == 0
        set(handles.tg4,'Value', 0);
        set(handles.tg4,'String', '<<');
        set(handles.tg4,'ForegroundColor', [0 0 1]);
        
        % Let's Cancel the Auto Save
        set(handles.c1,'Value',0);
        set(handles.c1,'ForegroundColor',[0 0 0]);
        AutoSaveFlag = get(handles.c1,'Value');
        
        return;
    end
end

% End of animation, reset tg4
set(handles.tg4,'Value', 0);
set(handles.tg4,'String', '<<');
set(handles.tg4,'ForegroundColor', [0 0 1]);

% Let's Cancel the Auto Save
set(handles.c1,'Value',0);
set(handles.c1,'ForegroundColor',[0 0 0]);
AutoSaveFlag = get(handles.c1,'Value');


% --- Executes on button press in p6.
function p6_Callback(hObject, eventdata, handles)
SetFrameProperties(handles);


% --- SetFrameProperties
function SetFrameProperties(handles)
global CurrentFrame output numCols;

HandlesArray = [handles.b2, handles.b1, handles.b3, handles.b5 handles.b6 ...
                handles.b7  handles.b8];
            
for i = 1:numCols
    % Overwrite any previous occurance
    CheckOverWrtie(HandlesArray(i));
    
    % Update the struct with new frame number
    UpdateOutputStruct(HandlesArray(i));
end

SaveOutputStruct(handles);


% --- SaveOutputStruct
function SaveOutputStruct(handles)
global output BackupTime;
global outputStructName outputStructFullPath backupDirectoryFullPath;
global saveOutputStructTimeFreq backupOutputStructTimeFreq;
global timeCheck;

if floor(toc(timeCheck)) >= saveOutputStructTimeFreq   % Save every saveOutputStructTimeFreq seconds
    save(outputStructFullPath, 'output');
    set(handles.t14,'String',strcat( 'Last Time Output File Saved:',{' '}, datestr(now,'dddd, mmmm dd, HH:MM:SS') ));
    
    % Add elapsedTime to BackupTime 
    BackupTime = BackupTime + floor(toc(timeCheck));
    
    % Reset timer
    timeCheck = tic;
end

if BackupTime >= backupOutputStructTimeFreq % Backup every backupOutputStructTimeFreq seconds
    
    if exist(outputStructFullPath,'file') && exist(backupDirectoryFullPath,'dir')
        
        if ismac || isunix
            copyfile(outputStructFullPath,strcat(backupDirectoryFullPath,'/Backup_',outputStructName,'_',datestr(now,'mmmm-dd-yyyy-HH-MM-SS'),'.mat'));
        elseif ispc
            copyfile(outputStructFullPath,strcat(backupDirectoryFullPath,'\Backup_',outputStructName,'_',datestr(now,'mmmm-dd-yyyy-HH-MM-SS'),'.mat'));
        end
    end
    
    % Reset the backup time
    BackupTime = 0;
    
end    
    

% --- UpdateOutputStruct
function UpdateOutputStruct(b)
global CurrentFrame output;

Data  = get(b,'String');
Value = get(b,'Value');


if Value > 2
    if isempty(find(output.(Data{1}).(Data{Value}) == CurrentFrame, 1))
        output.(Data{1}).(Data{Value})(end+1) = CurrentFrame;
    end
end


% --- CheckOverWrtie
function CheckOverWrtie(b)
global CurrentFrame output;

Data  = get(b,'String');

% Get the length of selected column
Len = length(Data);

% Check for occurance of CurrentFrame
for i = 3:Len
    
    % Get the location of CurrentFrame if any
    Loc = find( output.(Data{1}).(Data{i}) == CurrentFrame, 1);
    
    % If Loc is not empty, remove it from the struct.
    if ~isempty(Loc)
        output.(Data{1}).(Data{i})(Loc) = [];
        return;
    end
end

% --- Executes on button press in c1.
function c1_Callback(hObject, eventdata, handles)
global AutoSaveFlag;
AutoSaveFlag = get(handles.c1,'Value');

if AutoSaveFlag
    set(handles.c1,'ForegroundColor',[0 0 1]);
else
    set(handles.c1,'ForegroundColor',[0 0 0]);
end

% --- Executes on button press in p7.
function p7_Callback(hObject, eventdata, handles)
global numCols;
HandlesArray = [handles.b2, handles.b1, handles.b3, handles.b5 handles.b6 ...
                handles.b7  handles.b8];
            
for i = 1:numCols
    set( HandlesArray(i) , 'Value' , 3 );
end         
    


% --- Executes on button press in p8.
function p8_Callback(hObject, eventdata, handles)
global FrameArray LoadFrameFlag AutoSaveFlag;

% Let's Cancel the Auto Save
set(handles.c1,'Value',0);
set(handles.c1,'ForegroundColor',[0 0 0]);
AutoSaveFlag = get(handles.c1,'Value');

LoadFrameFlag = false;

var = uigetvariables({'Choose frame array variable from MATLAB workspace [Size > 1]'},'InputTypes',{'struct'});

if ~isempty(var)
    
    if ~iscell(var{1}) && ~isstruct(var{1}) && ~isempty(var{1}) && (numel(var{1})>1)
        
        FrameArray = var{1};

        LoadFrameFlag = true;
        
        AdjustSliderMinMaxUponFrameLoad(hObject, eventdata, handles);
        
        % Activate the "Go Back to Main Frames" button
        set(handles.p9,'Enable','On');
        
    end
    
else
    LoadFrameFlag = false;
    return;
end


% --- AdjustSliderMinMax
function AdjustSliderMinMaxUponFrameLoad(hObject, eventdata, handles)

global FrameArray;

set(handles.t12,'Visible','On');
set(handles.t13,'Visible','On');

numFrames = numel(FrameArray);

set(handles.t13,'String',num2str(numFrames));

minVal = 1;
maxVal = numFrames;

ResetSliderValues(hObject,eventdata,handles,minVal,maxVal);


% --- Executes on button press in p9.
function p9_Callback(hObject, eventdata, handles)
global AutoSaveFlag TotalFrameNo LoadFrameFlag;

set(handles.t12,'Visible','Off');
set(handles.t13,'Visible','Off');

LoadFrameFlag = false;

% Let's Cancel the Auto Save
set(handles.c1,'Value',0);
AutoSaveFlag = get(handles.c1,'Value');

minVal = 1;
maxVal = min(101, TotalFrameNo);

ResetSliderValues(hObject,eventdata,handles,minVal,maxVal);

% Let's remove "Go Back to Main Frames button"
set(handles.p9,'Enable','Off');


% --- ResetSliderValues
function ResetSliderValues(hObject,eventdata,handles,minVal,maxVal)

ArrowStep = max(minVal, 1 );
AreaStep =  min(maxVal, 10 );

CurrentRange = [];
CurrentRange(1) = ArrowStep/(maxVal-minVal);
CurrentRange(2) = AreaStep/(maxVal-minVal);

set(handles.es1,'String',num2str(minVal));
set(handles.es2,'String',num2str(maxVal));
set(handles.es3,'String',num2str(ArrowStep));
set(handles.es4,'String',num2str(AreaStep));

set(handles.s1,'Value',minVal);
set(handles.s1,'Min',minVal);
set(handles.s1,'Max',maxVal);
set(handles.s1,'SliderStep',CurrentRange);

s1_Callback(hObject, eventdata, handles);


% --- Executes when user attempts to close f1.
function f1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);


function outputStructQuickSave(handles)
global output;
global outputStructFullPath;
global timeCheck;

save(outputStructFullPath, 'output');
set(handles.t14,'String',strcat('Last Time Output File Saved:',{' '}, datestr(now,'dddd, mmmm dd, HH:MM:SS') ));

% Reset timer
timeCheck = tic;


function backupOutputStruct(handles)
global BackupTime;
global outputStructName outputStructFullPath backupDirectoryFullPath backupDirectoryName outputStructPath;

% Create Backup Directoy
if ~exist(strcat(outputStructPath,backupDirectoryName),'dir')
    mkdir(outputStructPath,backupDirectoryName);
end

if exist(outputStructFullPath,'file') && exist(backupDirectoryFullPath,'dir')
    
    if ismac || isunix
        copyfile(outputStructFullPath,strcat(backupDirectoryFullPath,'/',strrep(outputStructName,'.mat',''),'_',datestr(now,'mmmm-dd-yyyy-HH-MM-SS'),'.mat'));
    elseif ispc
        copyfile(outputStructFullPath,strcat(backupDirectoryFullPath,'\',strrep(outputStructName,'.mat',''),'_',datestr(now,'mmmm-dd-yyyy-HH-MM-SS'),'.mat'));
    end
end

% Reset the backup time
BackupTime = 0;



function eBase_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function eBase_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function SaveMenu_Callback(hObject, eventdata, handles)
outputStructQuickSave(handles);

% --------------------------------------------------------------------
function BackupMenu_Callback(hObject, eventdata, handles)
backupOutputStruct(handles);


% --------------------------------------------------------------------
function ExitMenu_Callback(hObject, eventdata, handles)

ansClose = questdlg('Exit FATGUI application?','FATGUI Application','No');

if strcmp(ansClose, 'Yes')
    outputStructQuickSave(handles);
    backupOutputStruct(handles);
    close(handles.f1);
end


% --------------------------------------------------------------------
function LoadMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function ExcelMenu_Callback(hObject, eventdata, handles)
global ErrorFlag1 output outputFieldNames numCols;
[FileName,FileDir] = uigetfile('*.xlsx','Select the Excel file');

if isnumeric(FileName) || isnumeric(FileDir)
    return;
end

FilePath = strcat(FileDir,FileName);

[FileDir,FileNameNoExt,FileNameExt] = fileparts(FilePath);

if length(FileNameExt) == 5
    if FileNameExt == '.xlsx'
        set(handles.t1,'String',strcat('Excel File: ',FileName));
        set(handles.t1,'ForegroundColor',[0,0,0]);
        ErrorFlag1 = false;
    else
        ErrorFlag1 = true;
    end
else
    ErrorFlag1 = true;
end

if ErrorFlag1
    set(handles.t1,'String','Error: load *.xlsx file only');
    set(handles.t1,'ForegroundColor',[1,0,0]);
    return;
end

ReadExcelFlag = false;

numCols = 7;

if ~ErrorFlag1
    [num,txt] = xlsread(FilePath);
    
    SizTXT = size(txt);
    
    if SizTXT(2) ~=numCols
        set(handles.t1,'String','Error: wrong number of columns');
        set(handles.t1,'ForegroundColor',[1,0,0]);
        return;
    end
    
    numRows = SizTXT(1);
    ReadExcelFlag = true;
    
end   

if ReadExcelFlag
    
    ErrorFlag1 = false;
    
    handlesMAT = [handles.b2 handles.b1 handles.b3 handles.b5 handles.b6 handles.b7 handles.b8];
    
    for j = 1:numCols
        
        % Clear column j
        set(handlesMAT(j),'String','');
        % Read column j
        col = txt(:,j);
        for i = 1:numRows
            if isempty(col{i})
                break; % goto next column, goto XX
            else
                data = get(handlesMAT(j),'String');
                data{end+1} = col{i};
                set(handlesMAT(j),'String',data);
                if i == 1
                    data = get(handlesMAT(j),'String');
                    data{end+1} = '-------';
                    set(handlesMAT(j),'String',data);
                    
                    % Build output struct variable
                    output.(col{i}) = [];
                else
                    % Build output struct variable
                    output.(col{1}).(col{i}) = [];
                end
            end
        end
        
        % XX: break comes here
        
        set(handlesMAT(j),'Value',3);
        
    end
    outputFieldNames = fieldnames(output);
   
else
    ErrorFlag1 = true;
end


% --------------------------------------------------------------------
function DataFileMenu_Callback(hObject, eventdata, handles)
global mocapstruct ErrorFlag2;

[FileName,FileDir] = uigetfile('*.mat','Select the motion capture data file');

if isnumeric(FileName) || isnumeric(FileDir)
    return;
end

set(handles.f1,'Pointer','watch');

FilePath = strcat(FileDir,FileName);

[FileDir,FileNameNoExt,FileNameExt] = fileparts(FilePath);

if length(FileNameExt) == 4
    if FileNameExt == '.mat'
        set(handles.t2,'String',strcat('Data File: ',FileName));
        set(handles.t2,'ForegroundColor',[0,0,0]);
        ErrorFlag2 = false;
    else
        ErrorFlag2 = true;
    end
else
    ErrorFlag2 = true;
end

if ErrorFlag2
    set(handles.t2,'String','Error: load *.mat file only');
    set(handles.t2,'ForegroundColor',[1,0,0]);
    return;
end

if ~ErrorFlag2
    set(handles.t2,'String','Loading *.mat file ...');
    set(handles.t2,'ForegroundColor',[0,0,1]);
    pause(1);
    DataFile = load(FilePath);
    waitfor(DataFile);
    set(handles.t2,'String',strcat('Data File: ',FileName));
    set(handles.t2,'ForegroundColor',[0,0,0]);
    DataFile = struct2cell(DataFile);
    mocapstruct = DataFile{1};
end


% --------------------------------------------------------------------
function DataWorkSpaceMenu_Callback(hObject, eventdata, handles)
global mocapstruct ErrorFlag3;

var = uigetvariables({'Choose motion capture struct variable from MATLAB workspace'},'InputTypes',{'struct'});
if ~isempty(var) && isstruct(var{1})
    mocapstruct = var{1};
    set(handles.t3,'String',strcat('Struct Name: ',varname(mocapstruct)));
    set(handles.t3,'ForegroundColor',[0,0,0]);
    ErrorFlag3 = false;
else
    set(handles.t3,'String','Error: load struct var again');
    set(handles.t3,'ForegroundColor',[1,0,0]);
    ErrorFlag3 = true;
    return;
end


% --------------------------------------------------------------------
function GoMenu_Callback(hObject, eventdata, handles)
global mocapstruct ErrorFlag1 ErrorFlag2 ErrorFlag3; 
global CurrentFrame TotalFrameNo;

% Initial motion capture data load check
if ~((~ErrorFlag1 && ~ErrorFlag2) || (~ErrorFlag1 && ~ErrorFlag3))
        errordlg('Error in input struct. Load the motion caption struct first', ...
             'Motion Capture Data Check');
    return;
end

% ----- Output file operations

% Create output struct
global output;

global outputStructName outputStructPath outputStructFullPath backupDirectoryName backupDirectoryFullPath;
[outputStructName,outputStructPath] = uiputfile('*.mat','Save Output Struct As');

if ~ischar(outputStructName) && ~ischar(outputStructPath)
    if outputStructName == 0 || outputStructPath == 0
        errordlg('Error in output struct name','File Name Check');
        set(handles.p5,'Enable','On');
        return;
    end
end

% Check for valid file name 
if regexp(outputStructName, '[/\*:?"<>|]', 'once')
    errordlg('Error in output struct name','File Name Check');
    return;
end

% Check file extension

if isempty(strfind(outputStructName,'.mat'))
    errordlg('Error in output struct file name','Output File Name Check');
    return;
end

outputStructFullPath = strcat(outputStructPath,outputStructName);
disp(outputStructFullPath);

if exist(outputStructFullPath,'file')

    outputTemp = load(outputStructFullPath);
    outputTemp = outputTemp.output;
    if ~isStructEmpty(outputTemp)
        while true
            ansAppend = questdlg('The output struct exists and it is NOT empty. Do you want to append data to this file?','Motion Caputre Struct Data Operation','No');
            
            if ~isempty(ansAppend) && ~strcmp(ansAppend, 'Cancel')
                break
            end
        end
        
        switch ansAppend
            case 'Yes'
                output = outputTemp;
            case 'No'
                movefile(outputStructFullPath,strcat(outputStructPath,'OLD_',strrep(outputStructName,'.mat',''),'_',datestr(now,'mmmm-dd-yyyy-HH-MM-SS'),'.mat'));
        end
    end
end

% Save the output struct file
outputStructQuickSave(handles);

% Reset the time (Used in SaveOutPutStruct function)
global timeCheck;
timeCheck = tic;
global BackupTime;
BackupTime = 0;

% Create Backup Directoy
backupDirectoryName = strcat('Backup_',strrep(outputStructName,'.mat',''));
if ~exist(strcat(outputStructPath,backupDirectoryName),'dir')
    mkdir(outputStructPath,backupDirectoryName);
end

backupDirectoryFullPath = strcat(outputStructPath,backupDirectoryName);

% Backup the output struct file
backupOutputStruct(handles);

% Check data struct load status and make other elements visible

if (~ErrorFlag1 && ~ErrorFlag2) || (~ErrorFlag1 && ~ErrorFlag3)
    set(handles.uip3,'Visible','On');
    
    % Plot 1st Frame
    CurrentFrame = 1;
    AnimateMarkerMovie(mocapstruct, CurrentFrame, hObject, eventdata, handles);
    
    % Get Total Frame Number
    TotalFrameNo = length(mocapstruct.markers_preproc.HeadF);
    
    % Write TotalFrameNo in GUI
    set(handles.t6,'String',num2str(TotalFrameNo));
    
    % Let's diable Load Excel adn Load Data (File & Workspace) Buttons
    % set(handles.p2,'Enable','Off'); % Decision left up to user
    set(handles.p3,'Enable','Off');
    set(handles.P4,'Enable','Off');
    set(handles.p5,'Enable','Off');
else
    errordlg('Error in input struct. Load the motion caption struct first', ...
             'Motion Capture Data Check');
    return;
end
