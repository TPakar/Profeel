function varargout = Profeel(varargin)

%%
% Copyright (C) 2020 Tomppa Pakarinen, tomppa.pakarinen@pshp.fi


% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/
%%
% PROFEEL MATLAB code for Profeel.fig
%      PROFEEL, by itself, creates a new PROFEEL or raises the existing
%      singleton*.
%
%      H = PROFEEL returns the handle to a new PROFEEL or the handle to
%      the existing singleton*.
%
%      PROFEEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROFEEL.M with the given input arguments.
%
%      PROFEEL('Property','Value',...) creates a new PROFEEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Profeel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Profeel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Profeel

% Last Modified by GUIDE v2.5 11-Jun-2020 19:41:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Profeel_OpeningFcn, ...
                   'gui_OutputFcn',  @Profeel_OutputFcn, ...
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

% --- Executes just before Profeel is made visible.
function Profeel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Profeel (see VARARGIN)

%charencoding = slCharacterEncoding();
%disp(['Character encoding: ', charencoding]);

% Get screen dimensions
handles.screendims = get(groot, 'ScreenSize');
screendims2(1) = handles.screendims(3)*0.15/2;
screendims2(2) = handles.screendims(4)*0.15/2;
screendims2(3:4) = handles.screendims(3:4)*0.85;
newscreenwidth = screendims2(4);
% Change the GUI positions and size (Size 85% from screen size and centered)
set(gcf, 'Position', screendims2);

% Aspect ratio for the screen used in design
defaultscreenwidth = 1224;


screenwidthratio = newscreenwidth/defaultscreenwidth;


% Choose default command line output for Profeel
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Profeel.

% Initializations
handles.normalizationtype = 'NON';
set(handles.popupmenuNormalisointi, 'Value', 2);
handles.residualcount = 0;
handles.residualState = 0;
handles.dtacriterion = 3; % 3mm as default
handles.dosecriterion = 3; % 3% as default
handles.order = [];
handles.interpres = [];
handles.meanminval = -14.3;
handles.meanmaxval = 14.3;
handles.exportables = struct;
handles.exportcount = 0;
handles.origlist = {};
handles.listboxArray = {};
handles.alldata = struct;
handles.datacount = 1;
handles.newlinecount = 1;
handles.nogammacalc = 1;
handles.caxcorrection = 1;
handles.showError = 1;
handles.redraw = 0;
handles.shiftaxis = 0;
handles.lowresvisualization = 1;
handles.depthscroll = 0;
handles.DTAbruteforce = 0;
handles.DTAanalytical = 1;
handles.gammaresults = struct;
handles.gammaresultcount = 1;
handles.gammaMinimumdose = [];
handles.isodosecount = 1;
handles.normalizationtoval = 20;
handles.omnicount = 1;
handles.reffixed = 0;
handles.pddColNames = {'R100 [cm]', 'R80 [cm]', 'R50 [cm]', 'D10 [%]', 'D20 [%]', 'J100/200', '', '', ''};
handles.profileColNames = {'CAX dev [cm]', 'Sym [%]', 'Homog. [%]', 'Dmax [%]', 'Dmin [%]', 'Dev[%]', 'FS [cm]', 'PenR [cm]', 'PenL [cm]'};
handles.symmetryThresHold = 300;
handles.penumbraThresHoldRatio = 2;
handles.PTWcount = 0;

if strcmp(get(hObject,'Visible'),'off')

end


% Change text fonts to normalized
txthandle = findall(gcf, '-property', 'FontUnits');
set(txthandle, 'FontUnits', 'normalized');

txthandle = findall(gcf, '-property', 'FontSize');

for i = 1:length(txthandle)
    currentFsize = get(txthandle(i), 'FontSize');
    set(txthandle(i), 'FontSize', currentFsize*screenwidthratio);
end

% UIWAIT makes Profeel wait for user response (see UIRESUME)
% uiwait(handles.figure1);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Profeel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on button press in pushbuttonOmni.
function pushbuttonOmni_Callback(hObject, eventdata, handles)

% Get Omnipro data folder
handles = guidata(hObject);
handles.omnidir = handles.mydir;

% Get OmniPro data. Collects structs inside the struct as: omni1, omni2, omni3, ... , omniN
%----->
%omnidata.omni3.name, omnidata.omni3.RadiationType, omnidata.omni3.Energy, 
%omnidata.omni3.SSD, omnidata.omni3.SID, omnidata.omni3.FieldSizeCr,
%omnidata.omni3.FieldSizeIn, omnidata.omni3.DataType, omnidata.omni3.DataFactor,
%omnidata.omni3.DataUnit, omnidata.omni3.LengthUnit, omnidata.omni3.Plane, 
%omnidata.omni3.xaxis, omnidata.omni3.yaxis, omnidata.omni3.data, omnidata.omni3.xpos,
%omnidata.omni3.ypos


handles.normalizationdist = str2double(get(handles.edit2, 'String'));
handles.normalizationperc = str2double(get(handles.edit6, 'String'));

handles.omnidata = collectOmnipro(handles.omnidir, handles.normalizationdist, handles.normalizationperc, handles.caxcorrection);

% Add names to the list
axes(handles.axes1);
cla reset;

if ~isempty(handles.listboxArray)
    for i = 1:length(fieldnames(handles.omnidata))
        handles.listboxArray{end+1} = ['omni', num2str(handles.omnicount), '_', handles.omnidata.(['omni',num2str(i)]).name];
        handles.origlist{end+1} = ['omni', num2str(handles.omnicount), '_', handles.omnidata.(['omni',num2str(i)]).name];
        handles.idx(end+1) = 0;
        handles.omnicount = handles.omnicount+1;
    end
    
else
    for i = 1:length(fieldnames(handles.omnidata))
        handles.listboxArray{i} = ['omni', num2str(handles.omnicount), '_', handles.omnidata.(['omni',num2str(i)]).name];
        handles.origlist = handles.listboxArray;
        handles.omnicount = handles.omnicount+1;
    end
    handles.idx = zeros(1,length(fieldnames(handles.omnidata)));
    handles.lastidx = handles.idx;
    
    
end


set(handles.listboxChooseData, 'String', handles.listboxArray);
%set(handles.Omnifolder, 'String', handles.omnidir(1).folder);

% Vector holding the order of chosen data
handles.order = zeros(1,length(handles.idx));


guidata(hObject, handles);







% --- Executes on button press in pushbuttonPTW.
function pushbuttonPTW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPTW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get PTW data folder
handles = guidata(hObject);
handles.PTWdir = handles.mydir;

set(handles.PTWfolder, 'String', handles.PTWdir(1).folder);
% Get manual normalization values
handles.normalizationdist = str2double(get(handles.edit2, 'String'));
handles.normalizationperc = str2double(get(handles.edit6, 'String'));
guidata(hObject,handles);
% Fetch PTW data
handles.PTWdata = collectPTW(handles.PTWdir, handles.normalizationdist, handles.normalizationperc);
% Reset GUI axes
axes(handles.axes1);
cla reset;

% Try if previous data has been imported (Replace try/catch)
if ~isempty(handles.origlist)
    for i = 1:length(fieldnames(handles.PTWdata))
        templist = ['PTW', num2str(i + handles.PTWcount), '_', handles.PTWdata.(['PTW', num2str(i)]).name];
        temporig = ['PTW', num2str(i + handles.PTWcount), '_', handles.PTWdata.(['PTW', num2str(i)]).name];
        handles.listboxArray{end+1} = strcat(templist);
        handles.origlist{end+1} = strcat(temporig);
        
        try
            [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,2),...
                handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,1), 'PTW', handles.caxcorrection);
            handles.PTWdata.(['PTW', num2str(i)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
            handles.PTWdata.(['PTW', num2str(i)]).datatype = 'profile';
            if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.PTWdata.(['PTW', num2str(i)]).params)) > 0
               [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,2),...
                handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,1));
                handles.PTWdata.(['PTW', num2str(i)]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
                handles.PTWdata.(['PTW', num2str(i)]).datatype = 'pdd';
            end
        catch
            [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,2),...
                handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,1));
                handles.PTWdata.(['PTW', num2str(i)]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
                handles.PTWdata.(['PTW', num2str(i)]).datatype = 'pdd';
        end
        handles.idx(end+1) = 0;
    end
    
% No previous data -> import all from start
else
    for i = 1:length(fieldnames(handles.PTWdata))
        templist = ['PTW', num2str(i + handles.PTWcount), '_', handles.PTWdata.(['PTW', num2str(i)]).name];
        handles.listboxArray{i} = templist;
        
        try
            [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,2), handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,1), 'PTW', handles.caxcorrection);
            handles.PTWdata.(['PTW', num2str(i)]).params = round([caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL],4);
            handles.PTWdata.(['PTW', num2str(i)]).datatype = 'profile';
            if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.PTWdata.(['PTW', num2str(i)]).params)) > 0
               [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,2),...
                handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,1));
                handles.PTWdata.(['PTW', num2str(i)]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
                handles.PTWdata.(['PTW', num2str(i)]).datatype = 'pdd';
            end
        catch
            [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,2),...
                handles.PTWdata.(['PTW', num2str(i)]).interpolated(:,1));
                handles.PTWdata.(['PTW', num2str(i)]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
            handles.PTWdata.(['PTW', num2str(i)]).datatype = 'pdd';
        end
    end
    handles.origlist = handles.listboxArray;
    handles.idx = zeros(1,length(fieldnames(handles.PTWdata)));
end
% Update the PTW dataset counter
handles.PTWcount = i + handles.PTWcount;

% Add files to listbox array
set(handles.listboxChooseData, 'String', handles.listboxArray);

handles.order = zeros(1,length(handles.idx));


% Update object and object handles
guidata(hObject, handles);



% --- Executes on button press in PTWX.
function PTWX_Callback(hObject, eventdata, handles)
% hObject    handle to PTWX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hint: get(hObject,'Value') returns toggle state of PTWX


% --- Executes on button press in PTWY.
function PTWY_Callback(hObject, eventdata, handles)
% hObject    handle to PTWY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PTWY


% --- Executes on button press in MatrixX.
function MatrixX_Callback(hObject, ~, handles)
% hObject    handle to MatrixX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MatrixX


% --- Executes on button press in MatrixY.
function MatrixY_Callback(hObject, eventdata, handles)
% hObject    handle to MatrixY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MatrixY


% --- Executes on selection change in popupmenuShowData.
function popupmenuShowData_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuShowData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
currentval_dropbox = contents{get(hObject,'Value')};

switch currentval_dropbox
    case 'X/Y profiilit'
        disp('X/Y profiles chosen');
        todelete = find(contains(handles.listboxArray, 'XY axis'));
        handles.listboxArray(todelete) = [];
        handles.origlist(todelete) = [];
        handles.idx = zeros(1,length(handles.listboxArray));
    case 'Koko matriisi'
        disp('Whole matrix chosen');
        todelete = find(contains(handles.listboxArray, 'X axis') | contains(handles.listboxArray, 'Y axis'));
        handles.listboxArray(todelete) = [];
        handles.origlist(todelete) = [];
        handles.idx = zeros(1,length(handles.listboxArray));
    case 'Matriisi diagonaali'
        disp('Diagonal matrix chosen');
        todelete = find(contains(handles.listboxArray, 'diag axis'));
        handles.listboxArray(todelete) = [];
        handles.origlist(todelete) = [];
        handles.idx = zeros(1,length(handles.listboxArray));
    otherwise
        disp('Errahh');
end

set(handles.listboxChooseData, 'String', handles.listboxArray);

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuShowData contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuShowData


% --- Executes during object creation, after setting all properties.
function popupmenuShowData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuShowData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxChooseData.
function listboxChooseData_Callback(hObject, ~, handles)
% hObject    handle to listboxChooseData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Don't update chose idx list if only re-drawing the data
if handles.redraw == 0
    contents = cellstr(get(hObject,'String'));

    handles.chosen = contents{get(hObject,'Value')};
    tempch = strsplit(handles.chosen, '_');
    handles.chosen = tempch{1};
    idx1 = double(contains(handles.listboxArray, [handles.chosen, '_']));
    tempidx = handles.idx + idx1;
    temparr = find(tempidx > 1);
    handles.legendlist = {};
    tempidx(temparr(:)) = 0;
    handles.idx = tempidx;
else
    handles.redraw = 0;
end

handles.legendlist = {};
axes(handles.axes1);
tempstr = get(handles.listboxChooseData, 'String');

try
    fieldnams = fieldnames(handles.alldata);
catch
    fieldnams = [];
end

% Check which value was picked first (for gamma computation).

if sum(handles.idx) > 0
    handles.order(find(~handles.order&handles.idx)) = sum(handles.idx);
else
    handles.order(1:end) = 0;
    set(handles.pushbuttonAddLine, 'Visible', 'Off');
    set(handles.pushbuttonAddPlane, 'Visible', 'Off');
    set(handles.sizeText, 'Visible', 'Off');
end

% Plot data
lastview = get(handles.axes1, 'View');
cla reset
hold on;
grid on;
pddrowcount = 1;
profilerowcount = 1;

pddrowdata = {'                                   '};
profilerowdata = {'                                   '};
pddtabledata = [];
profiletabledata = [];
handles.chosendata = {};
plotpos = [];
chosencount = 1;
% Reset gamma button visibility
set(handles.show2dgamma, 'Visible', 'Off');
set(handles.show2dDTA, 'Visible', 'Off');
set(handles.show2dDD, 'Visible', 'Off');
set(handles.quickprofile, 'Visible', 'Off');
set(handles.checkboxfixref, 'Visible', 'Off');

for i = 1:length(handles.idx)
    if handles.idx(i) == 1
        % Change the text indicating if data is chosen (in the list)
        if ~contains(handles.listboxArray{i},'*Chosen*')
            handles.listboxArray{i} = ['<HTML><B><font face="Comic Sans MS">', handles.listboxArray{i}, '    **Chosen** <B>'];
        end
        
        %3D object
            % Set gamma minimum dose edit and related text OFF as long as
            % not existing to 3D data
            set(handles.gammarefdoseEdit, 'Visible', 'Off');
            set(handles.mingammadosetxt, 'Visible', 'Off');
            set(handles.quickprofile, 'Visible', 'Off');
            set(handles.checkboxfixref, 'Visible', 'Off');
        if contains(handles.alldata.(['Data', num2str(i)]).datatype, '3D')
            ddisp = handles.alldata.(['Data',num2str(i)]).Depthdisplay;                 % Displayed depth
            draw = handles.alldata.(['Data',num2str(i)]).Depthraw;                      % Raw depth
            if handles.redraw == 0
                %set(gcf, 'WindowButtonMotionFcn', @mousepos, 'BusyAction', 'cancel');
                set(gcf, 'WindowScrollWheelFcn', @mousescroll, 'BusyAction', 'Queue')
                % Make 2D and 1D buttons visible
                set(handles.pushbuttonAddLine, 'Visible', 'On');
                set(handles.pushbuttonAddPlane, 'Visible', 'On');
            end
            % If 3D data -> plot initially mid axis XY-plane
            hold off;
            if handles.lowresvisualization == 1
                if ismember(['Displaydata', handles.normalizationtype], fieldnames(handles.alldata.(['Data', num2str(i)])))
                    plotthis = handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype])(:,:,ddisp);
                else
                    plotthis = handles.alldata.(['Data', num2str(i)]).Displaydata(:,:,ddisp);
                end
                % Display available range
                set(handles.textdepth,'String', ['Z: ' , num2str(handles.alldata.(['Data', num2str(i)]).Displayzpos(ddisp))]);
            else
                if ismember(['data', handles.normalizationtype], fieldnames(handles.alldata.(['Data', num2str(i)])))
                    plotthis = handles.alldata.(['Data', num2str(i)]).(['data', handles.normalizationtype])(:,:,draw);
                else
                    plotthis = handles.alldata.(['Data', num2str(i)]).data(:,:,draw);
                end

                % Display available range
                set(handles.textdepth,'String', ['Z: ' , num2str(handles.alldata.(['Data', num2str(i)]).zpos(draw))]);
            end
            
            
            if handles.redraw == 0
                set(handles.sizeText, 'Visible', 'On');
                
                set(handles.sizeText, 'String', ['Range [X,Y,Z] = [', num2str(round(min(handles.alldata.(['Data', num2str(i)]).xpos),2)), ':' , num2str(round(max(handles.alldata.(['Data', num2str(i)]).xpos),2)),...
                    ',', num2str(round(min(handles.alldata.(['Data', num2str(i)]).ypos),2)), ':', num2str(round(max(handles.alldata.(['Data', num2str(i)]).ypos),2)), ',', ...
                    num2str(round(min(handles.alldata.(['Data', num2str(i)]).zpos),2)), ':', num2str(round(max(handles.alldata.(['Data', num2str(i)]).zpos),2)), ']']);
            end
            
            if handles.lowresvisualization == 1
                if handles.redraw == 0
                    %mesh3d = meshgrid(handles.alldata.(['Data', num2str(i)]).Displayxpos, handles.alldata.(['Data', num2str(i)]).Displayypos);
                    handles.s = surf(handles.alldata.(['Data', num2str(i)]).Displayxpos, handles.alldata.(['Data', num2str(i)]).Displayypos,...
                        plotthis);
                else
                    set(handles.s, 'Zdata', plotthis);
                end
            else
                if handles.redraw == 0
                   % mesh3d = meshgrid(handles.alldata.(['Data', num2str(i)]).xpos, handles.alldata.(['Data', num2str(i)]).ypos);
                    handles.s = surf(handles.alldata.(['Data', num2str(i)]).xpos, handles.alldata.(['Data', num2str(i)]).ypos,...
                        plotthis);
                else
                    set(handles.s, 'Zdata', plotthis);
                end
            end
            if handles.redraw == 0
                handles.listboxArray{i} = handles.origlist{i};
                handles.listboxArray{i} = ['<HTML><B><font face="Comic Sans MS">', handles.listboxArray{i}, '    **Chosen** <B>'];
                handles.legendlist{end+1} = handles.origlist{i};

                xlabel('Distance X [cm]');
                ylabel('Distance Y [cm]');
                axis([-inf inf -inf inf]);
                if strcmp(handles.normalizationtype, 'NON')
                    Hc = colorbar;
                    if ismember('DataUnit', fieldnames(handles.alldata.(['Data', num2str(i)])))
                        title(Hc, handles.alldata.(['Data', num2str(i)]).DataUnit);
                    end
                else
                    Hc = colorbar;
                    title(Hc, 'Norm.');
                end
                set(handles.axes1, 'View', [0,90]);
            else
                set(handles.axes1, 'View', lastview);
            end
            
        
        % 2D plane
        elseif contains(handles.alldata.(['Data', num2str(i)]).datatype, '2D')
            % Disable mousescroll
            dimchar1 = '';
            dimchar2 = '';
            set(gcf, 'WindowButtonMotionFcn', []);
            set(gcf, 'WindowScrollWheelFcn', []);
            set(handles.quickprofile, 'Visible', 'On');
            % Set gamma minimum dose level edit-box and text visible
            set(handles.gammarefdoseEdit, 'Visible', 'On');
            set(handles.mingammadosetxt, 'Visible', 'On');
            set(handles.checkboxfixref, 'Visible', 'On');
            
            % Make 1D buttons visible
            set(handles.pushbuttonAddLine, 'Visible', 'On');
            set(handles.pushbuttonAddPlane, 'Visible', 'Off');
            % Reset chosen gamma analysis value
            if ismember('chosengam', fieldnames(handles.alldata.(['Data', num2str(i)])))
                handles.alldata.(['Data', num2str(i)]).chosengam = 0;
            end
            % Display available range
            set(handles.sizeText, 'Visible', 'On');
            set(handles.sizeText, 'String', ['Range [', handles.alldata.(['Data', num2str(i)]).Plane(1) ,',', handles.alldata.(['Data', num2str(i)]).Plane(2),']', '= [',...
                num2str(round(min(handles.alldata.(['Data', num2str(i)]).([lower(handles.alldata.(['Data', num2str(i)]).Plane(1)),'pos'])),2)), ':' , num2str(round(max(handles.alldata.(['Data', num2str(i)]).([lower(handles.alldata.(['Data', num2str(i)]).Plane(1)),'pos'])),2)),...
                ',', num2str(round(min(handles.alldata.(['Data', num2str(i)]).([lower(handles.alldata.(['Data', num2str(i)]).Plane(2)),'pos'])),2)), ':', num2str(round(max(handles.alldata.(['Data', num2str(i)]).([lower(handles.alldata.(['Data', num2str(i)]).Plane(2)),'pos'])),2)), ']']);
            
            % Show Gamma, DTA and DD buttons if existing in the data
           
            if sum(contains(fieldnames(handles.alldata.(['Data', num2str(i)])),'gamma')) > 0
                set(handles.show2dgamma, 'Visible', 'On');
                set(handles.show2dDTA, 'Visible', 'On');
                set(handles.show2dDD, 'Visible', 'On');
            end
            
            % If XY-plane data -> plot only the chosen matrix surface
            switch handles.residualcount
                case 0
                    % Don't use distance or mean normalization (not supported yet)
                    if ~contains(handles.normalizationtype, '%') && ~contains(lower(handles.normalizationtype), 'mean')
                        if handles.lowresvisualization == 1
                            plotthis = handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype]);
                        else
                            plotthis = handles.alldata.(['Data', num2str(i)]).(['data', handles.normalizationtype]);
                        end
                    else
                        if handles.lowresvisualization == 1
                            plotthis = handles.alldata.(['Data', num2str(i)]).DisplaydataNON;
                        else
                            plotthis = handles.alldata.(['Data', num2str(i)]).dataNON;
                        end
                        disp('Chosen normalization type does not exist for chosen data. Displaying non normalized data instead');
                    end
                case 1
                    handles.tempdata = handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype]);
                    plotthis = handles.tempdata;
                    handles.residualcount = handles.residualcount + 1;
                case 2
                    plotthis = tempdata  - handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype]);
                    handles.residualcount = 1;
            end
            if strcmp(handles.alldata.(['Data', num2str(i)]).Plane, 'XY')
                try
                    set(handles.textdepth, 'Visible', 'On', 'String', ['Z: ', handles.alldata.(['Data', num2str(i)]).depth]);
                catch
                    set(handles.textdepth, 'Visible', 'On', 'String', 'Z: 0');
                end
                
                if handles.lowresvisualization == 1
                    surf(handles.alldata.(['Data', num2str(i)]).Displayxpos, handles.alldata.(['Data', num2str(i)]).Displayypos,...
                    plotthis);
                    dimchar1 = 'Displayy';
                    dimchar2 = 'Displayx';
                else
                    surf(handles.alldata.(['Data', num2str(i)]).xpos, handles.alldata.(['Data', num2str(i)]).ypos,...
                    plotthis);
                    dimchar1 = 'y';
                    dimchar2 = 'x';
                end
                
                xlabel('Distance X [cm]');
                ylabel('Distance Y [cm]');
            elseif strcmp(handles.alldata.(['Data', num2str(i)]).Plane, 'XZ')
                if handles.lowresvisualization == 1
                    surf(handles.alldata.(['Data', num2str(i)]).Displayzpos, handles.alldata.(['Data', num2str(i)]).Displayxpos,...
                    plotthis);
                    dimchar1 = 'Displayx';
                    dimchar2 = 'Displayz';
                else
                    surf(handles.alldata.(['Data', num2str(i)]).zpos, handles.alldata.(['Data', num2str(i)]).xpos,...
                    plotthis);
                    dimchar1 = 'x';
                    dimchar2 = 'z';
                end
                
                xlabel('Distance Z [cm]');
                ylabel('Distance X [cm]');
            else
                if handles.lowresvisualization == 1
                    surf(handles.alldata.(['Data', num2str(i)]).Displayzpos, handles.alldata.(['Data', num2str(i)]).Displayypos,...
                    plotthis);
                    dimchar1 = 'Displayy';
                    dimchar2 = 'Displayz';
                else
                    surf(handles.alldata.(['Data', num2str(i)]).zpos, handles.alldata.(['Data', num2str(i)]).ypos,...
                    plotthis);  
                    dimchar1 = 'y';
                    dimchar2 = 'z';
                end
                
                xlabel('Distance Z [cm]');
                ylabel('Distance Y [cm]');
            end
            
            handles.chosendata2d{chosencount}{1} = plotthis;
            handles.chosendata2d{chosencount}{2} = handles.alldata.(['Data', num2str(i)]).([dimchar2, 'pos']);
            handles.chosendata2d{chosencount}{3} = handles.alldata.(['Data', num2str(i)]).([dimchar1, 'pos']);
            handles.chosendata{chosencount}(1,3) = i;
            
            chosencount = chosencount + 1;
            
            handles.listboxArray{i} = handles.origlist{i};
            handles.listboxArray{i} = ['<HTML><B><font face="Comic Sans MS">', handles.listboxArray{i}, '    **Chosen** <B>'];
            handles.legendlist{end+1} = handles.origlist{i};
            % Add colorbar and title as units
            try
                if strcmp(handles.normalizationtype, 'NON')
                    Hc = colorbar;
                    title(Hc, handles.alldata.(['Data', num2str(i)]).DataUnit);
                else
                    Hc = colorbar;
                    title(Hc, 'Norm.');
                end
            catch
                
            end
            
            if strcmp(handles.normalizationtype, 'NON')
                try
                    zlabel(['Value ', '[',  handles.alldata.(['Data', num2str(i)]).DataUnit, ']']);
                catch
                    zlabel(['Value ', '[Unknown unit]']);
                end
            else
                zlabel('Value [normalized]');
            end
            axis([-inf inf -inf inf]);
            
            
            
        elseif contains(handles.alldata.(['Data', num2str(i)]).datatype, 'Isodose')
            

            % Disable mousescroll
            dimchar1 = '';
            dimchar2 = '';
            set(gcf, 'WindowButtonMotionFcn', []);
            set(gcf, 'WindowScrollWheelFcn', []);
            % Set gamma minimum dose level edit-box and text visible
            set(handles.gammarefdoseEdit, 'Visible', 'On');
            set(handles.mingammadosetxt, 'Visible', 'on');
            
            % Make 1D buttons visible
            set(handles.pushbuttonAddLine, 'Visible', 'On');
            set(handles.pushbuttonAddPlane, 'Visible', 'Off');
            % Reset chosen gamma analysis value
            if ismember('chosengam', fieldnames(handles.alldata.(['Data', num2str(i)])))
                handles.alldata.(['Data', num2str(i)]).chosengam = 0;
            end
            
            % Show Gamma, DTA and DD buttons if existing in the data
           
            if sum(contains(fieldnames(handles.alldata.(['Data', num2str(i)])),'gamma')) > 0
                set(handles.show2dgamma, 'Visible', 'On');
                set(handles.show2dDTA, 'Visible', 'On');
                set(handles.show2dDD, 'Visible', 'On');
            end
            
            % If XY-plane data -> plot only the chosen matrix surface
            switch handles.residualcount
                case 0
                    % Don't use distance or mean normalization (not supported yet)
                    if ~contains(handles.normalizationtype, '%') && ~contains(lower(handles.normalizationtype), 'mean')
                        if handles.lowresvisualization == 1
                            plotthis = handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype]);
                        else
                            plotthis = handles.alldata.(['Data', num2str(i)]).(['data', handles.normalizationtype]);
                        end
                    else
                        if handles.lowresvisualization == 1
                            plotthis = handles.alldata.(['Data', num2str(i)]).DisplaydataNON;
                        else
                            plotthis = handles.alldata.(['Data', num2str(i)]).dataNON;
                        end
                        disp('Chosen normalization type does not exist for chosen data. Displaying non normalized data instead');
                    end
                case 1
                    handles.tempdata = handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype]);
                    plotthis = handles.tempdata;
                    handles.residualcount = handles.residualcount + 1;
                case 2
                    plotthis = tempdata  - handles.alldata.(['Data', num2str(i)]).(['Displaydata', handles.normalizationtype]);
                    handles.residualcount = 1;
            end
            if strcmp(handles.alldata.(['Data', num2str(i)]).Plane, 'XY')
                try
                    set(handles.textdepth, 'Visible', 'On', 'String', ['Z: ', handles.alldata.(['Data', num2str(i)]).depth]);
                catch
                    set(handles.textdepth, 'Visible', 'On', 'String', 'Z: 0');
                end
                
                if handles.lowresvisualization == 1
                    surf(handles.alldata.(['Data', num2str(i)]).Displayxpos, handles.alldata.(['Data', num2str(i)]).Displayypos,...
                    plotthis);
                    dimchar1 = 'Displayy';
                    dimchar2 = 'Displayx';
                else
                    surf(handles.alldata.(['Data', num2str(i)]).xpos, handles.alldata.(['Data', num2str(i)]).ypos,...
                    plotthis);
                    dimchar1 = 'y';
                    dimchar2 = 'x';
                end
                
                xlabel('Distance X [cm]');
                ylabel('Distance Y [cm]');
            elseif strcmp(handles.alldata.(['Data', num2str(i)]).Plane, 'XZ')
                if handles.lowresvisualization == 1
                    surf(handles.alldata.(['Data', num2str(i)]).Displayzpos, handles.alldata.(['Data', num2str(i)]).Displayxpos,...
                    plotthis);
                    dimchar1 = 'Displayx';
                    dimchar2 = 'Displayz';
                else
                    surf(handles.alldata.(['Data', num2str(i)]).zpos, handles.alldata.(['Data', num2str(i)]).xpos,...
                    plotthis);
                    dimchar1 = 'x';
                    dimchar2 = 'z';
                end
                
                xlabel('Distance Z [cm]');
                ylabel('Distance X [cm]');
            else
                if handles.lowresvisualization == 1
                    surf(handles.alldata.(['Data', num2str(i)]).Displayzpos, handles.alldata.(['Data', num2str(i)]).Displayypos,...
                    plotthis);
                    dimchar1 = 'Displayy';
                    dimchar2 = 'Displayz';
                else
                    surf(handles.alldata.(['Data', num2str(i)]).zpos, handles.alldata.(['Data', num2str(i)]).ypos,...
                    plotthis);  
                    dimchar1 = 'y';
                    dimchar2 = 'z';
                end
                
                xlabel('Distance Z [cm]');
                ylabel('Distance Y [cm]');
            end
            
            handles.chosendata2d{chosencount}{1} = plotthis;
            handles.chosendata2d{chosencount}{2} = handles.alldata.(['Data', num2str(i)]).([dimchar2, 'pos']);
            handles.chosendata2d{chosencount}{3} = handles.alldata.(['Data', num2str(i)]).([dimchar1, 'pos']);
            handles.chosendata{chosencount}(1,3) = i;
            
            chosencount = chosencount + 1;
            
            handles.listboxArray{i} = handles.origlist{i};
            handles.listboxArray{i} = ['<HTML><B><font face="Comic Sans MS">', handles.listboxArray{i}, '    **Chosen** <B>'];
            handles.legendlist{end+1} = handles.origlist{i};
            % Add colorbar and title as units
            try
                if strcmp(handles.normalizationtype, 'NON')
                    Hc = colorbar;
                    title(Hc, handles.alldata.(['Data', num2str(i)]).DataUnit);
                else
                    Hc = colorbar;
                    title(Hc, 'Norm.');
                end
            catch
                
            end
            
            if strcmp(handles.normalizationtype, 'NON')
                try
                    zlabel(['Value ', '[',  handles.alldata.(['Data', num2str(i)]).DataUnit, ']']);
                catch
                    zlabel(['Value ', '[Unknown unit]']);
                end
            else
                zlabel('Value [normalized]');
            end
            axis([-inf inf -inf inf]);
            
            
            
            
        % 1D line
        else
            % Disable mousescroll
            set(gcf, 'WindowButtonMotionFcn', []);
            set(gcf, 'WindowScrollWheelFcn', []);
            set(handles.quickprofile, 'Visible', 'Off');
            % Turn 1D buttons visibilty off
            set(handles.pushbuttonAddLine, 'Visible', 'Off');
            set(handles.pushbuttonAddPlane, 'Visible', 'Off');
            % Set gamma minimum dose edit and related text OFF as long as
            % not existing to 3D data
            set(handles.gammarefdoseEdit, 'Visible', 'Off');
            set(handles.mingammadosetxt, 'Visible', 'Off');
            set(handles.checkboxfixref, 'Visible', 'Off');
            
            try
                delete(colorbar); 
            catch
                disp('Colorbar error!');
            end
            try
                switch handles.residualcount
                    case 0
                        plotthis = handles.alldata.(['Data', num2str(i)]).(['data', handles.normalizationtype]);
                    case 1
                        handles.tempdata = handles.alldata.(['Data', num2str(i)]).(['data', handles.normalizationtype]);
                        plotthis = handles.tempdata;
                        handles.residualcount = handles.residualcount + 1;
                    case 2
                        plotthis = tempdata  - handles.alldata.(['Data', num2str(i)]).(['data', handles.normalizationtype]);
                        handles.residualcount = 1;
                end
                % Plot profiles, if chosen
                plotpos = handles.alldata.(['Data', num2str(i)]).pos;
                plot(plotpos, plotthis, 'Linewidth', 1.5);
                handles.legendlist{end+1} = handles.origlist{i};
  
                % Display error data if 3D->1D data error checkbox is chosen
                if handles.showError == 1 && (ismember('error', fieldnames(handles.alldata.(['Data', num2str(i)]))) || ismember('reference', fieldnames(handles.alldata.(['Data', num2str(i)]))))
                        if ismember('error', fieldnames(handles.alldata.(['Data', num2str(i)])))
                            errorup = plotthis + (handles.alldata.(['Data', num2str(i)]).error.*plotthis)/2;
                            errorlow = plotthis - (handles.alldata.(['Data', num2str(i)]).error.*plotthis)/2;
                            plot(plotpos, errorup, '--');
                            handles.legendlist{end+1} = 'Error+ (total error/2)';

                            plot(plotpos, errorlow, '--');
                            handles.legendlist{end+1} = 'Error- (total error/2)';
                        else
                            errorup = (handles.alldata.(['Data', num2str(i)]).reference);
                            plot(plotpos, errorup, '--');
                            handles.legendlist{end+1} = 'Reference (Raw/100)';
                        end
                end
                
                
                % Fill the parameters table (if parameters are existing)
                dtype = handles.alldata.(['Data', num2str(i)]).datatype;
                if ismember('params',fieldnames(handles.alldata.(['Data', num2str(i)]))) && ~isempty(fieldnames(handles.alldata.(['Data', num2str(i)])))
                    switch dtype
                        case 'profile'
                            if length(handles.origlist{i}) > 52
                                disp('Profile filename clipped to 52 characters');
                                profilerowdata{profilerowcount} = handles.origlist{i}(1:52);
                            else
                                profilerowdata{profilerowcount} = handles.origlist{i};
                            end
                            profiletabledata(profilerowcount,:) = handles.alldata.(['Data', num2str(i)]).params;
                            profilerowcount = profilerowcount+1;
                    
                        case 'pdd'
                            if length(handles.origlist{i}) > 64
                            disp('PDD filename clipped to 64 characters');
                            pddrowdata{pddrowcount} = handles.origlist{i}(1:64);
                            else
                                pddrowdata{pddrowcount} = handles.origlist{i};
                            end
                            pddtabledata(pddrowcount,:) = handles.alldata.(['Data', num2str(i)]).params;
                            pddrowcount = pddrowcount+1;
                    end
                else
                    disp('Data has no computed profile or pdd parameters');
                end
            catch
                % Display original data if current normalization data does
                % not exist
                plotthis = handles.alldata.(['Data', num2str(i)]).('dataNON');
                plotpos = handles.alldata.(['Data', num2str(i)]).pos;
                plot(plotpos, plotthis, 'Linewidth', 1.5);
                handles.legendlist{end+1} = handles.origlist{i};
                
                disp('Error in displaying data');
                disp(lasterror);
                 if ismember('params',fieldnames(handles.alldata.(['Data', num2str(i)])))
                     % Fill the parameters table (if parameters are existing)
                    dtype = handles.alldata.(['Data', num2str(i)]).datatype;
                    switch dtype
                        case 'profile'
                            if length(handles.origlist{i}) > 23
                                disp('Profile filename clipped to 23 characters');
                                profilerowdata{profilerowcount} = handles.origlist{i}(1:23);
                            else
                                profilerowdata{profilerowcount} = handles.origlist{i};
                            end
                            profiletabledata(profilerowcount,:) = handles.alldata.(['Data', num2str(i)]).params;
                            profilerowcount = profilerowcount+1;
                    
                        case 'pdd'
                            if length(handles.origlist{i}) > 23
                                disp('PDD filename clipped to 23 characters');
                                pddrowdata{pddrowcount} = handles.origlist{i}(1:23);
                            else
                                pddrowdata{pddrowcount} = handles.origlist{i};
                            end
                            pddtabledata(pddrowcount,:) = handles.alldata.(['Data', num2str(i)]).params;
                            pddrowcount = pddrowcount+1;
                    end
                else
                    disp('Data has no computed profile or pdd parameters');
                end
            end
            % Set view back to 1d view
            view(0,90);
            axis([-inf inf 0 inf]);
            
            xlabel('Distance [cm]');
            if strcmp(handles.normalizationtype, 'NON')
                try
                    ylabel(['Value ', '[',  handles.alldata.(['Data', num2str(i)]).DataUnit, ']']);
                catch
                    ylabel(['Value ', '[Unknown unit]']);
                end
            else
                ylabel('Value [normalized]');
            end
            handles.chosendata{chosencount}(:,1) = plotthis;
            handles.chosendata{chosencount}(:,2) = plotpos;
            handles.chosendata{chosencount}(1,3) = i;
            chosencount = chosencount + 1;
        end
        
    else
        handles.listboxArray{i} = handles.origlist{i};
    end
    
end


% Fill parameter tables
try
    set(handles.UItable, 'RowName', profilerowdata);
    set(handles.UItable, 'Data', profiletabledata);
    set(handles.UItablePDD, 'RowName', pddrowdata);
    set(handles.UItablePDD, 'Data', pddtabledata);
    
catch
    disp('No parameters computed'); 
end






legend(handles.legendlist, 'Interpreter', 'None', 'Location','northeast');

drawnow;

hold off;

set(handles.listboxChooseData, 'String', handles.listboxArray);
guidata(hObject, handles);


% Hints: contents = cellstr(get(hObject,'String')) returns listboxChooseData contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxChooseData


% --- Executes during object creation, after setting all properties.
function listboxChooseData_CreateFcn(hObject, ~, handles)
% hObject    handle to listboxChooseData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tempe = get(hObject,'String');

handles.customNomalization = str2double(tempe);

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuNormalisointi.
function popupmenuNormalisointi_Callback(hObject, ~, handles)
% hObject    handle to popupmenuNormalisointi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuNormalisointi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuNormalisointi

 contents = cellstr(get(hObject,'String'));
 Chosen =  contents{get(hObject,'Value')};

 
 switch Chosen
     
     case 'Distance/%'
         set(handles.edit2, 'Visible', 'On');
         set(handles.edit6, 'Visible', 'On');
         set(handles.meanmin, 'Visible', 'Off');
         set(handles.meanmax, 'Visible', 'Off');
         set(handles.text9, 'Visible', 'On');
         set(handles.text10, 'Visible', 'On');
         set(handles.editmeanperc, 'Visible', 'Off');
         set(handles.textmeanperc, 'Visible', 'Off');
         set(handles.text10, 'String', '%');
         % Get normalization values for manual normalization
         handles.normalizationdist = str2double(get(handles.edit2, 'String'));
         handles.normalizationperc = str2double(get(handles.edit6, 'String'));
         % Renormalize if has been already read
         % Omnipro
         try
            for normcount = 1:length(fieldnames(handles.alldata))
            if ~contains(handles.alldata.(['Data', num2str(normcount)]).name, "XY") % Ignore 2D data
                tempdata = [];
                if ~isempty(handles.alldata.(['Data', num2str(normcount)]).interpolated(:,1))
                    % Use interpolated data to find the new normalization
                    % constant
                    tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
                    % Find the normalization distance (/10 since position is in cm)
                    [~, idxdist1] = min(abs(tempdata(:,1) - handles.normalizationdist/10));
                    handles.alldata.(['Data', num2str(normcount)]).dataMAN = handles.normalizationperc/100.*handles.alldata.(['Data', num2str(normcount)]).dataNON/tempdata(idxdist1,2); % Note: Omnidata divided with 10 instead of 100, since interpolated data is multiplied with 10 before!
                else
                    % Use MAX normalized values for recalcomputation (avoids errors when PDD and profile data are mixed)
                    tempdata = handles.alldata.(['Data', num2str(normcount)]).dataMAX;
                    % Find the normalization distance (/10 since position is in cm)
                    [~, idxdist1] = min(abs(handles.alldata.(['Data', num2str(normcount)]).pos(:) - handles.normalizationdist/10));
                    handles.alldata.(['Data', num2str(normcount)]).dataMAN = handles.normalizationperc/100.*tempdata(:)/tempdata(idxdist1);
                end
            end
            end
         catch
            disp('Error in distance-% normalization! (Data may not exist) Normalization have not been changed');
         end
         
         
         % Change the normalization type keyword
         handles.normalizationtype = 'MAN';
     case 'Max 100%'
         set(handles.edit2, 'Visible', 'Off');
         set(handles.edit6, 'Visible', 'Off');
         set(handles.modifytxt, 'Visible', 'Off');
         set(handles.text9, 'Visible', 'Off');
         set(handles.text10, 'Visible', 'Off');
         set(handles.meanmin, 'Visible', 'Off');
         set(handles.meanmax, 'Visible', 'Off');
         set(handles.editmeanperc, 'Visible', 'Off');
         set(handles.textmeanperc, 'Visible', 'Off');
         handles.normalizationtype = 'MAX';
     case 'CAX 100%'
         set(handles.edit2, 'Visible', 'Off');
         set(handles.edit6, 'Visible', 'Off');
         set(handles.modifytxt, 'Visible', 'Off');
         set(handles.text9, 'Visible', 'Off');
         set(handles.text10, 'Visible', 'Off');
         set(handles.meanmin, 'Visible', 'Off');
         set(handles.meanmax, 'Visible', 'Off');
         set(handles.editmeanperc, 'Visible', 'Off');
         set(handles.textmeanperc, 'Visible', 'Off');
         handles.normalization = get(handles.edit2, 'Value');
         handles.normalizationtype = 'CAX';
     case 'No normalization'
         set(handles.edit2, 'Visible', 'Off');
         set(handles.modifytxt, 'Visible', 'Off');
         set(handles.edit6, 'Visible', 'Off');
         set(handles.text9, 'Visible', 'Off');
         set(handles.text10, 'Visible', 'Off');
         set(handles.meanmin, 'Visible', 'Off');
         set(handles.meanmax, 'Visible', 'Off');
         set(handles.editmeanperc, 'Visible', 'Off');
         set(handles.textmeanperc, 'Visible', 'Off');
         handles.normalization = get(handles.edit2, 'Value');
         handles.normalizationtype = 'NON';
     case 'Normalize to value'
         set(handles.edit2, 'Visible', 'On');
         set(handles.modifytxt, 'Visible', 'On');
         try
            set(handles.modifytxt, 'String', 'Dataunit');
         catch
            set(handles.modifytxt, 'String', 'Unknown unit');
         end
         set(handles.edit6, 'Visible', 'Off');
         set(handles.text9, 'Visible', 'Off');
         set(handles.modifytxt, 'Visible', 'Off');
         set(handles.text10, 'Visible', 'Off');
         set(handles.meanmin, 'Visible', 'Off');
         set(handles.meanmax, 'Visible', 'Off');
         set(handles.editmeanperc, 'Visible', 'Off');
         set(handles.textmeanperc, 'Visible', 'Off');
         handles.normalization = get(handles.edit2, 'Value');
         handles.normalizationtype = 'TOVAL';
         guidata(hObject, handles);
         %edit2_Callback(hObject, [], handles);
         %handles = guidata(hObject);
         
     case 'Mean value from range'
         set(handles.edit2, 'Visible', 'Off');
         set(handles.edit6, 'Visible', 'Off');
         set(handles.modifytxt, 'Visible', 'Off');
         set(handles.text9, 'Visible', 'On');
         set(handles.text10, 'Visible', 'On');
         set(handles.text10, 'String', 'cm');
         set(handles.meanmin, 'Visible', 'On');
         set(handles.meanmax, 'Visible', 'On');
         set(handles.editmeanperc, 'Visible', 'On');
         set(handles.textmeanperc, 'Visible', 'On');
         handles.normalizationtype = 'MEAN';
        
        handles.meanperc = str2double(get(handles.editmeanperc, 'String'))/100;
        % Try to fetch the minimum and maximum position values from edit
        try
            handles.minmeanval = str2double(get(handles.meanmin, 'Value'))/10;
        catch
           disp('Error! Minimum position have not been set: Dataset minimum chosen'); 
        end
        try
            handles.maxmeanval = str2double(get(handles.meanmax, 'Value'))/10;
        catch
           disp('Error! Maximum position have not been set: Dataset maximum chosen'); 
        end
        
        
        % See if the values are set for data
        if isempty(handles.meanminval)
           handles.meanminval = min(handles.alldata.Data1.pos);
           set(handles.meanmin, 'String', num2str(handles.meanminval));
        end
        if isempty(handles.meanmaxval)
           handles.meanmaxval = max(handles.alldata.Data1.pos);
           set(handles.meanmax, 'String', num2str(handles.meanmaxval));
        end
        
        % Data
        for normcount = 1:length(fieldnames(handles.alldata))
            if ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "2D") && ...
                    ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "3D") % Only done for 1D. Ignore 2D and 3D data
                tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
                % Find the closest distance-value for min and max from the data
                [~, idxmindist] = min(abs(tempdata(:,1) - handles.meanminval));
                [~, idxmaxdist] = min(abs(tempdata(:,1) - handles.meanmaxval));
                tempnormval = mean(tempdata(idxmindist:idxmaxdist,2));
                handles.alldata.(['Data', num2str(normcount)]).dataMEAN = handles.alldata.(['Data', num2str(normcount)]).dataNON/tempnormval;
            end
        end   
 end
 
 % Update the figure with only re-drawing chosen data
 handles.redraw = 1;
 guidata(hObject, handles);
try
    % Display modified data
    listboxChooseData_Callback(hObject, [], handles);
catch
    disp('Error in displaying data after changing normalization (Data may not exist)');
end
 

% --- Executes during object creation, after setting all properties.
function popupmenuNormalisointi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuNormalisointi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if strcmp(handles.normalizationtype,'TOVAL')
    handles.normalizationtoval = str2double(get(handles.edit2,'String'));
    for normcount = 1:length(fieldnames(handles.alldata))
        if strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, '3D')
            handles.alldata.(['Data', num2str(normcount)]).dataTOVAL = handles.alldata.(['Data', num2str(normcount)]).data/handles.normalizationtoval;
            handles.alldata.(['Data', num2str(normcount)]).DisplaydataTOVAL = handles.alldata.(['Data', num2str(normcount)]).Displaydata/handles.normalizationtoval;
        elseif strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, '2D')
            handles.alldata.(['Data', num2str(normcount)]).dataTOVAL = handles.alldata.(['Data', num2str(normcount)]).dataNON/handles.normalizationtoval;
            handles.alldata.(['Data', num2str(normcount)]).DisplaydataTOVAL = handles.alldata.(['Data', num2str(normcount)]).DisplaydataNON/handles.normalizationtoval;
        else
            handles.alldata.(['Data', num2str(normcount)]).dataTOVAL = handles.alldata.(['Data', num2str(normcount)]).dataNON/handles.normalizationtoval;
        end
    end
else
    

    handles.normalizationdist = str2double(get(hObject,'String'));
    handles.normalizationperc = str2double(get(handles.edit6, 'String'));
    % Renormalize if has been already read

    for normcount = 1:length(fieldnames(handles.alldata))
    if ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "2D") && ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "3D") % Ignore 2D and 3D data
        tempdata = [];
        if ~isempty(handles.alldata.(['Data', num2str(normcount)]).interpolated(:,1))
            % Use interpolated data to find the new normalization
            % constant
            tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
            % Find the normalization distance (/10 since position is in cm)
            [~, idxdist1] = min(abs(tempdata(:,1) - handles.normalizationdist));
            handles.alldata.(['Data', num2str(normcount)]).dataMAN = handles.normalizationperc/100.*handles.alldata.(['Data', num2str(normcount)]).dataNON/tempdata(idxdist1,2); % Note: Omnidata divided with 10 instead of 100, since interpolated data is multiplied with 10 before!
        else
            % Use MAX normalized values for recalcomputation (avoids errors when PDD and profile data are mixed)
            tempdata = handles.alldata.(['Data', num2str(normcount)]).dataMAX;
            % Find the normalization distance (/10 since position is in cm)
            [~, idxdist1] = min(abs(handles.alldata.(['Data', num2str(normcount)]).pos(:) - handles.normalizationdist));
            handles.alldata.(['Data', num2str(normcount)]).dataMAN = handles.normalizationperc/100.*tempdata(:)/tempdata(idxdist1);
        end
    end
    end

end
% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);
try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in displaying data after edit2 (Data may not exist)'); 
end


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


% --- Executes on button press in residualButton.
function residualButton_Callback(hObject, eventdata, handles)
% hObject    handle to residualButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.residualState = get(hObject, 'Value');

switch handles.residualState
    case 1
        handles.listboxArray = handles.origlist;
        handles.idx(:) = 0;
        set(handles.textconsole, 'String', 'Valitse 2 tiedostoa (samat dimensiot)');
        handles.redidualcount = 1;
    case 0
        handles.residualcount = 0;
end

guidata(hObject, handles);


% --- Executes on button press in Gammabutton.
function Gammabutton_Callback(hObject, eventdata, handles)

% Choose 2 files from the listbox and AFTER that press gamma-analysis
% First one chosen is the reference data, 2nd is the target data

% Data format for CalcGamma: Structure with fields .start, .width, data

% Swap data to right order according to which was chosen first
compare = find(handles.order);
dataord = 'orig';
if handles.order(compare(1)) > handles.order(compare(2))
   temp =  handles.order(compare(1));
   handles.order(compare(1)) = handles.order(compare(2));
   handles.order(compare(2)) = temp;
   tempdata = handles.chosendata{1};
   handles.chosendata{1} = handles.chosendata{2};
   handles.chosendata{2} = tempdata;
   disp('Data1 = ref, data2 = target, computing Gamma parameters...');
   dataord = 'reverse';
end



% 1D case
if ~strcmp(handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).datatype, '2D') && ~strcmp(handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).datatype, '3D')
    % Initialization, gamma calculation requires same sized input data
    interpres = 0.001;
    refdata = struct;
    targetdata = struct;
    refdata.data = handles.chosendata{1}(:,1);
    refdata.pos = handles.chosendata{1}(:,2);
    targetdata.data = handles.chosendata{2}(:,1);
    targetdata.pos = handles.chosendata{2}(:,2);

    % Sample intervals at CAX (use if interval not chosen)
    refinterval = abs(refdata.pos(floor(end/2)) - refdata.pos(floor(end/2)+1));
    targinterval = abs(targetdata.pos(floor(end/2)) - targetdata.pos(floor(end/2)+1));
    

    % Find the start position defined by the smaller data
    [startposition, idxRT] = max([min(refdata.pos), min(targetdata.pos)]);

    if ~isempty(handles.interpres)
        % Use manually set interval (from interpolation edit box)
        interval = handles.interpres/10;
    else
        switch idxRT
            case 1
                % Use reference cax sampling as interval
                interval = refinterval;
                disp(['Gamma interpolation interval set to reference data resolution: ', num2str(refinterval), 'cm']);
            case 2
                % Use target cax sampling as interval
                interval = targinterval;
                disp(['Gamma interpolation interval set to target data resolution: ', num2str(targinterval), 'cm']);
            otherwise
                % Set interval to 1mm
                interval = 0.01;
                disp('Interpolation interval set to 1mm');
        end
    end

    temp = [max(refdata.pos), max(targetdata.pos)];
            endposition = temp(idxRT);


    % Interpolate to same samplerate from same query points
    %if resampletarg == 1
        disp('Resampling target data...');
        targetdata.data = interp1(targetdata.pos, targetdata.data, startposition:interval:endposition);
        targetdata.pos = startposition:interval:endposition;
    %end
    %if resampleref == 1
        disp('Resampling reference data...');
        refdata.data = interp1(refdata.pos, refdata.data, startposition:interval:endposition);
        refdata.pos = startposition:interval:endposition;
    %end



    % Reference data (Chosen 1st)
    refdata.start = startposition;                                                            % Ref start position
    refdata.width = interval;                                                                 % Ref grid size

    % Target data (chosen 2nd)
    targetdata.start = startposition;                                                         % Target start position
    targetdata.width = interval;                                                              % Target grid size


    % Compute DTA and dose difference
    [dta, relativeDD, absoluteDD, dtaA, dtapassperc, dtaApassperc, DDpassperc] = DtaDComp(refdata, targetdata, handles.dosecriterion, handles.dtacriterion,...
        interpres, interval, [handles.DTAanalytical, handles.DTAbruteforce]);

    %dta = dta/10;
    % Compute gamma using function from: Author: Mark Geurts, mark.w.geurts@gmail.com
    % Comment this if you want to save time and only compute DTA and DD
    if handles.nogammacalc == 0
        gamma2 = CalcGamma(refdata, targetdata, handles.dosecriterion, handles.dtacriterion/10, 'local', 1);
        [passtemp,~] = find(gamma2 <= 1);
        gammapassperc = sum(passtemp)/length(gamma2)*100;
    end
    %save('gamma2.mat', 'gamma2');

    % Which was ref and which targ?:
    [~, exidx] = find(handles.idx);
    if strcmp(dataord,'reverse')
        exidx = flip(exidx);
    end

    % exportables = handles.listboxArray(exidx);
    % exporttype = contains(exportables, 'Matrix');
     targetname = handles.origlist{exidx(2)};
     referencename = handles.origlist{exidx(1)};

    % Add Gammadata only to the reference structure
    %handles.chosendata{1}(1,3)
    % Add dta and dose difference data to the reference structure
    if handles.DTAbruteforce == 1
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_',num2str(exidx(1)), '_', num2str(exidx(2))]).dta = dta;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).dtapassperc = dtapassperc;
    end
    if handles.DTAanalytical == 1
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_',num2str(exidx(1)), '_', num2str(exidx(2))]).dtaA = dtaA;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).dtaApassperc = dtaApassperc;
    end
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_' ,num2str(exidx(1)), '_', num2str(exidx(2))]).relativeDD = relativeDD;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).absoluteDD = absoluteDD;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).DDpassperc = DDpassperc;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).target = targetname;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).reference = referencename;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).GAMpos = refdata.pos;
    if handles.nogammacalc == 0
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).gamma = gamma2;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).gammapassperc = gammapassperc;
    end
    % Plot data
    hold on;
    grid on;
    if handles.DTAbruteforce == 1
        yyaxis right
        plot(refdata.pos, dta, 'LineStyle', '-', 'LineWidth', 1.2, 'color', 'b');
    end
    if handles.DTAanalytical == 1
        yyaxis right
        plot(refdata.pos, dtaA, 'LineStyle', '-', 'LineWidth', 1.2, 'color', 'c');
    end
    yyaxis right
    plot(refdata.pos, relativeDD, 'LineStyle', '-', 'LineWidth', 1.2, 'color', 'r');
    plot(refdata.pos, absoluteDD, 'LineStyle', '-', 'LineWidth', 1.2, 'color', 'g');
    if handles.nogammacalc == 0
        yyaxis right
        plot(refdata.pos, gamma2, 'LineStyle', '-','LineWidth', 1.2, 'color', 'k')

        %save('gamma.mat', 'gamma2');
    end
    axis([refdata.pos(1) refdata.pos(end) -inf inf]);
    
    % Add line names to legendlist
    if handles.DTAbruteforce == 1
        passchar = num2str(round(dtapassperc,2));
        % Make sure that displayed value has 1 decimal
        while length(passchar) < 4
            passchar(end+1) = ' ';
        end
        handles.legendlist{end+1} = ['DTA BF [cm] (pass: ', passchar(1:4), '%)'];
    end
    if handles.DTAanalytical == 1
        passchar = num2str(round(dtaApassperc,2));
        % Make sure that displayed value has 1 decimal
        while length(passchar) < 4
            passchar(end+1) = ' ';
        end
        handles.legendlist{end+1} = ['DTA AN [cm] (pass: ', passchar(1:4), '%)'];
    end

    passchar = num2str(round(DDpassperc,2));
    % Make sure that displayed value has 1 decimal
    while length(passchar) < 4
        passchar(end+1) = ' ';
    end
    handles.legendlist{end+1} = ['Ddiff rel [%] (pass: ', passchar(1:4), '%)'];
    handles.legendlist{end+1} = 'Ddiff abs [%]';
    if handles.nogammacalc == 0
        passchar = num2str(round(gammapassperc,2));
        % Make sure that displayed value has 1 decimal
        while length(passchar) < 4
            passchar(end+1) = ' ';
        end
        handles.legendlist{end+1} = ['Gamma (pass: ', passchar(1:4), '%)'];
    end
    legend(handles.legendlist, 'Interpreter', 'None');

    hold off;

    % 2D case
elseif strcmp(handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).datatype, '2D')
    ff = msgbox('Wait! Window will close automatically');
    % interpolate both data to same grid
    interpres = 0.001;
    refdata = struct;
    targetdata = struct;
    refdata.data = handles.chosendata2d{1}{1};
    refdata.posx = handles.chosendata2d{1}{2};
    refdata.posy = handles.chosendata2d{1}{3};
    targetdata.data = handles.chosendata2d{2}{1};
    targetdata.posx = handles.chosendata2d{2}{2};
    targetdata.posy = handles.chosendata2d{2}{3};
    
    
    % Find smallest dimensions and interpolate both data according to those
    % dimensions
    xstart = max([min(refdata.posx), min(targetdata.posx)]);
    ystart = max([min(refdata.posy), min(targetdata.posy)]);
    xend = min([max(refdata.posx), max(targetdata.posx)]);
    yend = min([max(refdata.posy), max(targetdata.posy)]);

    if ~isempty(handles.interpres)
        % Use manually set interval (from interpolation edit box)
        xstepsize = handles.interpres/10;
        ystepsize = handles.interpres/10;
    else
        % Smallest step size
        xstepsize = min([abs(targetdata.posx(2) - targetdata.posx(1)), abs(refdata.posx(2) - refdata.posx(1))]);
        ystepsize = min([abs(targetdata.posy(2) - targetdata.posy(1)), abs(refdata.posy(2) - refdata.posy(1))]);
    end
    
    gridx = xstart:xstepsize:xend;
    gridy = ystart:ystepsize:yend;
    
    % Interpolate
    [Xq,Yq] = meshgrid(gridx,gridy);
    [X,Y] = meshgrid(targetdata.posx,targetdata.posy);
    targetdata.data = interp2(X,Y,targetdata.data, Xq, Yq);
    targetdata.xpos = gridx;
    targetdata.ypos = gridy;
    
    [Xq,Yq] = meshgrid(gridx,gridy);
    [X,Y] = meshgrid(refdata.posx,refdata.posy);
    refdata.data = interp2(X,Y,refdata.data, Xq, Yq);
    refdata.xpos = gridx;
    refdata.ypos = gridy;
    
    targetdata.start = [gridy(1),gridx(1)];
    targetdata.width = [ystepsize,xstepsize];
    refdata.start = [gridy(1),gridx(1)];
    refdata.width = [ystepsize,xstepsize];

    % Find which values will be computed using the minimum dose threshold
    gammamask = ones(length(gridy), length(gridx));
    if ~isempty(handles.gammaMinimumdose)
        %
        handles.gammaMinimumdose
        [idxy, idxx] = find(refdata.data < handles.gammaMinimumdose);
        index = sub2ind([size(gammamask,1), size(gammamask,2)], idxy, idxx);
        gammamask(index) = 0;
    end
    
    refsave = refdata.data;
    
    [idxy, idxx] = find(gammamask);
    index = sub2ind([size(gammamask,1), size(gammamask,2)], idxy, idxx);
    
    % Compute gamma
    if handles.nogammacalc == 0
        testgamma2 = CalcGamma(refdata, targetdata, handles.dosecriterion, handles.dtacriterion/10, 'local', 1, 'limit', 2);
        
        [gidxydel, gidxxdel] = find(gammamask == 0);
        indexdel = sub2ind([size(gammamask,1), size(gammamask,2)], gidxydel, gidxxdel);
        
        passtemp = find(testgamma2(index) <= 1);
        gammapassperc = round(length(passtemp)/length(idxy)*100,2);
        testgamma2(indexdel) = NaN;
    end
    
    
    % Compute DTA and dose difference
    [dta, relativeDD, absoluteDD, dtaA, dtapassperc, dtaApassperc, DDpassperc] = DtaDComp(refdata, targetdata, handles.dosecriterion, handles.dtacriterion,...
        interpres, refdata.width, [handles.DTAanalytical, handles.DTAbruteforce], gammamask);
    

    % Which was ref and which targ?:
    [~, exidx] = find(handles.idx);
    if strcmp(dataord,'reverse')
        exidx = flip(exidx);
    end

    % exportables = handles.listboxArray(exidx);
    % exporttype = contains(exportables, 'Matrix');
     targetname = handles.origlist{exidx(2)};
     referencename = handles.origlist{exidx(1)};
    
    try 
    if handles.DTAbruteforce == 1
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_',num2str(exidx(1)), '_', num2str(exidx(2))]).dta = dta;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).dtapassperc = dtapassperc;
    end
    if handles.DTAanalytical == 1
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_',num2str(exidx(1)), '_', num2str(exidx(2))]).dtaA = dtaA;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).dtaApassperc = dtaApassperc;
    end
    catch
        
    end
    try
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_' ,num2str(exidx(1)), '_', num2str(exidx(2))]).relativeDD = relativeDD;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).absoluteDD = absoluteDD;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).DDpassperc = DDpassperc;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).reference = referencename;
        
    catch
    end
    if handles.nogammacalc == 0
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).gamma = testgamma2;
        handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).gammapassperc = gammapassperc;
    end 
    
    
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).GAMposx = [];
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).GAMposy = [];
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).GAMposx(:,1) = gridx;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).GAMposy(:,1) = gridy;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).(['gamma_', num2str(exidx(1)), '_', num2str(exidx(2))]).gammastruct = gammamask;
    handles.alldata.(['Data', num2str(handles.chosendata{1}(1,3))]).chosengam = 0;
    close(ff);
    set(handles.show2dgamma, 'Visible', 'On');
    set(handles.show2dDTA, 'Visible', 'On');
    set(handles.show2dDD, 'Visible', 'On');
end

%asd = handles.alldata;
%save('alladata.mat', 'asd');
% Update handles

guidata(hObject, handles);



function dosekrit_Callback(hObject, eventdata, handles)
% hObject    handle to dosekrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dosecriterion = str2double(get(hObject,'String'));

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of dosekrit as text
%        str2double(get(hObject,'String')) returns contents of dosekrit as a double


% --- Executes during object creation, after setting all properties.
function dosekrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dosekrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTAkrit_Callback(hObject, eventdata, handles)
% hObject    handle to DTAkrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dtacriterion = str2double(get(hObject,'String'));

guidata(hObject, handles);



% Hints: get(hObject,'String') returns contents of DTAkrit as text
%        str2double(get(hObject,'String')) returns contents of DTAkrit as a double


% --- Executes during object creation, after setting all properties.
function DTAkrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTAkrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function interpresolution_Callback(hObject, eventdata, handles)
% hObject    handle to interpresolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.interpres = str2double(get(hObject,'String'));

guidata(hObject, handles);


% Hints: get(hObject,'String') returns contents of interpresolution as text
%        str2double(get(hObject,'String')) returns contents of interpresolution as a double


% --- Executes during object creation, after setting all properties.
function interpresolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interpresolution (see GCBO)
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

handles.normalizationperc = str2double(get(hObject,'String'));
handles.normalizationdist = str2double(get(handles.edit2, 'String'));
% Renormalize if has been already read
% PTW data
% Renormalize if has been already read

for normcount = 1:length(fieldnames(handles.alldata))
if ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "2D") && ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "3D") % Ignore 2D and 3D data
    tempdata = [];
    if ~isempty(handles.alldata.(['Data', num2str(normcount)]).interpolated(:,1))
        % Use interpolated data to find the new normalization
        % constant
        tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
        % Find the normalization distance (/10 since position is in cm)
        [~, idxdist1] = min(abs(tempdata(:,1) - handles.normalizationdist));
        handles.alldata.(['Data', num2str(normcount)]).dataMAN = handles.normalizationperc/100.*handles.alldata.(['Data', num2str(normcount)]).dataNON/tempdata(idxdist1,2); % Note: Omnidata divided with 10 instead of 100, since interpolated data is multiplied with 10 before!
    else
        % Use MAX normalized values for recalcomputation (avoids errors when PDD and profile data are mixed)
        tempdata = handles.alldata.(['Data', num2str(normcount)]).dataMAX;
        % Find the normalization distance (/10 since position is in cm)
        [~, idxdist1] = min(abs(handles.alldata.(['Data', num2str(normcount)]).pos(:) - handles.normalizationdist));
        handles.alldata.(['Data', num2str(normcount)]).dataMAN = handles.normalizationperc/100.*tempdata(:)/tempdata(idxdist1);
    end
end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;
         
guidata(hObject, handles);

try
  listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Valitse data'); 
end


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


% --- Executes on button press in exportbutton.
function exportbutton_Callback(hObject, eventdata, handles)

% Check which data has been chosen from the listbox (contains("valittu") or by handles.idx)
% Export different parameters to different sheets (normalizations, dta, dosediff, ...)


% Column indexing for excel (A,B,C,...,AA,AB,AC...ZZ)
csvcols = char(65:90);
csvcols = num2cell(csvcols(1:length(csvcols)));
csvappend = csvcols;
for i = 1:length(csvcols)
   for j = 1:length(csvcols)
       csvappend{end+1} = [csvcols{i},csvcols{j}];
   end
end

% Collect data to array and then transfer to table

% Collect different normalizations, DTA, Gamma, DD to separate sheets
% CAX, MAX and no normalization + interpolated data are always computed (found from alldata struct)
% MEAN and MAN normalization computed if the user have chosen them before
% DTA, Gamma and DD are computed if user have performed the analysis

% Get filename
handles.myfile = get(handles.editexport, 'String');
if handles.exportcount > 0
   handles.myfile = [handles.myfile, '_', num2str(handles.exportcount), '.xls'];
end
myfile = handles.myfile;

% Loop through PTW and omnidata separately (simplest). Check which data
% exists -> create field names and copy the data
wcCONST = 1;
wcMAN = 1;
wcMEAN = 1;
wcINTRP = 1;
wcDTA = 1;
wcGAMMA = 1;
skipall = 0;

% ... Check if PTW data is empty and folder has been chosen...

try
    if ~isempty(handles.alldata) || ~isempty(get(handles.PTWfolder,'String'))
        alldata = handles.alldata;
        lenall = length(fieldnames(alldata));
    else
        lenall = 0;
        disp('No data to be exported')
    end
catch
    skipall = 1;
    lenall = 0;
    disp('Export error')
end

warning('off')

% Collect and export the data. Single column at a time, since the columns
% have different sizes.
if skipall == 0
% Waitbar
h = waitbar(0,'Exporting data...');
for i = 1:lenall
    % Export only data that has been chosen from listbox
    waitbar(i/length(fieldnames(alldata)),h)
    if contains(lower(handles.listboxArray{i}), 'chosen')

        % CAX, MAX and no normalization (always exist). Write data to 3rd row
        % and columnnames to 2nd

        % Measurement name
        measname = handles.origlist(i);
        xlswrite(myfile,measname, 'CAX normalization', [csvappend{wcCONST},'1']);
        xlswrite(myfile,measname, 'MAX normalization', [csvappend{wcCONST},'1']);
        xlswrite(myfile,measname, 'Not normalized', [csvappend{wcCONST},'1']);
        % Varnames
        varnames = {'Pos', 'Data'};
        xlswrite(myfile,varnames, 'CAX normalization', [[csvappend{wcCONST},'2'],':',[csvappend{wcCONST+1},'2']]);
        xlswrite(myfile,varnames, 'MAX normalization', [[csvappend{wcCONST},'2'],':',[csvappend{wcCONST+1},'2']]);
        xlswrite(myfile,varnames, 'Not normalized', [[csvappend{wcCONST},'2'],':',[csvappend{wcCONST+1},'2']]);
        % Position
        xlswrite(myfile,alldata.(['Data', num2str(i)]).pos, 'CAX normalization', [csvappend{wcCONST},'3']);
        xlswrite(myfile,alldata.(['Data', num2str(i)]).pos, 'MAX normalization', [csvappend{wcCONST},'3']);
        xlswrite(myfile,alldata.(['Data', num2str(i)]).pos, 'Not normalized', [csvappend{wcCONST},'3']);
        % Data
        xlswrite(myfile,alldata.(['Data', num2str(i)]).dataCAX, 'CAX normalization', [csvappend{wcCONST+1},'3']);
        xlswrite(myfile,alldata.(['Data', num2str(i)]).dataMAX, 'MAX normalization', [csvappend{wcCONST+1},'3']);
        xlswrite(myfile,alldata.(['Data', num2str(i)]).dataNON, 'Not normalized', [csvappend{wcCONST+1},'3']);
        % error and reference data for 3ddose and mcc
        % error 
        if isfield(alldata.(['Data', num2str(i)]), 'error')
            xlswrite(myfile,{'Error'}, 'Not normalized', [csvappend{wcCONST+2},'2']);
            xlswrite(myfile,alldata.(['Data', num2str(i)]).error, 'Not normalized', [csvappend{wcCONST+2},'3']);
        end
        % reference
        if isfield(alldata.(['Data', num2str(i)]), 'reference')
            xlswrite(myfile,{'Reference'}, 'Not normalized', [csvappend{wcCONST+2},'2']);
            xlswrite(myfile,alldata.(['Data', num2str(i)]).reference, 'Not normalized', [csvappend{wcCONST+2},'3']);
        end
        % Space between data
        wcCONST = wcCONST + 4;

        % Try manual data
        try
           xlswrite(myfile,measname, 'Manual normalization', [csvappend{wcMAN},'1']);
           xlswrite(myfile,varnames, 'Manual normalization', [[csvappend{wcMAN},'2'],':',[csvappend{wcMAN+1},'2']]);
           xlswrite(myfile,alldata.(['Data', num2str(i)]).pos, 'Manual normalization', [csvappend{wcMAN},'3']);
           xlswrite(myfile,alldata.(['Data', num2str(i)]).dataMAN, 'Manual normalization', [csvappend{wcMAN+1},'3']);
           % Space between data
           wcMAN = wcMAN + 3;
        catch
           disp(['No manually normalized dat to file:', alldata.(['Data', num2str(i)]).name]); 
        end
        % Try mean data
        try
           xlswrite(myfile,measname, 'Mean normalization', [csvappend{wcMEAN},'1']);
           xlswrite(myfile,varnames, 'Mean normalization', [[csvappend{wcMEAN},'2'],':',[csvappend{wcMEAN+1},'2']]);
           xlswrite(myfile,alldata.(['Data', num2str(i)]).pos, 'Mean normalization', [csvappend{wcMEAN},'3']);
           xlswrite(myfile,alldata.(['Data', num2str(i)]).dataMEAN, 'Mean normalization', [csvappend{wcMEAN+1},'3']);
           % Space between data
           wcMEAN = wcMEAN + 3;
        catch
           disp(['No manually normalized data:', alldata.(['Data', num2str(i)]).name]); 
        end

        % Try Gamma, DTA and DD
        gidx = find(contains(fieldnames(alldata.(['Data', num2str(i)])), 'gamma'));
        if ~isempty(gidx)
            for yc = 1:length(gidx)
                fields = fieldnames(alldata.(['Data', num2str(i)]));
                gammaname = fields{gidx(yc)};
                varnames = {'Pos'};
                gammafnames = fieldnames(alldata.(['Data', num2str(i)]).(gammaname));
                % Check which analyses exist
                if ismember('dtaA',gammafnames)
                    varnames{end+1} = ['dtaA(', num2str(round(alldata.(['Data', num2str(i)]).(gammaname).dtaApassperc),3),'%)'];
                end
                if ismember('dta',gammafnames)
                    varnames{end+1} = ['dtaB(', num2str(round(alldata.(['Data', num2str(i)]).(gammaname).dtapassperc),3),'%)'];
                end
                if ismember('absoluteDD',gammafnames)
                    varnames{end+1} = 'absDD';
                end
                if ismember('relativeDD',gammafnames)
                    varnames{end+1} = ['rDD (', num2str(round(alldata.(['Data', num2str(i)]).(gammaname).DDpassperc),3),'%)'];
                end
                if ismember('gamma', gammafnames)
                   varnames{end+1} = ['gamma (', num2str(round(alldata.(['Data', num2str(i)]).(gammaname).gammapassperc),3),'%)'];
                end
                
                
                 % File name
                xlswrite(myfile,{['Ref: ', char(measname)]}, 'Gammadata', [csvappend{wcDTA},'1']);
                % Target + ref
                xlswrite(myfile, {['Target: ', alldata.(['Data', num2str(i)]).(gammaname).target]}, 'Gammadata', [csvappend{wcDTA+1},'1']);
                % Variable names
                xlswrite(myfile,varnames, 'Gammadata', [csvappend{wcDTA},'2']);
                % Pos
                xlswrite(myfile,alldata.(['Data', num2str(i)]).(gammaname).GAMpos', 'Gammadata', [csvappend{wcDTA},'3']);
                
                
                % Data
                columnincrement = 1;
                if ismember('dtaA',gammafnames)
                    xlswrite(myfile,alldata.(['Data', num2str(i)]).(gammaname).dtaA', 'Gammadata', [csvappend{wcDTA+columnincrement},'3']);
                    columnincrement = columnincrement + 1;
                end
                if ismember('dta',gammafnames)
                    xlswrite(myfile,alldata.(['Data', num2str(i)]).(gammaname).dta, 'Gammadata', [csvappend{wcDTA+columnincrement},'3']);
                    columnincrement = columnincrement + 1;
                end
                if ismember('absoluteDD',gammafnames)
                    xlswrite(myfile,alldata.(['Data', num2str(i)]).(gammaname).absoluteDD', 'Gammadata', [csvappend{wcDTA+columnincrement},'3']);
                    columnincrement = columnincrement + 1;
                end
                if ismember('relativeDD',gammafnames)
                    xlswrite(myfile,alldata.(['Data', num2str(i)]).(gammaname).relativeDD', 'Gammadata', [csvappend{wcDTA+columnincrement},'3']);
                    columnincrement = columnincrement + 1;
                end
                if ismember('gamma', gammafnames)
                   xlswrite(myfile,alldata.(['Data', num2str(i)]).(gammaname).gamma', 'Gammadata', [csvappend{wcDTA+columnincrement},'3']);
                end
                
                % Space between data
                wcDTA = wcDTA + columnincrement+2;
            end
        end

    end
end
    disp('Export ready!');
    close(h);
end
warning('on')


handles.exportcount = handles.exportcount + 1;
guidata(hObject,handles);

function meanmin_Callback(hObject, eventdata, handles)
% hObject    handle to meanmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.meanminval = str2double(get(hObject,'String'));
if isempty(handles.meanmaxval)
    handles.meanmaxval = handles.meanminval;
end
handles.meanperc = str2double(get(handles.editmeanperc, 'String'))/100;


% Data
for normcount = 1:length(fieldnames(handles.alldata))
    if ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "2D") && ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "3D")
        tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
        [~, idxmindist] = min(abs(tempdata(:,1) - handles.meanminval));
        [~, idxmaxdist] = min(abs(tempdata(:,1) - handles.meanmaxval));
        tempnormval = mean(tempdata(idxmindist:idxmaxdist,2));
        handles.alldata.(['Data', num2str(normcount)]).dataMEAN = handles.alldata.(['Data', num2str(normcount)]).dataNON/tempnormval*handles.meanperc; % Omnidata multiplied with 10 before (collect function)
    end
end
% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

listboxChooseData_Callback(handles.listboxChooseData, [], handles);
% Hints: get(hObject,'String') returns contents of meanmin as text
%        str2double(get(hObject,'String')) returns contents of meanmin as a double


% --- Executes during object creation, after setting all properties.
function meanmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function meanmax_Callback(hObject, eventdata, handles)
% hObject    handle to meanmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.meanmaxval = str2double(get(hObject,'String'));
handles.meanperc = str2double(get(handles.editmeanperc, 'String'))/100;
if isempty(handles.meanminval)
    handles.meanminval = handles.meanmaxval;
end

% Data
for normcount = 1:length(fieldnames(handles.alldata))
    if ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "2D") && ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "3D") 
        tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
        [~, idxmindist] = min(abs(tempdata(:,1) - handles.meanminval));
        [~, idxmaxdist] = min(abs(tempdata(:,1) - handles.meanmaxval));
        tempnormval = mean(tempdata(idxmindist:idxmaxdist,2));
        handles.alldata.(['Data', num2str(normcount)]).dataMEAN = handles.alldata.(['Data', num2str(normcount)]).dataNON/tempnormval*handles.meanperc; % Omnidata multiplied with 10 before (collect function)
    end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);
listboxChooseData_Callback(handles.listboxChooseData, [], handles);
% Hints: get(hObject,'String') returns contents of meanmax as text
%        str2double(get(hObject,'String')) returns contents of meanmax as a double


% --- Executes during object creation, after setting all properties.
function meanmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editexport_Callback(hObject, eventdata, handles)
% hObject    handle to editexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = get(hObject, 'String');
if strcmp(temp, handles.myfile)
   % Do nothing, filename has not changed
else
   % Filename changed -> zero export counter. NOTE: If you re-use the same
   % filename twice after using another filename between them, the previous
   % files (with the same filename) will get overwritten
   handles.exportcount = 0;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of editexport as text
%        str2double(get(hObject,'String')) returns contents of editexport as a double


% --- Executes during object creation, after setting all properties.
function editexport_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editmeanperc_Callback(hObject, eventdata, handles)
% hObject    handle to editmeanperc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.meanmaxval)
    handles.meanmaxval = handles.meanminval;
end

handles.meanperc = str2double(get(hObject, 'String'))/100;

% Data
for normcount = 1:length(fieldnames(handles.alldata))
    if ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "2D") && ~strcmp(handles.alldata.(['Data', num2str(normcount)]).datatype, "3D")
        tempdata = handles.alldata.(['Data', num2str(normcount)]).interpolated;
        [~, idxmindist] = min(abs(tempdata(:,1) - handles.meanminval));
        [~, idxmaxdist] = min(abs(tempdata(:,1) - handles.meanmaxval));
        tempnormval = mean(tempdata(idxmindist:idxmaxdist,2));
        handles.alldata.(['Data', num2str(normcount)]).dataMEAN = handles.alldata.(['Data', num2str(normcount)]).dataNON/tempnormval*handles.meanperc; % Omnidata multiplied with 10 before (collect function)
    end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

listboxChooseData_Callback(handles.listboxChooseData, [], handles);
% Hints: get(hObject,'String') returns contents of editmeanperc as text
%        str2double(get(hObject,'String')) returns contents of editmeanperc as a double


% --- Executes during object creation, after setting all properties.
function editmeanperc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmeanperc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonAddData.
function pushbuttonAddData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic;
datatype = 'unknown';
try
    handles.mydir = dir(uigetdir);
catch
    handles.mydir = {};
end
% check if cancel is pressed
if ~isempty(handles.mydir)
%dirdir = handles.mydir;
%save('mydir.mat', 'dirdir');
guidata(hObject, handles);
skipped = 0;
filesindir = 0;

% Check datatype by reading 1st lines
for i = 1:length(handles.mydir)
   fid = fopen([handles.mydir(1).folder, '\', handles.mydir(i).name]);
   if handles.mydir(i).bytes > 1 
        testline = fgetl(fid);
        if contains(lower(testline), 'depth')
           % PTW.txt
           datatype = 'PTW';
        elseif contains(lower(testline), 'begin')
            datatype = 'mcc';
        elseif contains(testline, 'ascii')
            % Omni.txt
            handles.mydir(i).name
            datatype = 'omni';
        elseif contains(handles.mydir(i).name, 'ddose')
            % 3d dose file
            datatype = '3ddose';
            filesindir = filesindir + 1;
        elseif contains(handles.mydir(i).name, 'dcm')
            % 3d dose file
            datatype = 'dcm';
            handles.mydir(i).name
            filesindir = filesindir + 1;
        end
   else
       skipped = skipped + 1;
   end
end

fclose(fid);

% PTW data
if strcmp(datatype, 'PTW')
    pushbuttonPTW_Callback(hObject,eventdata,handles);
    handles = guidata(hObject);
    fields = fieldnames(handles.PTWdata);
    for i = 1:length(fields)
       handles.alldata.(['Data',num2str(handles.datacount)]) = handles.PTWdata.(fields{i});
       handles.alldata.(['Data',num2str(handles.datacount)]).directory = handles.mydir(1).folder;
       handles.datacount = handles.datacount + 1;
    end
    

 
elseif strcmp(datatype, 'mcc')
    handles.normalizationdist = str2double(get(handles.edit2, 'String'));
    handles.normalizationperc = str2double(get(handles.edit6, 'String'));
    mccdata = struct;
    mcclen = 0;
    % Collect file identifiers from files that are not directories
    filecount = 1;
    for j = 1:length(handles.mydir)
        if handles.mydir(j).bytes > 1 
            filepaths{filecount} = [handles.mydir(j).folder, '\', handles.mydir(j).name];
            filecount = filecount + 1;
        end
    end
    mccdata = readmccdata(filepaths, handles.normalizationdist, handles.normalizationperc);
    nams = fieldnames(mccdata);
    try
        datab4 = length(fieldnames(handles.alldata));
    catch
        datab4 = 0; 
    end
    for i = 1:length(nams)
        templist = ['mcc', num2str(i),'_',mccdata.(['mcc', num2str(i)]).name];
        temporig = templist;
        handles.alldata.(['Data', num2str(datab4 + i)]) = mccdata.(['mcc', num2str(i)]);
        handles.alldata.(['Data', num2str(datab4 + i)]).directory = handles.mydir(1).folder;
        handles.listboxArray{end+1} = strcat(templist);
        handles.origlist{end+1} = strcat(temporig);
        
        % Field params
        try
            [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(datab4 + i)]).interpolated(:,1), 'profile', handles.caxcorrection);
            handles.alldata.(['Data', num2str(datab4 + i)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
            handles.alldata.(['Data', num2str(datab4 + i)]).datatype = 'profile';
            % If symmetry of the profile is greater than 300% -> Try pdd
            if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.alldata.(['Data', num2str(datab4 + i)]).params)) > 0
                [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(datab4 + i)]).interpolated(:,1));
                handles.alldata.(['Data', num2str(datab4 + i)]).params = [R100, R80, R50, D100, D200, Ds, J1020];
                handles.alldata.(['Data', num2str(datab4 + i)]).datatype = 'pdd';
            end
        catch
            [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(datab4 + i)]).interpolated(:,1));
            handles.alldata.(['Data', num2str(datab4 + i)]).params = [R100, R80, R50, D100, D200, Ds, J1020];
            handles.alldata.(['Data', num2str(datab4 + i)]).datatype = 'pdd';
        end
        
        handles.datacount = handles.datacount + 1;
    end

    % Add files to listbox array
    handles.idx = zeros(1,length(fieldnames(handles.alldata)));
    set(handles.listboxChooseData, 'String', handles.listboxArray);
    handles.order = zeros(1,length(handles.idx));
    
    
    
% OmniPro
elseif strcmp(datatype, 'omni')
     ff = msgbox('Wait! Window will close automatically');
     pushbuttonOmni_Callback(hObject,eventdata,handles); 
     handles = guidata(hObject);
     fields = fieldnames(handles.omnidata);
     for i = 1:length(fieldnames(handles.omnidata))
       handles.alldata.(['Data',num2str(handles.datacount)]) = handles.omnidata.(fields{i});
       handles.alldata.(['Data',num2str(handles.datacount)]).directory = handles.mydir(1).folder;
       handles.datacount = handles.datacount + 1;
     end
     close(ff);
% 3ddose and DICOM data
elseif strcmp(datatype, '3ddose') || strcmp(datatype, 'dcm')
    ff = msgbox('Wait! Window will close automatically');
    % Get normalization values
    handles.normalizationdist = str2double(get(handles.edit2, 'String'));
    handles.normalizationperc = str2double(get(handles.edit6, 'String'));
    % Define maximum dimension size for data displaying
    maxdimsize = 120;
    
    for i = 1:filesindir
        truncate = 0;
        x = [];
        y = [];
        z = [];
        
        pathandname = [handles.mydir(1).folder, '\',handles.mydir(skipped + i).name];
        if strcmp(datatype, '3ddose')
            % 3Ddose files
            %[dose_3d, derror_3d, bound_x, bound_y, bound_z] = DOSXYZ3ddoseReader(pathandname);
            [dose_3d, derror_3d, bound_x, bound_y ,bound_z] = get3ddose(1,pathandname);
            disp('3ddose read');
            for j = 1:length(bound_x)-1
                x(j) = (bound_x(j+1) + bound_x(j))/2;
            end

            for j = 1:length(bound_y)-1
                y(j) = (bound_y(j+1) + bound_y(j))/2;
            end

            for j = 1:length(bound_z)-1
                z(j) = (bound_z(j+1) + bound_z(j))/2;
            end
        else
            % Dicom files
            [dose_3d, x, y, z, dataunit] = readdcmfiles(handles.mydir(1).folder, handles.mydir(skipped + i).name);
            handles.alldata.(['Data',num2str(handles.datacount)]).DataUnit = dataunit;
            disp('DICOM data read');
            
        end
        filename = handles.mydir(skipped + i).name;
%         x = bound_x;
%         y = bound_y;
%         z = bound_z;
        % Compute centers
        
        xtrunc = x;
        ytrunc = y;
        ztrunc = z;
        
        % Truncate X size for display, if data is too big
        if length(x) > maxdimsize
            xstep = floor(length(x)/maxdimsize);
            xtrunc = x(1:xstep:end);
%             while xtrunc(end) ~= x(end)
%                 xstep = xstep+1;
%                 xtrunc = x(1:xstep:end);
%             end
            truncate = 1;
        else
            xstep = 1;
        end

        % Truncate Y size for display, if data is too big
        if length(y) > maxdimsize
            ystep = floor(length(y)/maxdimsize);
            ytrunc = y(1:ystep:end);
%             while ytrunc(end) ~= y(end)
%                 ystep = ystep+1;
%                 ytrunc = y(1:ystep:end);
%             end
            truncate = 1;
        else
            ystep = 1;
        end
        
        % Truncate Z size for display, if data is too big
        if length(z) > maxdimsize
            zstep = floor(length(z)/maxdimsize);
            ztrunc = z(1:zstep:end);
%             while ztrunc(end) ~= z(end)
%                 zstep = zstep+1;
%                 ztrunc = z(1:zstep:end);
%             end
            truncate = 1;
        else
            zstep = 1;
        end
        
        % Add basic parameters
        handles.alldata.(['Data',num2str(handles.datacount)]).xpos = x;
        handles.alldata.(['Data',num2str(handles.datacount)]).ypos = y;
        handles.alldata.(['Data',num2str(handles.datacount)]).zpos = z;
        handles.alldata.(['Data',num2str(handles.datacount)]).data = dose_3d;
        try
            handles.alldata.(['Data',num2str(handles.datacount)]).error = derror_3d;
        catch 
        end
        handles.alldata.(['Data',num2str(handles.datacount)]).name = filename;
        handles.alldata.(['Data',num2str(handles.datacount)]).Plane = 'XYZ';
        handles.alldata.(['Data',num2str(handles.datacount)]).datatype = '3D';
        handles.alldata.(['Data',num2str(handles.datacount)]).Energy = 'Unknown';
        handles.alldata.(['Data',num2str(handles.datacount)]).FieldSize = 'Unknown';
        handles.alldata.(['Data',num2str(handles.datacount)]).directory = handles.mydir(1).folder;
        % Truncate size for display, if data is too big (e.g. > 100/dim)
        if truncate == 1
            disp('Data truncated to speedup display');
            
            
            Vdose = dose_3d(1:ystep:length(dose_3d(:,1,1)), 1:xstep:length(dose_3d(1,:,1)), 1:zstep:length(dose_3d(1,1,:)));
            try
                Verror = derror_3d(1:ystep:length(derror_3d(:,1,1)), 1:xstep:length(derror_3d(1,:,1)), 1:zstep:length(derror_3d(1,1,:)));
            catch
            end
            
            handles.alldata.(['Data',num2str(handles.datacount)]).Displaydata = Vdose;
            try
                handles.alldata.(['Data',num2str(handles.datacount)]).Displayerror = Verror;
            catch
            end
            handles.alldata.(['Data',num2str(handles.datacount)]).Displayxpos = xtrunc;
            handles.alldata.(['Data',num2str(handles.datacount)]).Displayypos = ytrunc;
            handles.alldata.(['Data',num2str(handles.datacount)]).Displayzpos = ztrunc;
        else
            handles.alldata.(['Data',num2str(handles.datacount)]).Displaydata = dose_3d;
            try
                handles.alldata.(['Data',num2str(handles.datacount)]).Displayerror = derror_3d;
            catch
            end
            handles.alldata.(['Data',num2str(handles.datacount)]).Displayxpos = x;
            handles.alldata.(['Data',num2str(handles.datacount)]).Displayypos = y;
            handles.alldata.(['Data',num2str(handles.datacount)]).Displayzpos = z;
        end
        % Find midaxis for 3D display
        handles.alldata.(['Data',num2str(handles.datacount)]).Depthraw = round(length(handles.alldata.(['Data',num2str(handles.datacount)]).data(1,1,:))/2);
        handles.alldata.(['Data',num2str(handles.datacount)]).Depthdisplay = round(length(handles.alldata.(['Data',num2str(handles.datacount)]).Displaydata(1,1,:))/2);
        
        handles.datacount = handles.datacount + 1;
    end

    
    % Try if previous data has been imported
    if ~isempty(handles.origlist)
        for i = handles.datacount-filesindir:handles.datacount-1
            templist = ['d3d', num2str(i), '_', handles.alldata.(['Data',num2str(i)]).name];
            temporig = ['d3d', num2str(i), '_', handles.alldata.(['Data',num2str(i)]).name];

            handles.listboxArray{end+1} = strcat(templist);
            handles.origlist{end+1} = strcat(temporig);

%             try
%                 [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(i)]).dataCAX, handles.alldata.(['Data', num2str(i)]).pos, '3ddata', handles.caxcorrection);
%                 handles.alldata.(['Data', num2str(i)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
%                 if sym > 300
%                    [R100, R80, R50, D100, D200, J1020] = pddparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,2),...
%                         handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,1));
%                     handles.alldata.(['Data', num2str(datab4 + i)]).params = [R100, R80, R50, D100, D200, J1020]; 
%                 end
%             catch
%                 [R100, R80, R50, D100, D200, J1020] = pddparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,2),...
%                 handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,1));
%                  handles.alldata.(['Data', num2str(datab4 + i)]).params = [R100, R80, R50, D100, D200, J1020];
%                 handles.alldata.(['Data', num2str(i)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
%                 disp(lasterror);
%             end
            handles.idx(end+1) = 0;
        end
    
    % No previous data -> import all from start
    else
        for i = 1:filesindir
            templist = ['d3d', num2str(i), '_', handles.alldata.(['Data',num2str(i)]).name];
            handles.listboxArray{i} = templist;
%             try
%                 [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.PTWdata.(['Data', num2str(i)]).dataCAX, handles.PTWdata.(['Data', num2str(i)]).pos, "3ddata", handles.caxcorrection);
%                 handles.alldata.(['Data', num2str(i)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
%                 if sym > 300
%                    [R100, R80, R50, D100, D200, J1020] = pddparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,2),...
%                     handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,1));
%                     handles.alldata.(['Data', num2str(datab4 + i)]).params = [R100, R80, R50, D100, D200, J1020]; 
%                 end
%             catch
%                 [R100, R80, R50, D100, D200, J1020] = pddparams(handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,2),...
%                 handles.alldata.(['Data', num2str(datab4 + i)]).interpolatedCAX(:,1));
%                 handles.alldata.(['Data', num2str(datab4 + i)]).params = [R100, R80, R50, D100, D200, J1020];
%                 handles.alldata.(['Data', num2str(i)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
%             end
        end
        
        handles.origlist = handles.listboxArray;
        handles.idx = zeros(1,filesindir);
    end
    % Add files to listbox array
    set(handles.listboxChooseData, 'String', handles.listboxArray);
    
    handles.order = zeros(1,length(handles.idx));
    close(ff);
    
elseif strcmp(datatype, 'dcm')
    
else
    disp('Unknown datatype');
end

toc;
%savestruct = handles.alldata;

%save('alladata.mat', 'savestruct');
end
guidata(hObject, handles);




% --- Executes on button press in pushbuttonAddLine.
function pushbuttonAddLine_Callback(hObject, eventdata, handles)

% Prompt XYZ boundaries
prompt = {'X [cm]:','Y [cm]:','Z [cm]:'};
dlgtitle = 'Inputs: Value or Empty = whole axis or a:b = from a to b.';
definput = {'','',''};
dims = [1 80];
newline = inputdlg(prompt,dlgtitle,dims,definput);
% check if cancel is pressed
if ~isempty(newline)
    
newlinematx = cell2mat(newline(1));
newlinematy = cell2mat(newline(2));
newlinematz = cell2mat(newline(3));

% Create lines from the chosen data
indices = find(handles.idx);

for i = 1:length(indices)
    currentnams = fieldnames(handles.alldata.(['Data', num2str(indices(i))]));
    % Direction and profile type
    if (isempty(newlinematx) || contains(newlinematx, ':')) && contains(handles.alldata.(['Data', num2str(indices(i))]).Plane,'X')
        datatype = 'profile';
        profdim = 'xpos';
        linedir = 1;
    elseif (isempty(newlinematy) || contains(newlinematy, ':')) && contains(handles.alldata.(['Data', num2str(indices(i))]).Plane,'Y')
        datatype = 'profile';
        profdim = 'ypos';
        linedir = 2;
    elseif (isempty(newlinematz) || contains(newlinematz, ':')) && contains(handles.alldata.(['Data', num2str(indices(i))]).Plane,'Z')
        profdim = 'zpos';
        datatype = 'profile';
        linedir = 3;
    end
    
   
   % Set final extracted position values as empty. If not chosen the empty
   % values won't get printed to the title
   valXpos = [];
   valYpos = [];
   valZpos = [];
   % Get boundaries (indices)
   switch linedir
       case 1
           % X profile
           if ismember('ypos', currentnams)
               [~ ,closestYval] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).ypos - str2double(newline{2})));
               valYpos = handles.alldata.(['Data', num2str(indices(i))]).ypos(closestYval);
           end
           
           if ismember('zpos', currentnams)
               [~ ,closestZval] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).zpos - str2double(newline{3})));
               valZpos = handles.alldata.(['Data', num2str(indices(i))]).zpos(closestZval);
           end
       case 2
           % Y profile
           if ismember('xpos', currentnams)
                [~ ,closestXval] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).xpos - str2double(newline{1})));
                valXpos = handles.alldata.(['Data', num2str(indices(i))]).xpos(closestXval);
           end
               
           if ismember('zpos', currentnams)
                [~ ,closestZval] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).zpos - str2double(newline{3})));
                valZpos = handles.alldata.(['Data', num2str(indices(i))]).zpos(closestZval);
           end
       case 3
           % Z profile
           if ismember('xpos', currentnams)
                [~ ,closestXval] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).xpos - str2double(newline{1})));
                valXpos = handles.alldata.(['Data', num2str(indices(i))]).xpos(closestXval);
           end
           if ismember('ypos', currentnams)
                [~ ,closestYval] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).ypos - str2double(newline{2})));
                valYpos = handles.alldata.(['Data', num2str(indices(i))]).ypos(closestYval);
           end
   end
   dataflipped = 0;
   % Get range
   if contains(newline{linedir}, ':')
       charpos = find(newline{linedir} == ':');
       firstbound = str2double(newline{linedir}(1:charpos-1));
       [~, firstidx] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).(profdim) - firstbound));
       secondbound = str2double(newline{linedir}(charpos+1:end));
       [~, secondidx] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).(profdim) - secondbound));
   else
       [firstbound, firstidx] = min(handles.alldata.(['Data', num2str(indices(i))]).(profdim));
       [secondbound, secondidx] = max(handles.alldata.(['Data', num2str(indices(i))]).(profdim));
   end
   
   % Check that the idx order is right
   if secondidx < firstidx
      temp = firstidx;
      firstidx = secondidx;
      secondidx = temp;
      dataflipped = 1;
      disp('Descending order swapped to ascending while creating the 1D line');
   end
   madata = [];
   errordata = [];
   madata(:,1) = handles.alldata.(['Data', num2str(indices(i))]).(profdim)(firstidx:secondidx);
   errordata(:,1) = madata(:,1);
   
   % Reduce dimensions 2D -> 1D
   if length(handles.alldata.(['Data', num2str(indices(i))]).Plane) == 2
       if strcmp(handles.alldata.(['Data', num2str(indices(i))]).Plane, 'XY') && strcmp(profdim, 'xpos')
            madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).dataNON(closestYval, firstidx:secondidx);
            disp(['Data not interpolated while getting the line. Line extracted from: Y: ', num2str(valYpos), 'cm']);
            disp(['Requested Y: ', newline{2}, 'cm']);
            if ismember('error',currentnams)
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(closestYval, firstidx:secondidx);
            end
       elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).Plane, 'XY') && strcmp(profdim, 'ypos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).dataNON(firstidx:secondidx,closestXval);
           disp(['Data not interpolated while getting the line. Line extracted from: X: ', num2str(valXpos), 'cm']);
           disp(['Requested X: ', newline{1}, 'cm']);
           if ismember('error',currentnams)
               errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(firstidx:secondidx, closestXval);
           end
       elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).Plane, 'YZ') && strcmp(profdim, 'ypos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).dataNON(firstidx:secondidx, closestZval);
           disp(['Data not interpolated while getting the line. Line extracted from: Z: ', num2str(valYZos), 'cm']);
           disp(['Requested Z: ', newline{3}, 'cm']);
           if ismember('error',currentnams)
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(firstidx:secondidx, closestZval);
           end
       elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).Plane, 'YZ') && strcmp(profdim, 'zpos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).dataNON(closestYval, firstidx:secondidx);
           disp(['Data not interpolated while getting the line. Line extracted from: Y: ', num2str(valYpos), 'cm']);
           disp(['Requested Y: ', newline{2}, 'cm']);
           if ismember('error',currentnams)
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(closestYval, firstidx:secondidx);
           end
       elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).Plane, 'XZ') && strcmp(profdim, 'zpos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).dataNON(closestXval, firstidx:secondidx);
           disp(['Data not interpolated while getting the line. Line extracted from: X: ', num2str(valXpos), 'cm']);
           disp(['Requested X: ', newline{1}, 'cm']);
           if ismember('error',currentnams)
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(closestXval, firstidx:secondidx);
           end
       elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).Plane, 'XZ') && strcmp(profdim, 'xpos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).dataNON(firstidx:secondidx, closestZval);
           disp(['Data not interpolated while getting the line. Line extracted from: Z: ', num2str(valZpos), 'cm']);
           disp(['Requested Z: ', newline{3}, 'cm']);
           if ismember('error',currentnams)
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(firstidx:secondidx, closestZval);
           end
       end
       
       
   % Reduce dimensions 3D -> 1D
   else
       if strcmp(profdim, 'xpos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).data(closestYval, firstidx:secondidx, closestZval);
           disp(['Data not interpolated while getting the line. Line extracted from: Y: ', num2str(valYpos), 'cm', ' Z: ', num2str(valZpos), 'cm']);
           disp(['Requested Y: ', newline{2}, 'cm', ' Z:', newline{3}, 'cm']);
            try
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(closestYval, firstidx:secondidx, closestZval);
            catch
                % No error data
            end
       elseif strcmp(profdim, 'ypos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).data(firstidx:secondidx, closestXval, closestZval);
           disp(['Data not interpolated while getting the line. Line extracted from: X: ', num2str(valXpos), 'cm', ' Z: ', num2str(valZpos), 'cm']);
           disp(['Requested X: ', newline{1}, 'cm', ' Z:', newline{3}, 'cm']);
            try
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(firstidx:secondidx, closestXval, closestZval);
            catch
                % No error data
            end
       elseif strcmp(profdim, 'zpos')
           madata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).data(closestYval, closestXval,  firstidx:secondidx);
           disp(['Data not interpolated while getting the line. Line extracted from: X: ', num2str(valXpos), 'cm', ' Y: ', num2str(valYpos), 'cm']);
           disp(['Requested X: ', newline{1}, 'cm', ' Y:', newline{2}, 'cm']);
            try
                errordata(:,2) = handles.alldata.(['Data', num2str(indices(i))]).error(closestYval, closestXval,  firstidx:secondidx);
            catch
                % No error data
            end
       end
   end
   if dataflipped == 1
      madata(:,2) = flip(madata(:,2));
      madata(:,1) = flip(madata(:,1));
   end
   
    % Create new line
    parentname = handles.alldata.(['Data',num2str(indices(i))]).name;
    templist = ['Line', num2str(handles.newlinecount), '_', handles.alldata.(['Data',num2str(indices(i))]).name,'_', 'X:', num2str(valXpos), '_Y:',...
        num2str(valYpos), '_Z:', num2str(valZpos)];
    templist = strrep(templist, 'Plane', '');
    temporig = templist;
    % Add name to listbox array
    handles.listboxArray{end+1} = strcat(templist);
    handles.origlist{end+1} = strcat(temporig);
    % Add info from the original data
    enrg = handles.alldata.(['Data',num2str(indices(i))]).Energy;
    try
        FS = handles.alldata.(['Data',num2str(indices(i))]).FieldSizeIn;
    catch
        FS = handles.alldata.(['Data',num2str(indices(i))]).FieldSize;
    end
    
    % Add the old directory
    datadir = handles.alldata.(['Data', num2str(indices(i))]).directory;
    % Get the line
    tempo = collect3ddata(templist, parentname,enrg, FS, madata, datatype, datadir, handles.newlinecount, handles.normalizationdist, handles.normalizationperc, handles.normalizationtoval);
    try
        tempo.(['NEW', num2str(handles.newlinecount)]).DataUnit = handles.alldata.(['Data',num2str(indices(i))]).DataUnit;
    catch
        disp('Error in getting old data unit (not existing). Set as Unknown');
    end
    
    handles.alldata.(['Data', num2str(length(handles.idx)+1)]) = tempo.(['NEW', num2str(handles.newlinecount)]);
    if ismember('error',currentnams)
        handles.alldata.(['Data', num2str(length(handles.idx)+1)]).error = errordata(:,2);
    end
    handles.newlinecount = handles.newlinecount + 1;
    %alladata = handles.alldata;
    %save('alladata.mat', 'alladata');
    % Compute line parameters
    try
        [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(length(handles.idx)+1)]).interpolated(:,2),...
            handles.alldata.(['Data', num2str(length(handles.idx)+1)]).interpolated(:,1), '3ddata', handles.caxcorrection);
        handles.alldata.(['Data', num2str(length(handles.idx)+1)]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
        handles.alldata.(['Data', num2str(length(handles.idx)+1)]).datatype = 'profile';
        if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.alldata.(['Data', num2str(length(handles.idx)+1)]).params)) > 0
           [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(length(handles.idx)+1)]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(length(handles.idx)+1)]).interpolated(:,1));
            handles.alldata.(['Data', num2str(length(handles.idx)+1)]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
            handles.alldata.(['Data', num2str(length(handles.idx)+1)]).datatype = 'pdd';
            disp('Symmetry too low -> Try PDD');
        end
    catch
        disp('Trying PDD data');
        [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(length(handles.idx)+1)]).interpolatedCAX(:,2),...
                handles.alldata.(['Data',num2str(length(handles.idx)+1)]).interpolatedCAX(:,1));
            handles.alldata.(['Data', num2str(length(handles.idx)+1)]).params = [R100, R80, R50, D100, D200, Ds, J1020];
            handles.alldata.(['Data', num2str(length(handles.idx)+1)]).datatype = 'pdd';
    end
    
   % Increase indices
   handles.idx(end+1) = 0;
   handles.datacount = handles.datacount + 1;
   handles.order(end+1) = 0;
end
%alladata = handles.alldata;
%save('alladata.mat', 'alladata');
set(handles.listboxChooseData, 'String', handles.listboxArray);
end
% Update handles
guidata(hObject, handles);







% --- Executes on button press in pushbuttonAddPlane.
function pushbuttonAddPlane_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddPlane (see GCBO)

prompt = {'Plane (XY, XZ or YZ)', 'Plane depth [cm]:'};
dlgtitle = 'Choose plane and the depth As default X = horizontal, Y = vertical and Z- into the screen';
definput = {'XY', ''};
dims = [1 150];
newline = inputdlg(prompt,dlgtitle,dims,definput);

% check if cancel is pressed
if ~isempty(newline)
newplane = upper(newline{1});
planedepth = newline{2};
indices = find(handles.idx);

for i = 1:length(indices)
    
   currentnams = fieldnames(handles.alldata.(['Data', num2str(indices(i))]));
   currentdata = handles.alldata.(['Data', num2str(indices(i))]).data;
   goaldepth = str2double(planedepth);
   
   % Get boundaries (indices)
   % XY profile
   if strcmp(newplane,'ZY') || strcmp(newplane,'YZ')
       newplane = 'YZ';
        
       % Interpolate to get also planes between plane depth values
       
       dimdepth = handles.alldata.(['Data', num2str(indices(i))]).xpos;
       lend = length(dimdepth);
       [~,closestdepthidx] = min(abs(dimdepth - goaldepth));
       closestdepth = dimdepth(closestdepthidx);
       depthdiff = closestdepth - goaldepth;
       % Check if data ascending or descending
       if dimdepth(1) - dimdepth(2) < 0 % Ascending
           datadir = 'asc';
       else    % Descending
           datadir = 'desc';
       end
       if closestdepthidx < lend && ((((depthdiff < 0  && strcmp(datadir, 'asc'))) ||  (depthdiff > 0 && strcmp(datadir, 'desc'))))
          % value between this and next index
          disp('Requested point does not exist in data -> plane interpolated (message 1)')
          interval = abs(dimdepth(closestdepthidx) - dimdepth(closestdepthidx+1));
          lincoeff = abs(interval - depthdiff)/interval;
          closestplane = squeeze(currentdata(:,closestdepthidx,:));
          madata = closestplane + (squeeze(currentdata(:,closestdepthidx+1,:)) - closestplane)*lincoeff;
          % Error values
          if ismember('error', currentnams)
              error = handles.alldata.(['Data', num2str(indices(i))]).error;
              closesterror = squeeze(error(:,closestdepthidx,:));
              errordata = closesterror + (closesterror - squeeze(error(:,closestdepthidx+1,:)))*lincoeff;
          end
       elseif closestdepthidx > 1 && ((((depthdiff > 0 && strcmp(datadir, 'asc'))) ||  (depthdiff < 0  && strcmp(datadir, 'desc'))))
            % value between this and next index
          disp('Requested point does not exist in data -> plane interpolated (message 2)')
          interval = abs(dimdepth(closestdepthidx) - dimdepth(closestdepthidx-1));
          lincoeff = abs(interval - depthdiff)/interval;
          closestplane = squeeze(currentdata(:,closestdepthidx-1,:));
          madata = closestplane + (squeeze(currentdata(:,closestdepthidx,:)) - closestplane)*lincoeff;
          % Error values
          if ismember('error', currentnams)
              disp('Closest plane chosen');
              error = handles.alldata.(['Data', num2str(indices(i))]).error;
              closesterror = squeeze(error(:,closestdepthidx-1,:));
              errordata = closesterror + (squeeze(error(:,closestdepthidx,:)) - closestplane)*lincoeff;
          end
       else
           % Find closest depth plane from existing data
           [~ ,closestdepth] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).xpos - str2double(planedepth)));
           madata = squeeze(handles.alldata.(['Data', num2str(indices(i))]).data(:, closestdepth, :));
           if ismember('error', currentnams)
                errordata = squeeze(handles.alldata.(['Data', num2str(indices(i))]).error(:, closestdepth, :));
           end
       end
       
       % Reduce dimensions 3D -> 2D in YZ plane
       
       madatahor = handles.alldata.(['Data', num2str(indices(i))]).zpos;
       madatavert = handles.alldata.(['Data', num2str(indices(i))]).ypos;
       
       
   
   elseif strcmp(newplane,'XZ') || strcmp(newplane,'ZX')
       newplane = 'XZ';
       % Interpolate to get also planes between plan values
       
       dimdepth = handles.alldata.(['Data', num2str(indices(i))]).ypos;
       lend = length(dimdepth);
       [~,closestdepthidx] = min(abs(dimdepth - goaldepth));
       closestdepth = dimdepth(closestdepthidx);
       depthdiff = closestdepth - goaldepth;
       % Check if data ascending or descending
       if dimdepth(1) - dimdepth(2) < 0 % Ascending
           datadir = 'asc';
       else    % Descending
           datadir = 'desc';
       end
       if closestdepthidx < lend && ((((depthdiff < 0  && strcmp(datadir, 'asc'))) ||  (depthdiff > 0 && strcmp(datadir, 'desc'))))
          % value between this and next index
          disp('Requested point does not exist in data -> plane interpolated (message 1.1)')
          interval = abs(dimdepth(closestdepthidx) - dimdepth(closestdepthidx+1));
          lincoeff = abs(interval - depthdiff)/interval;
          closestplane = squeeze(currentdata(closestdepthidx,:,:));
          madata = closestplane + (squeeze(currentdata(closestdepthidx+1,:,:)) - closestplane)*lincoeff;
          % Error values
          if ismember('error', currentnams)                                         % Do same for the error data
              error = handles.alldata.(['Data', num2str(indices(i))]).error;
              closesterror = squeeze(error(closestdepthidx,:,:));
              errordata = closesterror + (closesterror - squeeze(error(closestdepthidx+1,:,:)))*lincoeff;
          end
       elseif closestdepthidx > 1 && ((((depthdiff > 0 && strcmp(datadir, 'asc'))) ||  (depthdiff < 0  && strcmp(datadir, 'desc'))))
            % value between this and next index
            disp('Requested point does not exist in data -> plane interpolated (message 2.1)')
          interval = abs(dimdepth(closestdepthidx) - dimdepth(closestdepthidx-1));
          lincoeff = abs(interval - depthdiff)/interval;
          closestplane = squeeze(currentdata(closestdepthidx-1,:,:));
          madata = closestplane + (squeeze(currentdata(closestdepthidx,:,:)) - closestplane)*lincoeff;
          % Error values
          if ismember('error', currentnams)
              disp('Closest plane chosen');
              error = handles.alldata.(['Data', num2str(indices(i))]).error;
              closesterror = squeeze(error(closestdepthidx-1,:,:));
              errordata = closesterror + (squeeze(error(closestdepthidx,:,:)) - closestplane)*lincoeff;
          end
       else
           % Find closest depth plane from existing data
           [~ ,closestdepth] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).ypos - str2double(planedepth)));
           madata = squeeze(handles.alldata.(['Data', num2str(indices(i))]).data(closestdepth, :, :));
           if ismember('error', currentnams)
                errordata = squeeze(handles.alldata.(['Data', num2str(indices(i))]).error(closestdepth,:,:));
           end
       end
       
       % Reduce dimensions 3D -> 2D in XZ plane
       
       madatahor = handles.alldata.(['Data', num2str(indices(i))]).zpos;
       madatavert = handles.alldata.(['Data', num2str(indices(i))]).xpos;
       
   
       
   else
       newplane = 'XY';
       % Interpolate to get also planes between plan values
       
       dimdepth = handles.alldata.(['Data', num2str(indices(i))]).zpos;
       lend = length(dimdepth);
       [~,closestdepthidx] = min(abs(dimdepth - goaldepth));
       closestdepth = dimdepth(closestdepthidx);
       depthdiff = closestdepth - goaldepth;
       % Check if data ascending or descending
       if dimdepth(1) - dimdepth(2) < 0 % Ascending
           datadir = 'asc';
       else    % Descending
           datadir = 'desc';
       end
       if closestdepthidx < lend && ((((depthdiff < 0  && strcmp(datadir, 'asc'))) ||  (depthdiff > 0 && strcmp(datadir, 'desc'))))
          % value between this and next index
          disp('Requested point does not exist in data -> plane interpolated (message 1.2)')
          interval = abs(dimdepth(closestdepthidx) - dimdepth(closestdepthidx+1));
          lincoeff = abs(interval - depthdiff)/interval;
          closestplane = squeeze(currentdata(:,:,closestdepthidx));
          madata = closestplane + (currentdata(:,:,closestdepthidx+1) - closestplane)*lincoeff;
          % Error values
          if ismember('error', currentnams)
              error = handles.alldata.(['Data', num2str(indices(i))]).error;
              closesterror = squeeze(error(:,:,closestdepthidx));
              errordata = closesterror + (closesterror - (error(:,:,closestdepthidx+1)))*lincoeff;
          end
       elseif closestdepthidx > 1 && ((((depthdiff > 0 && strcmp(datadir, 'asc'))) ||  (depthdiff < 0  && strcmp(datadir, 'desc'))))
            % value between this and next index
            disp('Requested point does not exist in data -> plane interpolated (message 2.2)')
          interval = abs(dimdepth(closestdepthidx) - dimdepth(closestdepthidx-1));
          lincoeff = abs(interval - depthdiff)/interval;
          closestplane = currentdata(:,:,closestdepthidx-1);
          madata = closestplane + (currentdata(:,:,closestdepthidx) - closestplane)*lincoeff;
          % Error values
          if ismember('error', currentnams)
              error = handles.alldata.(['Data', num2str(indices(i))]).error;
              closesterror = error(:,:,closestdepthidx-1);
              errordata = closesterror + (closesterror - error(:,:,closestdepthidx))*lincoeff;
          end
       else
           disp('Closest plane from data chosen')
           % Find closest depth plane from existing data
           [~ ,closestdepth] = min(abs(handles.alldata.(['Data', num2str(indices(i))]).zpos - str2double(planedepth)));
           madata = squeeze(handles.alldata.(['Data', num2str(indices(i))]).data(:,:,closestdepth));
           
           if ismember('error', currentnams)
                errordata = (handles.alldata.(['Data', num2str(indices(i))]).error(:,:,closestdepth));
           end
       end
       
       % Reduce dimensions 3D -> 2D in XY plane
       % Order reverser (hor -> vert) because XY plane is named in
       % different order with relative to the order in data compared to YZ
       % and XZ
       madatahor = handles.alldata.(['Data', num2str(indices(i))]).xpos;
       madatavert = handles.alldata.(['Data', num2str(indices(i))]).ypos;
       % Error values
   end

    % Create new plane
    templist = ['Plane', num2str(handles.newlinecount), '_', [handles.alldata.(['Data',num2str(indices(i))]).name,'mod'],'_',...
        newplane, '_', planedepth, '_cm'];
%     try
%         templist = strrep(cell2mat(templist), '.', '_');
%     catch
%         templist = strrep(templist, '.','_');
%     end
    temporig = templist;
    % Add name to listarray
    handles.listboxArray{end+1} = strcat(templist);
    handles.origlist{end+1} = strcat(temporig);
    % Add info from the original data
    enrg = handles.alldata.(['Data',num2str(indices(i))]).Energy;
    parentname = handles.alldata.(['Data',num2str(indices(i))]).name;
    try
        FS = handles.alldata.(['Data',num2str(indices(i))]).FieldSizeIn;
    catch
        FS = handles.alldata.(['Data',num2str(indices(i))]).FieldSize;
    end
    % Add the old directory
    datadir = handles.alldata.(['Data', num2str(indices(i))]).directory;
    tempo = modifynewplane(madata, madatahor, madatavert, templist, parentname, datadir, enrg, FS, handles.newlinecount, newplane, handles.normalizationtoval);
    % Get the line
    
    try
        tempo.(['NEW', num2str(handles.newlinecount)]).DataUnit = handles.alldata.(['Data',num2str(indices(i))]).DataUnit;
    catch
        disp('Error in finding the DataUnit');
    end
    
    %tempo = collect3ddata(templist, enrg, FS, madata, datatype, handles.newlinecount, handles.normalizationdist, handles.normalizationperc);
    handles.alldata.(['Data', num2str(length(handles.idx)+1)]) = tempo.(['NEW', num2str(handles.newlinecount)]);
    handles.newlinecount = handles.newlinecount + 1;
    
   if ismember('error', currentnams)
        handles.alldata.(['Data', num2str(length(handles.idx)+1)]).error = errordata;
   end
   handles.alldata.(['Data', num2str(length(handles.idx)+1)]).depth = planedepth;
   % Increase indices
   handles.idx(end+1) = 0;
   handles.datacount = handles.datacount + 1;
   handles.order(end+1) = 0;
end
%alladata = handles.alldata;
%save('alladata.mat', 'alladata');
set(handles.listboxChooseData, 'String', handles.listboxArray);
end
guidata(hObject, handles);




% --- Executes on button press in checkboxGamma.
function checkboxGamma_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nogammacalc = get(hObject,'Value');

guidata(hObject, handles);



% --- Executes on button press in checkboxCaxCorr.
function checkboxCaxCorr_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCaxCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.caxcorrection = get(hObject,'Value');

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkboxCaxCorr


% --- Executes on button press in checkboxShowError.
function checkboxShowError_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showError = get(hObject,'Value');

% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Kuvaajaa ei p�ivitetty'); 
end
% Hint: get(hObject,'Value') returns toggle state of checkboxShowError


% --- Executes on button press in pushbuttonInvertAxis.
function pushbuttonInvertAxis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonInvertAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLABfor i = 1:length(indices)
    
indices = find(handles.idx); 

for i = 1:length(indices)
    if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D') &&  ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
        % Flip all data relative to its axis
        handles.alldata.(['Data', num2str(indices(i))]).dataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataCAX);
        handles.alldata.(['Data', num2str(indices(i))]).dataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAX);
        handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2) = flip(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2));
        handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,2) = flip(handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,2));
        % Try carch for datatypes which may not exist
        try
            handles.alldata.(['Data', num2str(indices(i))]).dataMEAN = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMEAN);
        catch
            disp('Mean normalized data not flipped (does not exist)');
        end
        try
            handles.alldata.(['Data', num2str(indices(i))]).dataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAN);
        catch
            disp('Manually normalized data not flipped (does not exist)');
        end
        
        handles.alldata.(['Data', num2str(indices(i))]).dataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).dataNON);
        
        try
            % Switch penumbras
%             tempPenumbra = handles.alldata.(['Data', num2str(indices(i))]).params(end);
%             handles.alldata.(['Data', num2str(indices(i))]).params(end) = handles.alldata.(['Data', num2str(indices(i))]).params(end-1);
%             handles.alldata.(['Data', num2str(indices(i))]).params(end-1) = tempPenumbra;
%             % Change sign for CAXdev
%             handles.alldata.(['Data', num2str(indices(i))]).params(1) = handles.alldata.(['Data', num2str(indices(i))]).params(1);
              [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1), 'profile', handles.caxcorrection);
              handles.alldata.(['Data', num2str(indices(i))]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
              handles.alldata.(['Data', num2str(indices(i))]).datatype = 'profile';
              if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.alldata.(['Data', num2str(indices(i))]).params)) > 0
                 [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                    handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1));
                handles.alldata.(['Data', num2str(indices(i))]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
                handles.alldata.(['Data', num2str(indices(i))]).datatype = 'pdd';
              end
        catch
            disp('No parameters to modify');
        end
    end
end


% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end


% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonShift.
function pushbuttonShift_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indices = find(handles.idx); 

for i = 1:length(indices)
    if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D') &&  ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
        handles.alldata.(['Data', num2str(indices(i))]).pos = handles.alldata.(['Data', num2str(indices(i))]).pos + handles.shiftaxis;
        handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1) = handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1) + handles.shiftaxis;
        handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,1) = handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,1) + handles.shiftaxis;
    end
    
    if ismember('params', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
        % Update cax deviation, CAX existing
        %try
        try
            [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                    handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1), 'profile', handles.caxcorrection);
            handles.alldata.(['Data', num2str(indices(i))]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
            handles.alldata.(['Data', num2str(indices(i))]).datatype = 'profile';
            if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.alldata.(['Data', num2str(indices(i))]).params)) > 0
               [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1));
            handles.alldata.(['Data', num2str(indices(i))]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
            handles.alldata.(['Data', num2str(indices(i))]).datatype = 'pdd';
            end
        catch
           disp('Error in computing field parameters after shift! Parameters not computed'); 
        end
        %handles.alldata.(['Data', num2str(indices(i))]).params(1) = handles.alldata.(['Data', num2str(indices(i))]).params(1) + handles.shiftaxis;
        %catch
         %  disp('Error while computing field parameters during shift');
        %   disp(lasterror);
        %end
    end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;
guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch ME
   disp('Error in shifting data. Figure not updated'); 
   disp(lasterror);
end



function editShift_Callback(hObject, eventdata, handles)
% hObject    handle to editShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.shiftaxis = str2double(get(hObject,'String'));

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of editShift as text
%        str2double(get(hObject,'String')) returns contents of editShift as a double


% --- Executes during object creation, after setting all properties.
function editShift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuFilters.
function popupmenuFilters_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFilters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
chosenfilt = contents{get(hObject,'Value')};
indices = find(handles.idx);

if strcmp(chosenfilt, 'Savitzky-Golay filter')
    prompt = {'Polynomial order','frame length (>= order)'};
    dlgtitle = 'Set required values for filtering';
    dims = [1 80];
    filtvalues = inputdlg(prompt,dlgtitle,dims);
elseif strcmp(chosenfilt, 'Scale with value (multiply)')
    prompt = {'Scaling factor'};
    dlgtitle = 'Choose the scaling factor';
    dims = [1 80];
    filtvalues = inputdlg(prompt,dlgtitle,dims);
end


% check if cancel is pressed
if ~isempty(filtvalues)
for i = 1:length(indices)
    tic;
    if strcmp(chosenfilt, 'Savitzky-Golay filter')
        if strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D')
            try
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
            catch
                
            end
            try
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
            catch
                
            end
            try
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
            catch
                
            end
            try
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
            catch
                
            end
            try
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
            catch
                
            end
                
            try
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataTOVAL = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataTOVAL,...
                    str2double(filtvalues{1}),str2double(filtvalues{2}));
            catch    
                
            end
            
        end
        try
            handles.alldata.(['Data', num2str(indices(i))]).dataCAX = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).dataCAX,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
        catch
            disp('No dataCAX');
        end
        try
        handles.alldata.(['Data', num2str(indices(i))]).dataMAX = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).dataMAX,...
            str2double(filtvalues{1}),str2double(filtvalues{2}));
        catch
            disp('No dataMAX');
        end
        try
        handles.alldata.(['Data', num2str(indices(i))]).dataNON = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).dataNON,...
            str2double(filtvalues{1}),str2double(filtvalues{2}));
        catch
            disp('No dataNON');
        end
        try
        handles.alldata.(['Data', num2str(indices(i))]).dataMAN = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).dataMAN,...
            str2double(filtvalues{1}),str2double(filtvalues{2}));
        catch
            disp('No dataMAN');
        end
        try
            handles.alldata.(['Data', num2str(indices(i))]).dataMEAN = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).dataMEAN,...
                str2double(filtvalues{1}),str2double(filtvalues{2}));
        catch
            disp('dataMEAN does not exist');
        end
        
        try
                handles.alldata.(['Data', num2str(indices(i))]).dataTOVAL = sgolayfilt(handles.alldata.(['Data', num2str(indices(i))]).dataTOVAL,...
                    str2double(filtvalues{1}),str2double(filtvalues{2}));
        catch    
           disp('dataTOVAL does not exist');
        end
        
    elseif strcmp(chosenfilt, 'Scale with value (multiply)')
        if strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D')
            
                handles.alldata.(['Data', num2str(indices(i))]).(['Displaydata', handles.normalizationtype]) = handles.alldata.(['Data', num2str(indices(i))]).(['Displaydata', handles.normalizationtype])*str2double(filtvalues{1});
                handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype]) = handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype])*str2double(filtvalues{1});
                
        elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
                if ismember(['Displaydata', handles.normalizationtype], fieldnames(handles.alldata.(['Data', num2str(indices(i))]))) && ~strcmp(handles.normalizationtype, 'NON')
                    handles.alldata.(['Data', num2str(indices(i))]).(['Displaydata', handles.normalizationtype]) = handles.alldata.(['Data', num2str(indices(i))]).(['Displaydata', handles.normalizationtype])*str2double(filtvalues{1});
                    handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype]) = handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype])*str2double(filtvalues{1});
                else
                    handles.alldata.(['Data', num2str(indices(i))]).data = handles.alldata.(['Data', num2str(indices(i))]).data*str2double(filtvalues{1});
                    handles.alldata.(['Data', num2str(indices(i))]).Displaydata = handles.alldata.(['Data', num2str(indices(i))]).Displaydata*str2double(filtvalues{1});
                end
                
        elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'profile') || strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'PDD')
                handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype]) = handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype])*str2double(filtvalues{1});
                handles.alldata.(['Data', num2str(indices(i))]).interpolated = handles.alldata.(['Data', num2str(indices(i))]).interpolated*str2double(filtvalues{1});
                
                if ismember('interpolatedCAX', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX = handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX*str2double(filtvalues{1});
                end
        end
            
    end
    disp(['Data', num2str(indices(i)), ' filtered in ', num2str(toc), ' seconds']);
end

%savethis = handles.alldata;
%save('alladata.mat', 'savethis');
% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in mirroring data. Figure not updated'); 
end

end

% --- Executes during object creation, after setting all properties.
function popupmenuFilters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFilters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mirrordataButton.
function mirrordataButton_Callback(hObject, eventdata, handles)
% hObject    handle to mirrordataButton (see GCBO)
% hObject    handle to pushbuttonInvertAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLABfor i = 1:length(indices)
    
indices = find(handles.idx);

for i = 1:length(indices)
    if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D') &&  ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
        handles.alldata.(['Data', num2str(indices(i))]).pos = handles.alldata.(['Data', num2str(indices(i))]).pos - handles.alldata.(['Data', num2str(indices(i))]).pos*2;
        handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,1) = handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,1) - handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,1)*2;
        handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1) = handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1) - handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1)*2;
         try
            % Switch penumbras
%             tempPenumbra = handles.alldata.(['Data', num2str(indices(i))]).params(end);
%             handles.alldata.(['Data', num2str(indices(i))]).params(end) = handles.alldata.(['Data', num2str(indices(i))]).params(end-1);
%             handles.alldata.(['Data', num2str(indices(i))]).params(end-1) = tempPenumbra;
            [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
            handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1), 'profile', handles.caxcorrection);
            handles.alldata.(['Data', num2str(indices(i))]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
            handles.alldata.(['Data', num2str(indices(i))]).datatype = 'profile';
            if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.alldata.(['Data', num2str(indices(i))]).params)) > 0
               [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1));
            handles.alldata.(['Data', num2str(indices(i))]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
            handles.alldata.(['Data', num2str(indices(i))]).datatype = 'pdd';
            end
        catch
            disp('No parameters to modify');
         end
        if ismember('params',lower(fieldnames(handles.alldata.(['Data', num2str(indices(i))]))))
            % Update cax deviation
            %handles.alldata.(['Data', num2str(indices(i))]).params(1) = handles.alldata.(['Data', num2str(indices(i))]).params(1)*-1;
            try
                [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                    handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1), 'profile', handles.caxcorrection);
                handles.alldata.(['Data', num2str(indices(i))]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
                handles.alldata.(['Data', num2str(indices(i))]).datatype = 'profile';
                if (abs(sym) > handles.symmetryThresHold && abs(penR-penL) > min([penR,penL])) || sum(isnan(handles.alldata.(['Data', num2str(indices(i))]).params)) > 0
                   [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                        handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1));
                    handles.alldata.(['Data', num2str(indices(i))]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
                    handles.alldata.(['Data', num2str(indices(i))]).datatype = 'pdd';
                end
            catch
                disp('Error in calculating parameters during mirroring');
            end
        end
    end
end



% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in mirroring data. Figure not updated'); 
end



% --- Executes on button press in checkboxTruncate.
function checkboxTruncate_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTruncate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.lowresvisualization = get(hObject,'Value');

if sum(handles.idx) > 0
    handles.redraw = 1;
end
guidata(hObject, handles);

%try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
%catch
%   disp('Error in displaying data during data truncation'); 
%end



function mousepos(hObject, ~,handles)


handles = guidata(hObject);

% Get mouse position
handles.mousepoint = get(gcf, 'Currentpoint');
handles.axespos = get(handles.axes1,'Position');
% Normalize mouseposition
handles.mousepoint(1) = handles.mousepoint(1)/handles.screendims(3);
handles.mousepoint(2) = handles.mousepoint(2)/handles.screendims(4);

% Set up button press function if the two overlap

guidata(hObject, handles);


function mousescroll(hObject, eventdata, handles)

handles = guidata(hObject);
datachsn = find(handles.idx);
% Get mouse position
handles.mousepoint = get(gcf, 'Currentpoint');
handles.axespos = get(handles.axes1,'Position');
% Normalize mouseposition
handles.mousepoint(1) = handles.mousepoint(1)/handles.screendims(3);
handles.mousepoint(2) = handles.mousepoint(2)/handles.screendims(4);
for i = 1:length(datachsn)
    % Set up button press function if the two overlap
    if handles.mousepoint(1) > handles.axespos(1) && handles.mousepoint(1) < handles.axespos(3)...
            && handles.mousepoint(2) > handles.axespos(2) && handles.mousepoint(2) < handles.axespos(4)  % Check that the mouse is placed over the figure
        if eventdata.VerticalScrollCount < 0
            if handles.lowresvisualization == 1 && handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay > 1
                handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay = handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay - 1;
            elseif handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw > 1 && handles.lowresvisualization == 0
                handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw = handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw - 1;
            end
        elseif eventdata.VerticalScrollCount > 0
            if handles.lowresvisualization == 1 && handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay < length(handles.alldata.(['Data',num2str(datachsn(i))]).Displaydata(1,1,:))
                handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay = handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay + 1;
            elseif handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw < length(handles.alldata.(['Data',num2str(datachsn(i))]).data(1,1,:))
                handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw = handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw + 1;
            end
        end
    end

% Update the figure with only re-drawing chosen data (excecuting in this function much faster than calling listbox display function)
    handles.redraw = 1;

    if contains(handles.alldata.(['Data', num2str(datachsn(i))]).datatype, '3D')
        ddisp = handles.alldata.(['Data',num2str(datachsn(i))]).Depthdisplay;
        draw = handles.alldata.(['Data',num2str(datachsn(i))]).Depthraw;

        % If 3D data -> plot initially mid axis XY-plane
        hold off;
        if handles.lowresvisualization == 1
            if ismember(['Displaydata', handles.normalizationtype], fieldnames(handles.alldata.(['Data', num2str(datachsn(i))])))
                plotthis = handles.alldata.(['Data', num2str(datachsn(i))]).(['Displaydata', handles.normalizationtype])(:,:,ddisp);
            else
                plotthis = handles.alldata.(['Data', num2str(datachsn(i))]).Displaydata(:,:,ddisp);
            end
            % Display available range
            set(handles.textdepth,'String', ['Z: ' , num2str(handles.alldata.(['Data', num2str(datachsn(i))]).Displayzpos(ddisp))]);
        else
            if ismember(['data', handles.normalizationtype], fieldnames(handles.alldata.(['Data', num2str(datachsn(i))])))
                plotthis = handles.alldata.(['Data', num2str(datachsn(i))]).(['data', handles.normalizationtype])(:,:,draw);
            else
                plotthis = handles.alldata.(['Data', num2str(datachsn(i))]).data(:,:,draw);
            end

            % Display available range
            set(handles.textdepth,'String', ['Z: ' , num2str(handles.alldata.(['Data', num2str(datachsn(i))]).zpos(draw))]);
        end

        set(handles.s, 'Zdata', plotthis);
        
    end

end

handles.redraw = 0;

guidata(hObject, handles);



% --- Executes on button press in addfunction.
function addfunction_Callback(hObject, eventdata, handles)
% hObject    handle to addfunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[myfcnfile, fcnpath] = uigetfile;

% check if cancel is pressed
if fcnpath~=0
fid = fopen([fcnpath, '/', myfcnfile]);

% 1st line = the function
% 2nd line = constants found in the function
% 3rd line = constants that will be prompted from the user
% Use matlab syntax without comments
% End the file with 'end'!!!

popupvalues = {};
myfcns = struct;
myfcns.myfcn = fgetl(fid);
myconstants = strsplit(regexprep(fgetl(fid), ' ',''),',');
myprompts = strsplit(regexprep(fgetl(fid),' ',''), ',');
excount = 1;
set(handles.operatorspopup, 'Visible', 'On');

line = '';
% Rest of the file contains pre-defined cases. Start of an example case
% marked with 'id:' before the name
while ~strcmp(line,'end')
line = fgetl(fid);
if contains(line,'id:')
    name = regexprep(line(4:end), ' ','_');
    myfcns.(name).name = name;
    popupvalues{end+1} = name;
    
    
    myfcns.(name).prompts = myprompts;
    
    while ~isempty(line) && ~strcmp(line, -1)
        line = fgetl(fid);
        tempconst = strsplit(regexprep(line, ' ',''),'=');
        tempconstcmp = tempconst{1};
        if ismember(tempconstcmp,myconstants) && ~ismember(tempconstcmp,myprompts)
            [~, idxc] = ismember(tempconstcmp,myconstants);
            % Save the line as it is (and use eval later to evaluate it)
            myfcns.(name).myconstants{idxc} = line;
        end
    end
    excount = excount + 1;
end

end


% Special case
fnams = fieldnames(myfcns);
%if contains(lower(myfcnfile), 'jarkko')
    handles.specialcase = true;
    for i = 1:length(fnams)
       if ~contains(fnams{i},'myfcn')
            
            myfcns.(fnams{i}).special = true;
       end
    end
%else
%    handles.specialcase = false;
%    for i = 1:length(fnams)-1
%        if ~contains(fnams{i},'myfcn')
%            myfcns.(fnams{i}).special = false;
%        end
%    end
%end

% Add a handles and update
set(handles.operatorspopup, 'String', popupvalues);
handles.myfcns = myfcns;
guidata(hObject, handles);

end




% --- Executes on selection change in operatorspopup.
function operatorspopup_Callback(hObject, eventdata, handles)
% hObject    handle to operatorspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




contents = cellstr(get(hObject,'String'));
chosenpopup = contents(get(hObject,'Value'));

tempchosen = find(handles.idx);

if ~isempty(tempchosen)
% Get filename. If data is extracted from other data, get parent name
if isfield(handles.alldata.(['Data', num2str(tempchosen(1))]), 'parentname')
    handles.fcnfilename = handles.alldata.(['Data', num2str(tempchosen(1))]).parentname;
else
    handles.fcnfilename = handles.alldata.(['Data', num2str(tempchosen(1))]).name;
end
handles.fcnpath = handles.alldata.(['Data', num2str(tempchosen(1))]).directory;

% Try to search MU and Fraction files in a case of the special case
if handles.specialcase == true || handles.specialcase == false
    MUs = ''; 
    fractions = '';
    operatorsdir = dir(handles.fcnpath);
    handles.fcnfilename = removeFiletype(handles.fcnfilename);
    for i = 1:length(operatorsdir)
        tempname = operatorsdir(i).name;
        if contains(lower(tempname), lower(handles.fcnfilename))
           if contains(tempname, 'MU')
                fidMU = fopen([handles.fcnpath, '\', operatorsdir(i).name]);
                MUs = fgetl(fidMU);
                fclose(fidMU);
                disp('MU file found'); 
                
           elseif contains(lower(tempname), 'fraction')
                fidFractions = fopen([handles.fcnpath, '\', operatorsdir(i).name]);
                fractions = fgetl(fidFractions);
                fclose(fidFractions);
                disp('Fractions file found');
           end
        end
    end 
    for j = 1:length(handles.myfcns.(chosenpopup{1}).prompts)
       if strcmp(handles.myfcns.(chosenpopup{1}).prompts{j}, 'N') || strcmp(handles.myfcns.(chosenpopup{1}).prompts{j}, 'U')
            switch handles.myfcns.(chosenpopup{1}).prompts{j}
                case 'N'
                    prevals{j} = fractions;
                case 'U'
                    prevals{j} = MUs;
            end
       else
           prevals{j} = '';
       end
 
    end
else
    for j = 1:length(handles.myfcns.(chosenpopup{1}).prompts)
        prevals{j} = '';
    end
end

% Prompt values
prompt = [handles.myfcns.(chosenpopup{1}).prompts, 'Data Unit'];
prevals{end+1} = 'mGy';
dlgtitle = 'Set constants';
dims = [1 80];
try
   newline = inputdlg(prompt, dlgtitle, dims, prevals);
catch
   newline = inputdlg(prompt, dlgtitle, dims);
end


% check if cancel is pressed
if ~isempty(newline)
    % Add prompted values to constants
    for i = 1:length(prompt)-1
        if ~isempty(cell2mat(newline(i)))
            handles.myfcns.(chosenpopup{1}).myconstants{end+1} = [prompt{i}, ' = ',cell2mat(newline(i))];
        % Special case
        elseif contains(prompt{i}, 'DbackCH')
            handles.myfcns.(chosenpopup{1}).myconstants{end+1} = 'DbackCH = DbackCH10x10';
        % Otherwise to 0
        else
            handles.myfcns.(chosenpopup{1}).myconstants{end+1} = [prompt{i}, ' = 0'];
        end
    end

    % Compute the parameters to memory
    for i = 1:length(handles.myfcns.(chosenpopup{1}).myconstants)
        if ~isempty(handles.myfcns.(chosenpopup{1}).myconstants{i})
            eval(handles.myfcns.(chosenpopup{1}).myconstants{i});
        end
    end

    % Compute new data
    indices = find(handles.idx);

    for i = 1:length(indices)
       try
            % Original data
            if ismember('data', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).data;
                eval(handles.myfcns.myfcn);
                myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
                handles.alldata.(['Data', num2str(indices(i))]).data = eval(myvar{1});
            elseif ismember('dataNON', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).dataNON;
                eval(handles.myfcns.myfcn);
                myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
                handles.alldata.(['Data', num2str(indices(i))]).dataNON = eval(myvar{1});
            end

            % Displayed data
            if ismember('Displaydata', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).Displaydata;
                eval(handles.myfcns.myfcn);
                myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
                handles.alldata.(['Data', num2str(indices(i))]).Displaydata = eval(myvar{1});
            elseif ismember('DisplaydataNON', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                Fdose = handles.alldata.(['Data', num2str(i)]).DisplaydataNON; 
                eval(handles.myfcns.myfcn);
                myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = eval(myvar{1});
            end

            % Interpolated data
            if ismember('interpolated', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2);
                eval(handles.myfcns.myfcn);
                myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
                handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2) = eval(myvar{1});
            end
            if ismember('interpolatedCAX', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,2);
                eval(handles.myfcns.myfcn);
                myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
                handles.alldata.(['Data', num2str(indices(i))]).interpolatedCAX(:,2) = eval(myvar{1});
                disp('Data operated successfully');
            end


       catch
            disp('2D or 3D data');
            % Original data
            try
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).data;
            catch
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).dataNON;
            end
            eval(handles.myfcns.myfcn);
            myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
            handles.alldata.(['Data', num2str(indices(i))]).data = eval(myvar{1});
            % Displayed data
            try
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).Displaydata;
            catch
                Fdose = handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON;
            end

            eval(handles.myfcns.myfcn);
            myvar = regexprep(strsplit(handles.myfcns.myfcn, '='),' ','');
            if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype,'3D')
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = eval(myvar{1});
            else
                handles.alldata.(['Data', num2str(indices(i))]).Displaydata = eval(myvar{1});
            end
        end
    end

    if ~isempty(newline(end))
        handles.alldata.(['Data', num2str(indices(i))]).DataUnit = cell2mat(newline(end));
    else
        handles.alldata.(['Data', num2str(indices(i))]).DataUnit = 'Unknown';
    end

    % Update the figure with only re-drawing chosen data
    handles.redraw = 1;
    %testsave = handles.alldata;
    %save('alladata.mat', 'testsave');

    guidata(hObject, handles);

    try
        listboxChooseData_Callback(handles.listboxChooseData, [], handles);
    catch
       disp('Error while operating data. Figure not updated'); 
    end

    % Hints: contents = cellstr(get(hObject,'String')) returns operatorspopup contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from operatorspopup
end
else
    disp('No data chosen...');
end
% --- Executes during object creation, after setting all properties.
function operatorspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to operatorspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rotateX.
function rotateX_Callback(hObject, eventdata, handles)
% hObject    handle to rotateX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indices = find(handles.idx);

for i = 1:length(indices)
    if  ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'profile') && ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'pdd')
            % Swap position axes
            if strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
                
                
                temp = handles.alldata.(['Data', num2str(indices(i))]).ypos;
                handles.alldata.(['Data', num2str(indices(i))]).ypos = handles.alldata.(['Data', num2str(indices(i))]).zpos;
                handles.alldata.(['Data', num2str(indices(i))]).zpos = temp;
                temp = handles.alldata.(['Data', num2str(indices(i))]).Displayypos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayypos = handles.alldata.(['Data', num2str(indices(i))]).Displayzpos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayzpos = temp;

                % Rotate data along X by 90 degrees
                handles.alldata.(['Data', num2str(indices(i))]).data = permute(handles.alldata.(['Data', num2str(indices(i))]).data, [3,2,1]);
                handles.alldata.(['Data', num2str(indices(i))]).Displaydata = permute(handles.alldata.(['Data', num2str(indices(i))]).Displaydata, [3,2,1]);


                % Rotate error data if existing

                if ismember('error', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).error = permute(handles.alldata.(['Data', num2str(indices(i))]).error, [3,2,1]);
                end
                
                handles.alldata.(['Data',num2str(indices(i))]).Depthdisplay = round(length(handles.alldata.(['Data',num2str(indices(i))]).Displaydata(1,1,:))/2);
                handles.alldata.(['Data',num2str(indices(i))]).Depthraw = round(length(handles.alldata.(['Data',num2str(indices(i))]).data(1,1,:))/2);
            end
    end
end


% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in rotating data. Figure not updated'); 
end







% --- Executes on button press in rotateY.
function rotateY_Callback(hObject, eventdata, handles)
% hObject    handle to rotateY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indices = find(handles.idx);

for i = 1:length(indices)
    if  ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'profile') && ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'pdd')
            % Swap position axes
            if strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
                

                
                temp = handles.alldata.(['Data', num2str(indices(i))]).xpos;
                handles.alldata.(['Data', num2str(indices(i))]).xpos = handles.alldata.(['Data', num2str(indices(i))]).zpos;
                handles.alldata.(['Data', num2str(indices(i))]).zpos = temp;
                temp = handles.alldata.(['Data', num2str(indices(i))]).Displayxpos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayxpos = handles.alldata.(['Data', num2str(indices(i))]).Displayzpos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayzpos = temp;

                % Rotate data along Y by 90 degrees
                handles.alldata.(['Data', num2str(indices(i))]).data = permute(handles.alldata.(['Data', num2str(indices(i))]).data, [1,3,2]);
                handles.alldata.(['Data', num2str(indices(i))]).Displaydata = permute(handles.alldata.(['Data', num2str(indices(i))]).Displaydata, [1,3,2]);


                % Rotate error data if existing

                if ismember('error', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).error = permute(handles.alldata.(['Data', num2str(indices(i))]).error, [1,3,2]);
                end
                handles.alldata.(['Data',num2str(indices(i))]).Depthdisplay = round(length(handles.alldata.(['Data',num2str(indices(i))]).Displaydata(1,1,:))/2);
                handles.alldata.(['Data',num2str(indices(i))]).Depthraw = round(length(handles.alldata.(['Data',num2str(indices(i))]).data(1,1,:))/2);
            end
    end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in rotating data. Figure not updated'); 
end




% --- Executes on button press in rotateZ.
function rotateZ_Callback(hObject, eventdata, handles)
% hObject    handle to rotateZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indices = find(handles.idx);

for i = 1:length(indices)
    if  ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'profile') && ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'pdd')
            % Swap position axes
            if strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
                if ~ismember('zangle', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).zangle = 90; 
                elseif handles.alldata.(['Data', num2str(indices(i))]).zangle >= 270
                   handles.alldata.(['Data', num2str(indices(i))]).zangle = 0; 
                else
                    handles.alldata.(['Data', num2str(indices(i))]).zangle = handles.alldata.(['Data', num2str(indices(i))]).zangle + 90;
                end
%                 
%                 if handles.alldata.(['Data', num2str(indices(i))]).zangle == 90
%                     handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
%                 elseif handles.alldata.(['Data', num2str(indices(i))]).zangle == 180
%                     handles.alldata.(['Data', num2str(indices(i))]).ypos = flip(handles.alldata.(['Data', num2str(indices(i))]).ypos);
%                 elseif handles.alldata.(['Data', num2str(indices(i))]).zangle == 270
%                     handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
%                 else
%                     handles.alldata.(['Data', num2str(indices(i))]).ypos = flip(handles.alldata.(['Data', num2str(indices(i))]).ypos);
%                 end
                
                temp = handles.alldata.(['Data', num2str(indices(i))]).ypos;
                handles.alldata.(['Data', num2str(indices(i))]).ypos = handles.alldata.(['Data', num2str(indices(i))]).xpos;
                handles.alldata.(['Data', num2str(indices(i))]).xpos = temp;
                temp = handles.alldata.(['Data', num2str(indices(i))]).Displayypos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayypos = handles.alldata.(['Data', num2str(indices(i))]).Displayxpos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayxpos = temp;

                % Rotate data along Z by 90 degrees
                handles.alldata.(['Data', num2str(indices(i))]).data = rot90(handles.alldata.(['Data', num2str(indices(i))]).data);
                handles.alldata.(['Data', num2str(indices(i))]).Displaydata = rot90(handles.alldata.(['Data', num2str(indices(i))]).Displaydata);

                % Rotate error data if existing

                if ismember('error', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).error = rot90(handles.alldata.(['Data', num2str(indices(i))]).error);
                end
                
            elseif strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D')
                try
                    
                if ~ismember('zangle', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).zangle = 90; 
                elseif handles.alldata.(['Data', num2str(indices(i))]).zangle >= 270
                   handles.alldata.(['Data', num2str(indices(i))]).zangle = 0; 
                else
                    handles.alldata.(['Data', num2str(indices(i))]).zangle = handles.alldata.(['Data', num2str(indices(i))]).zangle + 90;
                end
                
%                 if handles.alldata.(['Data', num2str(indices(i))]).zangle == 90
%                     handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
%                     handles.alldata.(['Data', num2str(indices(i))]).Displayxpos = flip(handles.alldata.(['Data', num2str(indices(i))]).Displayxpos);
%                 elseif handles.alldata.(['Data', num2str(indices(i))]).zangle == 180
%                     handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
%                     handles.alldata.(['Data', num2str(indices(i))]).Displayxpos = flip(handles.alldata.(['Data', num2str(indices(i))]).Displayxpos);
%                 elseif handles.alldata.(['Data', num2str(indices(i))]).zangle == 270
%                     handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
%                     handles.alldata.(['Data', num2str(indices(i))]).Displayxpos = flip(handles.alldata.(['Data', num2str(indices(i))]).Displayxpos);
%                 else
%                     handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
%                     handles.alldata.(['Data', num2str(indices(i))]).Displayypos = flip(handles.alldata.(['Data', num2str(indices(i))]).Displayypos);
%                 end
                %handles.alldata.(['Data', num2str(indices(i))]).xpos = flip(handles.alldata.(['Data', num2str(indices(i))]).xpos);
               % handles.alldata.(['Data', num2str(indices(i))]).Displayypos = flip(handles.alldata.(['Data', num2str(indices(i))]).Displayypos);    
                    
                temp = handles.alldata.(['Data', num2str(indices(i))]).ypos;
                handles.alldata.(['Data', num2str(indices(i))]).ypos = handles.alldata.(['Data', num2str(indices(i))]).xpos;
                handles.alldata.(['Data', num2str(indices(i))]).xpos = temp;
                temp = handles.alldata.(['Data', num2str(indices(i))]).Displayypos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayypos = handles.alldata.(['Data', num2str(indices(i))]).Displayxpos;
                handles.alldata.(['Data', num2str(indices(i))]).Displayxpos = temp;

                % Rotate data along Z by 90 degrees
                handles.alldata.(['Data', num2str(indices(i))]).dataCAX = rot90(handles.alldata.(['Data', num2str(indices(i))]).dataCAX);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX = rot90(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX);
                handles.alldata.(['Data', num2str(indices(i))]).dataMAX = rot90(handles.alldata.(['Data', num2str(indices(i))]).dataMAX);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX = rot90(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX);
                handles.alldata.(['Data', num2str(indices(i))]).dataNON = rot90(handles.alldata.(['Data', num2str(indices(i))]).dataNON);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = rot90(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON);
                try
                    handles.alldata.(['Data', num2str(indices(i))]).dataMAN = rot90(handles.alldata.(['Data', num2str(indices(i))]).dataMAN);
                    handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN = rot90(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN);
                catch
                end
                try
                    handles.alldata.(['Data', num2str(indices(i))]).dataMEAN = rot90(handles.alldata.(['Data', num2str(indices(i))]).dataMEAN);
                    handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN = rot90(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN);
                catch
                end
                % Rotate error data if existing

                if ismember('error', fieldnames(handles.alldata.(['Data', num2str(indices(i))])))
                    handles.alldata.(['Data', num2str(indices(i))]).error = rot90(handles.alldata.(['Data', num2str(indices(i))]).error);
                end
                
                catch
                   disp('Error in 2D rotation along Z-axis'); 
                end
            end
    end
end



% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
  disp('Error in rotating data. Figure not updated'); 
end


% % Flip the axes also in the figure
%  if handles.alldata.(['Data', num2str(indices(i))]).zangle == 90
%      set ( gca, 'ydir', 'reverse' )
%  elseif handles.alldata.(['Data', num2str(indices(i))]).zangle == 180
%      set ( gca, 'xdir', 'reverse' )
%  elseif handles.alldata.(['Data', num2str(indices(i))]).zangle == 270
%      set ( gca, 'ydir', 'normal' )
%  else
%      set ( gca, 'ydir', 'normal' )
%      set ( gca, 'xdir', 'normal' )
%  end



% --- Executes on button press in bruteforcedta.
function bruteforcedta_Callback(hObject, eventdata, handles)
% hObject    handle to bruteforcedta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DTAbruteforce = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of bruteforcedta
guidata(hObject, handles);


% --- Executes on button press in analyticaldta.
function analyticaldta_Callback(hObject, eventdata, handles)
% hObject    handle to analyticaldta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DTAanalytical = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of analyticaldta
guidata(hObject, handles);


% --- Executes on button press in removedatabutton.
function removedatabutton_Callback(hObject, eventdata, handles)

% Go through chosen data and remove it from the list
% Need to be removed from alldata structure, handles.idx, handles.origlist,
% handles.listboxarray, handles.order, reduce datacount, listbox object

tempstruct = struct;

for i = 1:length(handles.listboxArray)
   if contains(handles.listboxArray{i},'Chosen')
       handles.alldata.(['Data', num2str(i)]) = [];
       handles.origlist{i} = [];
       handles.listboxArray{i} = [];
       handles.datacount = handles.datacount - 1;
       % Update the listbox
   end
end

tempfields = fieldnames(handles.alldata);
delcount = 0;
for i = 1:length(tempfields)
    if ~isempty(handles.alldata.(tempfields{i}))
        % Collect the non-empty fields
        tempcount = str2double(tempfields{i}(5:end)) - delcount;
        tempstruct.(['Data', num2str(tempcount)]) = handles.alldata.(tempfields{i});
    else
        % If structure field was deleted, decrease the Data count from the
        % fieldnames
        handles.idx(i-delcount) = [];
        handles.order(i-delcount) = [];
        delcount = delcount + 1;
    end
end
handles.alldata = struct;
handles.alldata = tempstruct;

handles.listboxArray = handles.listboxArray(~cellfun('isempty',handles.listboxArray));
handles.origlist = handles.origlist(~cellfun('isempty',handles.origlist));
% Update the figure with only re-drawing chosen data
handles.redraw = 1;
set(handles.listboxChooseData, 'String', handles.listboxArray);
set(handles.listboxChooseData, 'Value', 1);
guidata(hObject, handles);
%testdata = handles.alldata;
%save('alladata.mat', 'testdata');

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in rotating data. Figure not updated'); 
end

% hObject    handle to removedatabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in mirror2DX.
function mirror2DX_Callback(hObject, eventdata, handles)
% hObject    handle to mirror2DX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



indices = find(handles.idx);
for i = 1:length(indices)
    if strcmp('2D', handles.alldata.(['Data', num2str(indices(i))]).datatype)
            % Mirror X
            handles.alldata.(['Data', num2str(indices(i))]).dataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataCAX,2);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX,2);
            
            handles.alldata.(['Data', num2str(indices(i))]).dataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAX,2);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX,2);
            
            handles.alldata.(['Data', num2str(indices(i))]).dataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).dataNON,2);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON,2);
            
            try
                handles.alldata.(['Data', num2str(indices(i))]).dataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAN,2);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN,2);
            catch
            end
            
            try
                handles.alldata.(['Data', num2str(indices(i))]).dataMENA = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMEAN,2);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN,2);
            catch
            end
    elseif strcmp('3D', handles.alldata.(['Data', num2str(indices(i))]).datatype)
        % Mirror X
        handles.alldata.(['Data', num2str(indices(i))]).data = flip(handles.alldata.(['Data', num2str(indices(i))]).data,2);
        handles.alldata.(['Data', num2str(indices(i))]).Displaydata = flip(handles.alldata.(['Data', num2str(indices(i))]).Displaydata,2);
        
    end
end
% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end



% --- Executes on button press in mirror2DY.
function mirror2DY_Callback(hObject, eventdata, handles)
% hObject    handle to mirror2DY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indices = find(handles.idx);
for i = 1:length(indices)
    if strcmp('2D', handles.alldata.(['Data', num2str(indices(i))]).datatype)
            % Mirror Y
            handles.alldata.(['Data', num2str(indices(i))]).dataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataCAX,1);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX,1);
            
            handles.alldata.(['Data', num2str(indices(i))]).dataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAX,1);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX,1);
            
            handles.alldata.(['Data', num2str(indices(i))]).dataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).dataNON,1);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON,1);
            
            try
                handles.alldata.(['Data', num2str(indices(i))]).dataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAN,1);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN,1);
            catch
            end
            
            try
                handles.alldata.(['Data', num2str(indices(i))]).dataMENA = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMEAN,1);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN,1);
            catch
            end
    elseif strcmp('3D', handles.alldata.(['Data', num2str(indices(i))]).datatype)
        % Mirror X
        handles.alldata.(['Data', num2str(indices(i))]).data = flip(handles.alldata.(['Data', num2str(indices(i))]).data,1);
        handles.alldata.(['Data', num2str(indices(i))]).Displaydata = flip(handles.alldata.(['Data', num2str(indices(i))]).Displaydata,1);
        
    end
end
% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end




% --- Executes on button press in mirror2DZ.
function mirror2DZ_Callback(hObject, eventdata, handles)
% hObject    handle to mirror2DZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


indices = find(handles.idx);
for i = 1:length(indices)
    if strcmp('2D', handles.alldata.(['Data', num2str(indices(i))]).datatype) && strcmp('XZ', handles.alldata.(['Data', num2str(indices(i))]).Plane)
            % Mirror Z
            handles.alldata.(['Data', num2str(indices(i))]).dataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataCAX,2);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataCAX,2);
            
            handles.alldata.(['Data', num2str(indices(i))]).dataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAX,2);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAX,2);
            
            handles.alldata.(['Data', num2str(indices(i))]).dataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).dataNON,2);
            handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataNON,2);
            
            try
                handles.alldata.(['Data', num2str(indices(i))]).dataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMAN,2);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMAN,2);
            catch
            end
            
            try
                handles.alldata.(['Data', num2str(indices(i))]).dataMENA = flip(handles.alldata.(['Data', num2str(indices(i))]).dataMEAN,2);
                handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN = flip(handles.alldata.(['Data', num2str(indices(i))]).DisplaydataMEAN,2);
            catch
            end
    elseif strcmp('3D', handles.alldata.(['Data', num2str(indices(i))]).datatype)
        % Mirror X
        handles.alldata.(['Data', num2str(indices(i))]).data = flip(handles.alldata.(['Data', num2str(indices(i))]).data,3);
        handles.alldata.(['Data', num2str(indices(i))]).Displaydata = flip(handles.alldata.(['Data', num2str(indices(i))]).Displaydata,3);
        
    end
end
% Update the figure with only re-drawing chosen data
handles.redraw = 1;

guidata(hObject, handles);

try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end



% --- Executes on button press in show2dgamma.
function show2dgamma_Callback(hObject, eventdata, handles)
% hObject    handle to show2dgamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

donothing = 0;
indx = 1;
indice = find(handles.idx);
useprompt = false;
try
    fnames = fieldnames(handles.alldata.(['Data', num2str(indice)]));
catch
    disp('Error, too many datasets chosen! Choose only one reference dataset at a time!');
    return;
end

try
    gammaidx = find(contains(fnames, 'gamma'));
    if handles.alldata.(['Data', num2str(indice)]).chosengam == 0 || ...
            strcmp(handles.alldata.(['Data', num2str(indice)]).currentgaman, 'gamma')
        useprompt = true;
    else
        indx = handles.alldata.(['Data', num2str(indice)]).chosengam;
        useprompt = false;
    end
catch
   disp('No 2D gamma, DTA or DD analysis existing for chosen data.'); 
end


if length(gammaidx) > 1 && useprompt == true
   % Prompt wanted analysis
    list = {};
    for i = 1:length(gammaidx)
        list{end+1} = fnames{gammaidx(i)};
    end
    [indx,~] = listdlg('PromptString',{'Multiple results found.',...
    'Only one result can be selected at a time.',''},...
    'SelectionMode','single','ListString',list);
    
    % check if cancel is pressed
    if isempty(indx)
        donothing = 1;
    else
        donothing = 0;
    end
end

if donothing == 0
    tempdata = handles.alldata.(['Data', num2str(indice)]).(fnames{gammaidx(indx)});
    if ismember('gamma', fieldnames(tempdata))
        axes(handles.axes1);
        cla;
        surf(tempdata.GAMposx(:,1), tempdata.GAMposy(:,1), tempdata.gamma);
        view(0,90);
        axis([tempdata.GAMposx(1) tempdata.GAMposx(end) tempdata.GAMposy(1) tempdata.GAMposy(end)]);
        legend(['gamma map (', num2str(tempdata.gammapassperc), '%)']);
        colorbar;
    end
    handles.alldata.(['Data', num2str(indice)]).chosengam = indx;
    handles.alldata.(['Data', num2str(indice)]).chosengamname = fnames{gammaidx(indx)};
    handles.alldata.(['Data', num2str(indice)]).currentgaman = 'gamma';
end

guidata(hObject, handles);


% --- Executes on button press in show2dDTA.
function show2dDTA_Callback(hObject, eventdata, handles)
% hObject    handle to show2dDTA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

donothing = 0;
indx = 1;
indice = find(handles.idx);
useprompt = false;

try
    fnames = fieldnames(handles.alldata.(['Data', num2str(indice)]));
catch
    disp('Error, too many datasets chosen! Choose only one reference dataset at a time!');
    return;
end

try
    gammaidx = find(contains(fnames, 'gamma'));
    if handles.alldata.(['Data', num2str(indice)]).chosengam == 0 || ...
            strcmp(handles.alldata.(['Data', num2str(indice)]).currentgaman, 'DTA')
        useprompt = true;
    else
        indx = handles.alldata.(['Data', num2str(indice)]).chosengam;
        useprompt = false;
    end
catch
   disp('No 2D gamma, DTA or DD analysis existing for chosen data.'); 
end



if length(gammaidx) > 1 && useprompt == true
   % Prompt wanted analysis
    list = {};
    for i = 1:length(gammaidx)
        list{end+1} = fnames{gammaidx(i)};
    end
    [indx,~] = listdlg('PromptString',{'Multiple results found.',...
    'Only one result can be selected at a time.',''},...
    'SelectionMode','single','ListString',list);
    
    % check if cancel is pressed
    if isempty(indx)
        donothing = 1;
    else
        donothing = 0;
    end
end
if donothing == 0
    tempdata = handles.alldata.(['Data', num2str(indice)]).(fnames{gammaidx(indx)});
    if ismember('dta', fieldnames(tempdata))
        axes(handles.axes1);
        cla;

        surf(tempdata.GAMposx(:,1), tempdata.GAMposy(:,1), tempdata.dta);

        view(0,90);
        axis([tempdata.GAMposx(1) tempdata.GAMposx(end) tempdata.GAMposy(1) tempdata.GAMposy(end)]);
        legend(['DTA (accuracy == resolution) (', num2str(tempdata.dtapassperc), '%)']);
        colorbar;
    end
    handles.alldata.(['Data', num2str(indice)]).chosengam = indx;
    handles.alldata.(['Data', num2str(indice)]).chosengamname = fnames{gammaidx(indx)};
    handles.alldata.(['Data', num2str(indice)]).currentgaman = 'DTA';
end

guidata(hObject, handles);



% --- Executes on button press in show2dDD.
function show2dDD_Callback(hObject, eventdata, handles)
% hObject    handle to show2dDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

donothing = 0;
indx = 1;
indice = find(handles.idx);
useprompt = false;

try
    fnames = fieldnames(handles.alldata.(['Data', num2str(indice)]));
catch
    disp('Error, too many datasets chosen! Choose only one reference dataset at a time!');
    return;
end

try
    gammaidx = find(contains(fnames, 'gamma'));
    if handles.alldata.(['Data', num2str(indice)]).chosengam == 0 || ...
            strcmp(handles.alldata.(['Data', num2str(indice)]).currentgaman, 'DD')
        useprompt = true;
    else
        indx = handles.alldata.(['Data', num2str(indice)]).chosengam;
        useprompt = false;
    end
catch
   disp('No 2D gamma, DTA or DD analysis existing for chosen data.'); 
end



if length(gammaidx) > 1 && useprompt == true
   % Prompt wanted analysis
    list = {};
    for i = 1:length(gammaidx)
        list{end+1} = fnames{gammaidx(i)};
    end
    [indx,~] = listdlg('PromptString',{'Multiple results found.',...
    'Only one result can be selected at a time.',''},...
    'SelectionMode','single','ListString',list);
    
    % check if cancel is pressed
    if isempty(indx)
        donothing = 1;
    else
        donothing = 0;
    end
end
if donothing == 0
    tempdata = handles.alldata.(['Data', num2str(indice)]).(fnames{gammaidx(indx)});
    if ismember('absoluteDD', fieldnames(tempdata))
        axes(handles.axes1);
        cla;
        tempdata = handles.alldata.(['Data', num2str(indice)]).(fnames{gammaidx(indx)});
        surf(tempdata.GAMposx(:,1), tempdata.GAMposy(:,1), tempdata.absoluteDD);

        view(0,90);
        axis([tempdata.GAMposx(1) tempdata.GAMposx(end) tempdata.GAMposy(1) tempdata.GAMposy(end)]);
        legend(['absolute Dose difference (', num2str(tempdata.DDpassperc), '%)']);
        colorbar;
        
    end
    handles.alldata.(['Data', num2str(indice)]).chosengam = indx;
    handles.alldata.(['Data', num2str(indice)]).chosengamname = fnames{gammaidx(indx)};
    handles.alldata.(['Data', num2str(indice)]).currentgaman = 'DD';
end

guidata(hObject, handles);




% --- Executes on button press in savematfile.
function savematfile_Callback(hObject, eventdata, handles)
% hObject    handle to savematfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tempdata = handles.alldata;
save('alldata.mat', 'tempdata');

guidata(hObject, handles);




% --- Executes on button press in quickprofile.
function quickprofile_Callback(hObject, eventdata, handles)
% hObject    handle to quickprofile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);

% Create line so it will surface the uistack
handles.myline = line(0:0.1,0:0.1);
handles.fixedline = line(0:0.1,0:0.1);
%handles.legendlist{end+1} = 
legend(handles.legendlist, 'Interpreter', 'None', 'Location','northeast');

% Draw line (use the created l)
handles.myline = drawline('Color', 'r');
handles.myline.Label = 'myline';
set(handles.myline, 'LabelVisible', 'Off');


% If reference is fixed, draw a fixed line to the figure
if handles.reffixed == 1
    handles.fixedline = drawline('Position', [handles.myline.Position(1,1), handles.myline.Position(1,2); handles.myline.Position(2,1), handles.myline.Position(2,2)], 'Color', 'g');
    bringToFront(handles.fixedline);
    handles.fixedline.Label = 'fixed';
    set(handles.fixedline, 'LabelVisible', 'Off');
    disp('Fixed referenceline created');
end

% Make sure that the line show in front (rotating the figure seems to help)
bringToFront(handles.myline);

% Get the ref idx
[~,handles.refidx] = find(handles.order == 1);
% Find chosen data
handles.currentidx = find(handles.idx);

figure(2);
cla;
grid on
hold on
diagonal = [];


for i = 1:length(handles.currentidx)
    gammaon = 0;
    switch handles.alldata.(['Data', num2str(handles.currentidx(i))]).datatype
        case '2D'
            % Create a line mask. Matrix which is bound by the line ends in current axes (minimum dim size == 1)
            % Find incides from current file
            if strcmp(handles.alldata.(['Data', num2str(handles.currentidx(i))]).Plane, 'XY')
                [~,x1] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).xpos - handles.myline.Position(1,1)));
                [~,x2] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).xpos - handles.myline.Position(2,1)));
                [~,y1] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).ypos - handles.myline.Position(1,2)));
                [~,y2] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).ypos - handles.myline.Position(2,2)));
            elseif strcmp(handles.alldata.(['Data', num2str(handles.currentidx(i))]).Plane, 'YZ')
                [~,x1] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).zpos - handles.myline.Position(1,1)));
                [~,x2] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).zpos - handles.myline.Position(2,1)));
                [~,y1] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).ypos - handles.myline.Position(1,2)));
                [~,y2] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).ypos - handles.myline.Position(2,2)));
            else
                [~,x1] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).zpos - handles.myline.Position(1,1)));
                [~,x2] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).zpos - handles.myline.Position(2,1)));
                [~,y1] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).xpos - handles.myline.Position(1,2)));
                [~,y2] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).xpos - handles.myline.Position(2,2)));
            end
            % Print DTA/DD/GAMMA data if chosen
            if isfield(handles.alldata.(['Data', num2str(handles.currentidx(i))]), 'chosengam') && handles.alldata.(['Data', num2str(handles.currentidx(i))]).chosengam ~= 0
                gamname = handles.alldata.(['Data', num2str(handles.currentidx(i))]).chosengamname;
                [~,x11] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).GAMposx - handles.myline.Position(1,1)));
                [~,x22] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).GAMposx - handles.myline.Position(2,1)));
                [~,y11] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).GAMposy - handles.myline.Position(1,2)));
                [~,y22] = min(abs(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).GAMposy - handles.myline.Position(2,2)));
                gammaon = 1;
            end
            % Line ends are always in the corners of the datamatrix mask (except horizontal or vertical cases) -> get
            % diagonal values.
            
            % Define the initial data bound by the created box
            if isfield(handles.alldata.(['Data', num2str(handles.currentidx(i))]), 'chosengam') && handles.alldata.(['Data', num2str(handles.currentidx(i))]).chosengam ~= 0
                switch handles.alldata.(['Data', num2str(handles.currentidx(i))]).currentgaman
                    case 'DTA'
                        datamask2 = handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).dta(min([y11,y22]):max([y11,y22]), min([x11,x22]):max([x11,x22]));
                        handles.minval = min(min(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).dta));
                    case 'DD'
                        datamask2 = handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).absoluteDD(min([y11,y22]):max([y11,y22]), min([x11,x22]):max([x11,x22]));
                        handles.minval = min(min(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).absoluteDD));
                    case 'gamma'
                        datamask2 = handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).gamma(min([y11,y22]):max([y11,y22]), min([x11,x22]):max([x11,x22]));
                        handles.minval = min(min(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(gamname).gamma));
                end
                if x11 > x22
                 datamask2 = flip(datamask,2);
                 end
                 if y11 > y22
                    datamask2 = flip(datamask2,1);
                 end
                 
                 % Get size
                 [ysiz2, xsiz2] = size(datamask2);
                 
                 % Check if diagonal can be computed (dimsize > 1),
                 % otherwise take the horizontal/vertical vector
                 if ysiz2 > 1 && xsiz2 > 1

                     [X, Y] = meshgrid(1:xsiz2, 1:ysiz2);
                     if xsiz2 <= ysiz2
                         xax = linspace(1,xsiz2,ysiz2);
                         [Xq, Yq] = meshgrid(xax,1:ysiz2);
                         datainterpmask = interp2(X,Y,datamask2,Xq,Yq);
                         diagonal2 = diag(datainterpmask);
                     elseif ysiz2 <= xsiz2
                         yax = linspace(1,ysiz2,xsiz2);
                         [Xq, Yq] = meshgrid(1:xsiz2, yax);
                         datainterpmask = interp2(X,Y,datamask2,Xq,Yq);
                         diagonal2 = diag(datainterpmask);
                     end     

                 elseif ysiz < 2
                     diagonal2 = datamask2(1, 1:xsiz2);
                 elseif xsiz < 2
                     diagonal2 = datamask2(1:ysiz2, 1);
                 end

                 xposvect2 = linspace(0,sqrt((handles.myline.Position(1,2) - handles.myline.Position(2,2))^2 + (handles.myline.Position(2,1) - handles.myline.Position(1,1))^2),length(diagonal2));
            end
            datamask = handles.alldata.(['Data', num2str(handles.currentidx(i))]).(['data', handles.normalizationtype])(min([y1,y2]):max([y1,y2]), min([x1,x2]):max([x1,x2]));
                handles.minval = min(min(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(['data', handles.normalizationtype])));    
            % Flip data so that the line defines a diagonal from [0,0] to
            % [end, end]
            
            % Fix axes so the diagonal is taken correctly
             if x1 > x2
                 datamask = flip(datamask,2);
             end
             if y1 > y2
                datamask = flip(datamask,1);
             end
             % Get size
             [ysiz, xsiz] = size(datamask);

             % Check if diagonal can be computed (dimsize > 1),
             % otherwise take the horizontal/vertical vector
             if ysiz > 1 && xsiz > 1

                 [X, Y] = meshgrid(1:xsiz, 1:ysiz);
                 if xsiz <= ysiz
                     xax = linspace(1,xsiz,ysiz);
                     [Xq, Yq] = meshgrid(xax,1:ysiz);
                     datainterpmask = interp2(X,Y,datamask,Xq,Yq);
                     diagonal = diag(datainterpmask);
                 elseif ysiz <= xsiz
                     yax = linspace(1,ysiz,xsiz);
                     [Xq, Yq] = meshgrid(1:xsiz, yax);
                     datainterpmask = interp2(X,Y,datamask,Xq,Yq);
                     diagonal = diag(datainterpmask);
                 end     
            
             elseif ysiz < 2
                 diagonal = datamask(1, 1:xsiz);
             elseif xsiz < 2
                 diagonal = datamask(1:ysiz, 1);
             end
             
             xposvect = linspace(0,sqrt((handles.myline.Position(1,2) - handles.myline.Position(2,2))^2 + (handles.myline.Position(2,1) - handles.myline.Position(1,1))^2),length(diagonal));
             
        case '3D'
            %plot(handles.alldata.(['Data', num2str(handles.currentidx(i))]).(['data', handles.normalizationtype]));
            % Take the current Z depth and do the same as in 2D. 
            disp('Profiles existing only for 2D data. Extract 2D plane first from 3D!');
            return
    end
    
    
    
     handles.lineplot(i) = plot(xposvect,diagonal, 'LineWidth', 2);
     if gammaon == 1
         handles.lineplot2(i) = plot(xposvect2,diagonal2, 'LineWidth', 2, 'Color', 'r');
         yyaxis right
     end
     
     axis([0 inf handles.minval inf]);
     if isfield(handles.alldata.(['Data', num2str(handles.currentidx(i))]), 'chosengam') && handles.alldata.(['Data', num2str(handles.currentidx(i))]).chosengam ~= 0
        ylabel(handles.alldata.(['Data', num2str(handles.currentidx(i))]).currentgaman);
     elseif strcmp(handles.normalizationtype, 'NON') 
        ylabel('Dose');
     else
        ylabel('Value [normalized]'); 
     end
     xlabel('Pos [cm]');
    
end

% Plot the line
%posvectstarts = [min([y1,y2]), min([x1,x2])];
%posvectcomponents = [posvectstarts(1):1/leny:leny];
hold off
if gammaon == 0
    title('blue = ref');
else
    title(['blue = ref data, red = ', handles.alldata.(['Data', num2str(handles.currentidx(i))]).currentgaman]);
end

% Add listenerclc
if handles.reffixed == 1
    addlistener([handles.myline, handles.fixedline],'MovingROI',@(src,evt) allevents(src,evt,handles));
else
    addlistener(handles.myline,'MovingROI',@(src,evt) allevents(src,evt,handles));
end
% Create an isodose structure
% Get the dose value
% answer = inputdlg('Isodosevalue [depends on chosen normalization]');
% % check if cancel is pressed
% if ~isempty(answer)
    
%     refvalue = str2num(cell2mat(answer));
%     
% 
% 
%     indices = find(handles.idx);
% 
%     % Go through all chosen data
%     for i = 1:length(indices)
%         if strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
%             [sizey, sizex, sizez] = size(handles.alldata.(['Data', num2str(indices(i))]).data);
%             % Create an mask image with 1s for indices abose the ref value
%             [iy, ix, iz] = find(handles.alldata.(['Data', num2str(indices(i))]).data < refvalue);
%             % Change the indices to subscripts
%             subidx = sub2ind([sizey, sizex, sizez], iy, ix, iz);
%             mask = ones(sizey, sizex, sizez);
%         else
%             [sizey, sizex] = size(handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype]));
%             % Create an mask image with 1s for indices abose the ref value
%             [iy, ix] = find(handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype]) < refvalue);
%             % Change the indices to subscripts
%             subidx = sub2ind([sizey, sizex], iy, ix);
%             mask = ones(sizey, sizex);
%         end
%         
%         
%         isodosestruct = handles.alldata.(['Data', num2str(indices(i))]).(['data', handles.normalizationtype]);
%         
%         % Get the masked data, set other values to NaN
%         isodosestruct(subidx) = NaN;
%         % Binary mask image to find the boundaries
%         
%         mask(subidx) = 0;
%         % Find outline and create a line structure (smoothed for visualization) from the isodose
%         edgedata = smooth(edge(mask));
%              
%         % Create the new isodose structure and add it to the listbox
%         newisostruct = struct;
%         newisostruct.dataNON = isodosestruct;
%         newisostruct.edge = edgedata;
%         newisostruct.name = [num2str(refvalue), '_' , handles.alldata.(['Data', num2str(indices(i))]).name];
%         newisostruct.datatype = 'Isodose';
%         newisostruct.Plane = handles.alldata.(['Data', num2str(indices(i))]).Plane;
%         if contains(newisostruct.Plane, 'X') && contains(newisostruct.Plane, 'Y') && contains(newisostruct.Plane, 'Z') 
%             newisostruct.xpos = handles.alldata.(['Data', num2str(indices(i))]).xpos; 
%             newisostruct.ypos = handles.alldata.(['Data', num2str(indices(i))]).ypos;
%             newisostruct.zpos = handles.alldata.(['Data', num2str(indices(i))]).zpos;
%         elseif contains(newisostruct.Plane, 'X') && contains(newisostruct.Plane, 'Y')
%             newisostruct.xpos = handles.alldata.(['Data', num2str(indices(i))]).xpos; 
%             newisostruct.ypos = handles.alldata.(['Data', num2str(indices(i))]).ypos; 
%         elseif contains(newisostruct.Plane, 'X') && contains(newisostruct.Plane, 'Z')
%             newisostruct.xpos = handles.alldata.(['Data', num2str(indices(i))]).xpos; 
%             newisostruct.zpos = handles.alldata.(['Data', num2str(indices(i))]).zpos; 
%         else
%             newisostruct.zpos = handles.alldata.(['Data', num2str(indices(i))]).zpos; 
%             newisostruct.ypos = handles.alldata.(['Data', num2str(indices(i))]).ypos; 
%         end
%             
%         newisostruct.DataUnit = handles.alldata.(['Data', num2str(indices(i))]).DataUnit; 
%                 
%         templist = ['Isodose', num2str(handles.isodosecount), '_', newisostruct.name];
%         handles.isodosecount = handles.isodosecount + 1;
%         %templist = strrep(templist, 'Plane', '');
%         temporig = templist;
%         % Add name to listbox array
%         handles.listboxArray{end+1} = strcat(templist);
%         handles.origlist{end+1} = strcat(temporig);
%         
%        handles.idx(end+1) = 0;
%        handles.alldata.(['Data', num2str(length(handles.idx))]) = newisostruct;
%        handles.datacount = handles.datacount + 1;
%        handles.order(end+1) = 0;
%     end
%     %alladata = handles.alldata;
%     %save('alladata.mat', 'alladata');
%     set(handles.listboxChooseData, 'String', handles.listboxArray);
%     
% end


% Update the figure with only re-drawing chosen data
% handles.redraw = 1;
% 
 guidata(hObject, handles);
% 
% try
%     listboxChooseData_Callback(handles.listboxChooseData, [], handles);
% catch
%    disp('Error in axis invertion. Figure not updated'); 
% end


function gammarefdoseEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gammarefdoseEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.gammaMinimumdose = str2double(get(hObject,'String'));

guidata(hObject, handles);




% Hints: get(hObject,'String') returns contents of gammarefdoseEdit as text
%        str2double(get(hObject,'String')) returns contents of gammarefdoseEdit as a double


% --- Executes during object creation, after setting all properties.
function gammarefdoseEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gammarefdoseEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function allevents(scr,evt, myHandles)
    fixval = get(myHandles.checkboxfixref, 'Value');
    % If reference is fixed, draw a fixed line to the figure
    
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
           % Compute diagonal and update figure. Do not compute, nor update the reference data if chosen so
           if strcmp(scr.Label, 'myline') 
                myHandles.myline.Position = evt.CurrentPosition;
           else
               myHandles.fixedline.Position = evt.CurrentPosition;
           end
            
            figure(2);
            hold on
            % Get position
            diagonal = [];
            % Find chosen data
            myHandles.currentidx = find(myHandles.idx);
            
            for i = 1:length(myHandles.currentidx)
                    gammaon = 0;
                    switch myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).datatype
                        case '2D'
                            if fixval == 1 && (myHandles.currentidx(i) == myHandles.refidx)
                                linestr = 'fixedline';
                            else
                                linestr = 'myline';
                            end
                            
                            % Create a line mask. Matrix which is bound by the line ends in current axes (minimum dim size == 1)
                            % Find incides from current file
                            if strcmp(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).Plane, 'XY')
                                [~,x1] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).xpos - myHandles.(linestr).Position(1,1)));
                                [~,x2] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).xpos - myHandles.(linestr).Position(2,1)));
                                [~,y1] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).ypos - myHandles.(linestr).Position(1,2)));
                                [~,y2] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).ypos - myHandles.(linestr).Position(2,2)));
                            elseif strcmp(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).Plane, 'YZ')
                                [~,x1] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).zpos - myHandles.(linestr).Position(1,1)));
                                [~,x2] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).zpos - myHandles.(linestr).Position(2,1)));
                                [~,y1] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).ypos - myHandles.(linestr).Position(1,2)));
                                [~,y2] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).ypos - myHandles.(linestr).Position(2,2)));
                            else
                                [~,x1] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).zpos - myHandles.(linestr).Position(1,1)));
                                [~,x2] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).zpos - myHandles.(linestr).Position(2,1)));
                                [~,y1] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).xpos - myHandles.(linestr).Position(1,2)));
                                [~,y2] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).xpos - myHandles.(linestr).Position(2,2)));
                            end
                         if isfield(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]), 'chosengam') && myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).chosengam ~= 0
                            gammaon = 1;
                            gamname = myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).chosengamname;
                            [~,x11] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).GAMposx - myHandles.myline.Position(1,1)));
                            [~,x22] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).GAMposx - myHandles.myline.Position(2,1)));
                            [~,y11] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).GAMposy - myHandles.myline.Position(1,2)));
                            [~,y22] = min(abs(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).GAMposy - myHandles.myline.Position(2,2)));
                         end
                            % Line ends are always in the corners of the datamatrix -> get
                            % diagonal values.

                            % Define the initial data bound by the created box
                            if isfield(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]), 'chosengam') && myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).chosengam ~= 0
                                switch myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).currentgaman
                                    case 'DTA'
                                        datamask2 = myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).dta(min([y11,y22]):max([y11,y22]), min([x11,x22]):max([x11,x22]));
                                        myHandles.minval = min(min(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).dta));
                                    case 'DD'
                                        datamask2 = myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).absoluteDD(min([y11,y22]):max([y11,y22]), min([x11,x22]):max([x11,x22]));
                                        myHandles.minval = min(min(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).absoluteDD));
                                    case 'gamma'
                                        datamask2 = myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).gamma(min([y11,y22]):max([y11,y22]), min([x11,x22]):max([x11,x22]));
                                        myHandles.minval = min(min(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(gamname).gamma));
                                end
                                
                                % Flip data so that the line defines a diagonal from [0,0] to
                                % [end, end]
                                 if x11 > x22
                                     datamask2 = flip(datamask2,2);
                                 end
                                 if y11 > y22
                                    datamask2 = flip(datamask2,1);
                                 end
                                 % Get size
                                 [ysiz2, xsiz2] = size(datamask2);

                                 % Check if diagonal can be computed (dimsize > 1),
                                 % otherwise get the horizontal/vertical vector
                                 if ysiz2 > 1 && xsiz2 > 1

                                     [X, Y] = meshgrid(1:xsiz2, 1:ysiz2);
                                     if xsiz2 <= ysiz2
                                         xax = linspace(1,xsiz2,ysiz2);
                                         [Xq, Yq] = meshgrid(xax,1:ysiz2);
                                         datainterpmask = interp2(X,Y,datamask2,Xq,Yq);
                                         diagonal2 = diag(datainterpmask);
                                     elseif ysiz2 <= xsiz2
                                         yax = linspace(1,ysiz2,xsiz2);
                                         [Xq, Yq] = meshgrid(1:xsiz2, yax);
                                         datainterpmask = interp2(X,Y,datamask2,Xq,Yq);
                                         diagonal2 = diag(datainterpmask);
                                     end
                                 elseif ysiz2 < 2
                                        diagonal2 = datamask2(1, 1:xsiz2);
                                 elseif xsiz2 < 2
                                        diagonal2 = datamask2(1:ysiz2, 1);
                                 end
                                 
                                 xposvect2 = linspace(0,sqrt((myHandles.myline.Position(1,2) - myHandles.myline.Position(2,2))^2 + (myHandles.myline.Position(2,1) - myHandles.myline.Position(1,1))^2),length(diagonal2));
                            end
                            % Draw the original data
                            datamask = myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(['data', myHandles.normalizationtype])(min([y1,y2]):max([y1,y2]), min([x1,x2]):max([x1,x2]));
                                myHandles.minval = min(min(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(['data', myHandles.normalizationtype])));
                            % Flip data so that the line defines a diagonal from [0,0] to
                            % [end, end]
                             if x1 > x2
                                 datamask = flip(datamask,2);
                             end
                             if y1 > y2
                                datamask = flip(datamask,1);
                             end
                             % Get size
                             [ysiz, xsiz] = size(datamask);

                             % Check if diagonal can be computed (dimsize > 1),
                             % otherwise get the horizontal/vertical vector
                             if ysiz > 1 && xsiz > 1

                             [X, Y] = meshgrid(1:xsiz, 1:ysiz);
                             if xsiz <= ysiz
                                 xax = linspace(1,xsiz,ysiz);
                                 [Xq, Yq] = meshgrid(xax,1:ysiz);
                                 datainterpmask = interp2(X,Y,datamask,Xq,Yq);
                                 diagonal = diag(datainterpmask);
                             elseif ysiz <= xsiz
                                 yax = linspace(1,ysiz,xsiz);
                                 [Xq, Yq] = meshgrid(1:xsiz, yax);
                                 datainterpmask = interp2(X,Y,datamask,Xq,Yq);
                                 diagonal = diag(datainterpmask);
                             end

                             %diagmat = diag(diag(datainterpmask));

    %                          prof2d = datainterpmask(1:size(diagmat,1), 1:size(diagmat,2)) - diagmat;
    %                          figure(3);
    %                         surf(prof2d);
    %                         view(0,90);
                             elseif ysiz < 2
                                 diagonal = datamask(1, 1:xsiz);
                             elseif xsiz < 2
                                 diagonal = datamask(1:ysiz, 1);
                             end


                             xposvect = linspace(0,sqrt((myHandles.myline.Position(1,2) - myHandles.myline.Position(2,2))^2 + (myHandles.myline.Position(2,1) - myHandles.myline.Position(1,1))^2),length(diagonal));

                        case '3D'
                            %plot(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).(['data', myHandles.normalizationtype]));
                            disp('Error. Cannot create a profile for 3D data. Extract 2D plane first');
                        case '1D'
                            disp('Error. Cannot create a profile for 1D data...');
                    end

                    
                    
                    set(myHandles.lineplot(i), 'XData', xposvect, 'YData', diagonal); 
                    if gammaon == 1
                        set(myHandles.lineplot2(i), 'XData', xposvect2, 'YData', diagonal2);
                        yyaxis right
                    end
                    
                    axis([0 inf myHandles.minval inf]);
                    if isfield(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]), 'chosengam') && myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).chosengam ~= 0
                        ylabel(myHandles.alldata.(['Data', num2str(myHandles.currentidx(i))]).currentgaman);
                     elseif strcmp(myHandles.normalizationtype, 'NON') 
                       ylabel('Dose');
                    else
                       ylabel('Value [normalized]'); 
                    end
                    xlabel('Pos [cm]');
            end
            
            % Plot the line
            %posvectstarts = [min([y1,y2]), min([x1,x2])];
            %posvectcomponents = [posvectstarts(1):1/leny:leny];

            %guidata(hObject, myHandles);
            
    end






% --- Executes on button press in checkboxfixref.
function checkboxfixref_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxfixref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.reffixed == 1
   try
      delete(handles.fixedline);
   catch
      disp('No reference line to be deleted');
   end
end
handles.reffixed = get(hObject,'Value');


guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of checkboxfixref


% --- Executes on button press in selectallbutton.
function selectallbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectallbutton (see GCBO)
mask = ones(1,length(handles.idx));
if sum(handles.idx) ~= length(handles.idx)
    handles.idx = handles.idx|mask;
else
    handles.idx = ~(handles.idx|mask);
end


% Update the figure with only re-drawing chosen data
handles.redraw = 1;
guidata(hObject, handles);


try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fileMenuSave_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fileMenuExit_Callback(hObject, eventdata, handles)
fclose all;
close all;
clc;
closereq(); 


% hObject    handle to fileMenuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function editMenu_Callback(hObject, eventdata, handles)
% hObject    handle to editMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function viewMenu_Callback(hObject, eventdata, handles)
% hObject    handle to viewMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function toolsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to toolsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function layoutMenu_Callback(hObject, eventdata, handles)
% hObject    handle to layoutMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function toolsMenuAssignCurve_Callback(hObject, eventdata, handles)
% hObject    handle to toolsMenuAssignCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function assignProfile_Callback(hObject, eventdata, handles)
% hObject    handle to assignProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Find chosen data
indices = find(handles.idx);

for i = 1:length(indices)
    if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D') || ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
        % Check that curvetype is not already a profile
        if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'profile')
            try
                [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL] = fieldparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                    handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1), 'Data', handles.caxcorrection);
                handles.alldata.(['Data', num2str(indices(i))]).params = [caxdev, sym, hom, dmax, dmin, dev, FW, penR, penL];
                handles.alldata.(['Data', num2str(indices(i))]).datatype = 'profile';
                disp([handles.alldata.(['Data', num2str(indices(i))]).name, ' assigned as profile succesfully']);
                % Compute profile parameters
            catch
                disp('Error in assignProfile function. Could not compute profile parameters.');
            end
            
        else
            disp('Data already assigned as profile');
        end
    else
        disp('Only 1D data can be assigned as Profile or PDD');
    end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;
guidata(hObject, handles);


try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end


% --------------------------------------------------------------------
function assignPDD_Callback(hObject, eventdata, handles)
% hObject    handle to assignPDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Find chosen data
indices = find(handles.idx);

for i = 1:length(indices)
    if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '2D') || ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, '3D')
        % Check that curvetype is not already a profile
        if ~strcmp(handles.alldata.(['Data', num2str(indices(i))]).datatype, 'pdd')
            try
                [R100, R80, R50, D100, D200, Ds, J1020] = pddparams(handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,2),...
                    handles.alldata.(['Data', num2str(indices(i))]).interpolated(:,1));
                handles.alldata.(['Data', num2str(indices(i))]).params = [R100, R80, R50, D100, D200, Ds, J1020]; 
                handles.alldata.(['Data', num2str(indices(i))]).datatype = 'pdd';
                disp([handles.alldata.(['Data', num2str(indices(i))]).name, ' assigned as PDD succesfully']);
                % Compute profile parameters
            catch
                disp('Error in assignProfile function. Could not compute PDD parameters.');
            end
            
        else
            disp('Data already assigned as profile');
        end
    else
        disp('Only 1D data can be assigned as Profile or PDD');
    end
end

% Update the figure with only re-drawing chosen data
handles.redraw = 1;
guidata(hObject, handles);


try
    listboxChooseData_Callback(handles.listboxChooseData, [], handles);
catch
   disp('Error in axis invertion. Figure not updated'); 
end



