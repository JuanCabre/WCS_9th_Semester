function varagout = afc(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AFC :: APNet FDTD Code %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% APNet 3rd generation FDTD code
%% (c) 2006-- Ondrej Franek, Mauro Pelosi

version = 'v3.4'; % 24-11-2011

%% default setting of variables
DefaultVars = struct('eps0'            ,8.854e-12, ...
                     'mu0'             ,1.257e-6, ...
                     'c'               ,1/sqrt(8.854e-12*1.257e-6), ...
                     'SpaceDepSwitch'  ,[0 0 1], ...
                     'Dx'              ,5e-3, ...
                     'Dt_mode'         ,'auto', ...
                     'Dt'              ,1e-12, ...
                     'CFL_coef'        ,0.999, ...
                     'Df'              ,1e6, ...
                     'Nt'              ,1000, ...
                     'Drop_dB'         ,-30, ...
                     'Drop2_dB'        ,-30, ...
                     'drop_flag'       ,false, ...
                     'drop2_flag'      ,false, ...
                     'limits'          ,[0 0 0;10 10 10], ...
                     'pml'             ,true, ...
                     'field_flag'      ,true, ...
                     'nfff_flag'       ,true, ...
                     'nfff_faces'      ,[1 1 1 1 1 1], ...
                     'energy_flag'     ,true, ...
                     'e_step'          ,100, ...
                     'pmax'            ,8, ...
                     'm_s'             ,4, ...
                     'm_k'             ,4, ...
                     'm_a'             ,1, ...
                     'sig_max'         ,(4 + 1)/150/pi/5e-3, ...
                     'kappa_max'       ,15, ...
                     'a_max'           ,0.04, ...
                     'freq_nf'         ,900e6, ...
                     'freq_ff'         ,900e6, ...
                     'x_nfff'          ,[1 1], ...
                     'y_nfff'          ,[1 1], ...
                     'z_nfff'          ,[1 1], ...
                     'phcenter'        ,[0 0 0], ...
                     'resolution'      ,[90 180], ...
                     'fc'              ,900e6, ...
                     'fs'              ,300e6, ...
                     'offset'          ,370, ...
                     'mag'             ,1, ...
                     'fmin'            ,600e6, ...
                     'fmax'            ,1200e6, ...
                     'fd_thres'        ,-20, ...
                     'td_thres'        ,-60, ...
                     'GaussSineShape'  ,'auto', ...
                     'ExcitationPulse' ,'gauss_sine', ...
                     'Geometry'        ,{{}}, ...
                     'GeometryString'  ,[], ...
                     'Engine'          ,struct('Name','Matlab'), ...
                     'JobDir'          ,[], ...
                     'MuAveraging'     ,1, ...
                     'FieldNormalizing',1, ...
                     'locked_flag'     ,false, ...
                     'Vs'              ,[], ...
                     'V'               ,[], ...
                     'I'               ,[], ...
                     'Vs_fft'          ,[], ...
                     'V_fft'           ,[], ...
                     'I_fft'           ,[], ...
                     'Z_fft'           ,[], ...
                     'Rs'              ,[], ...
                     'P_in'            ,[], ...
                     'P_rad_nf'        ,[], ...
                     'P_rad_ff'        ,[], ...
                     'P_rad_ff_theta'  ,[], ...
                     'P_rad_ff_phi'    ,[], ...
                     'Probes'          ,[], ...
                     'Eth'             ,[], ...
                     'Eph'             ,[], ...
                     'Ex_field'        ,[], ...
                     'Ey_field'        ,[], ...
                     'Ez_field'        ,[], ...
                     'energy'          ,[], ...
                     'SAR'             ,[], ...
                     'SAR_Aver_1g'     ,[], ...
                     'SAR_Aver_10g'    ,[], ...
                     'E_field'         ,[], ...
                     'H_field'         ,[], ...
                     'J_field'         ,[], ...
                     'overall'         ,[], ...
                     'time_start'      ,[], ...
                     'time_end'        ,[], ...
                     'partial_times'   ,[] );

%% default time step computed automatically
DefaultVars.Dt = 5e-3/((1/sqrt(8.854e-12*1.257e-6))*sqrt(3))*0.999;
DefaultVars.offset = 370*DefaultVars.Dt;

%% variable initialization - these have to be global
v = [];
xmax = []; ymax = []; zmax = [];
eps_rx = []; eps_ry = []; eps_rz = [];
mu_rx = [];  mu_ry = [];  mu_rz = [];
sig_x = [];  sig_y = [];  sig_z = [];
kappa_max = []; a_max = []; m_s = []; m_k = []; m_a = [];
SourceSubs = []; SourceDir = [];
x_nfff = []; y_nfff = []; z_nfff = []; f_patt = [];
phase_center = []; resolution = [];
nfff_flag = []; field_flag = []; energy_flag = []; e_step = [];
Rs = []; hardsource_flag = [];
Cax = []; Cay = []; Caz = [];
Cbx = []; Cby = []; Cbz = [];
Dbx = []; Dby = []; Dbz = [];
s_x = []; s_y = []; s_z = [];
in = []; out = []; overall = []; limits = [];
ActualFileName = []; ActualPathName = [];
GeomAxes = []; GeometryFig = []; Engine = [];
SuccessFlag = false;
plot_axes = []; probe_pointer = []; freq_pointer = []; pol_label = [];

%% string lists
Units      = {'cells'};
Directions = {'+x','-x','+y','-y','+z','-z'};
Components = {'Ex','Ey','Ez','Hx','Hy','Hz'};
Renderers  = {'painters','zbuffer','opengl'};

%% some more defaults
DefaultSourceResistance = 50;
DefaultMassDensity      = 1e3;

%% tasklist filename
TasklistName  = 'tasklist.txt';

%% Figure title
FigTitle = ['AFC ',version];

%% Geometry list strings
InvisibilityChar = '# ';
SeparationString = '  ::  ';

%% initial setting of directories
AFCDir    = pwd;
WorkDir   = pwd;
ExtObjDir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default font properties and colors
if isunix, %% default setting for Fyrkat
    DefaultFont.FontName = 'Verdana';
    DefaultFont.FontSize = 10;
    DefaultColor = [220 220 220]./255;
elseif ispc, %% default setting for Windows XP
    DefaultFont.FontName = 'Tahoma';
    DefaultFont.FontSize = 8;
    DefaultColor = get(0,'DefaultUIControlBackgroundColor');
else
end;
DefaultFont.FontUnits  = 'points';
DefaultFont.FontAngle  = 'normal';
DefaultFont.FontWeight = 'normal';

%% default engines
DefaultEngines = {{}};
DefaultEngines{1} = struct('Name','Matlab'); %% the only default is Matlab

DefaultPrefs = struct('Font',DefaultFont, ...
                      'Color',DefaultColor, ...
                      'Renderer','zbuffer', ...
                      'FortranDir',fullfile(AFCDir,'fortran'), ...
                      'DisplayGrid',1, ...
                      'CheckBrickLimits',1, ...
                      'NewSaving',1, ...
                      'NewInOutFiles',1, ...
                      'CacheOptim',0, ...
                      'Engines',DefaultEngines, ...
                      'ActiveEngine',1);

%% preferences file
if exist(fullfile(AFCDir,'preferences.mat'),'file')
    
    pref = load(fullfile(AFCDir,'preferences.mat'));
    
    %% fill non-existent fields with default values
    PrefFields = fieldnames(DefaultPrefs);
    for i = 1:length(PrefFields),
        if ~isfield(pref,PrefFields{i}),
            pref.(PrefFields{i}) = DefaultPrefs.(PrefFields{i});
        end;
    end;
    
else
    
    pref = DefaultPrefs;
    save(fullfile(AFCDir,'preferences.mat'),'-struct','pref');
    
end;

if pref.DisplayGrid, PatchLineStyle = '-';
else                 PatchLineStyle = 'none';
end;

if pref.NewSaving, SavingVersion = '-v7.3';
else               SavingVersion = '-v7';
end;

%% setting the GUI font properties
for name = fieldnames(pref.Font).',
    set(0,['DefaultUIPanel'  ,name{1}],pref.Font.(name{1}));
    set(0,['DefaultUIControl',name{1}],pref.Font.(name{1}));
end;

%% setting the GUI background color
set(0,'DefaultUIPanelBackgroundColor'  ,pref.Color);
set(0,'DefaultUIControlBackgroundColor',pref.Color);

%% shorter name for the active engine
Engine = pref.Engines(pref.ActiveEngine);

%% lock image :-)
RGB = [ pref.Color ; 0 0 0 ; .8 .8 1 ];
R_Data = RGB(:,1); G_Data = RGB(:,2); B_Data = RGB(:,3);
string = [ '000011100111222002212221121122221101212211011112' ; ...
           '000111000112221012222220211222212002122120021111' ];
data = zeros(16,8); data(33:end) = hex2dec(string(:));
data = [ data data(:,end:-1:1) ] + 1;
LockImageData = cat(3,R_Data(data),G_Data(data),B_Data(data));
LockImageAxes = [];

%% obtaining the default screen resolution
ScreenSize = get(0,'ScreenSize');
ScreenCenter = ScreenSize(3:4)./2;

%% figure dimensions in pixels
FigWidth = 265;
FigHeight = 605;

%% creating the main figure
MainFig = figure('Visible','off', ...
                 'Color',pref.Color, ...
                 'MenuBar','none', ...
                 'Name',FigTitle, ...
                 'NumberTitle','off', ...
                 'Units','pixels', ...
                 'Position',[10 ScreenSize(4)-FigHeight-60 FigWidth FigHeight], ...
                 'Renderer','painters', ...
                 'Resize','off', ...
                 'Toolbar','none', ...
                 'DockControls','off', ...
                 'CloseRequestFcn',@Exit);

             
             
%%%%%%%%%%%%%%
%% MENU BAR %%
%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
MenuFile   = uimenu(MainFig,'Label','File');

    MenuNew  = uimenu(MenuFile,'Label','New...', ...
                               'Accelerator','N', ...
                               'Callback',@Open); %% this calls Open too
    MenuOpen = uimenu(MenuFile,'Label','Open...', ...
                               'Accelerator','O', ...
                               'Callback',@Open);
    MenuClose = uimenu(MenuFile,'Label','Close', ...
                                'Accelerator','W', ...
                                'Callback',@Close, ...
                                'Enable','off');
    % --------------------------------------------------------------------
    MenuSave = uimenu(MenuFile,'Label','Save', ...
                               'Accelerator','S', ...
                               'Separator','on', ...
                               'Callback',@Save, ...
                               'Enable','off');
    MenuSaveAs = uimenu(MenuFile,'Label','Save As...', ...
                                 'Callback',@SaveAs, ...
                                 'Enable','off');
    % --------------------------------------------------------------------
    MenuExport = uimenu(MenuFile,'Label','Export...', ...
                                 'Separator','on', ...
                                 'Callback',@Export, ...
                                 'Enable','off');
    % --------------------------------------------------------------------
    MenuDeployMulti = uimenu(MenuFile,'Label','Deploy Multiple...', ...
                                      'Separator','on', ...
                                      'Callback',@SimulationDeploy);
    MenuCollectMulti = uimenu(MenuFile,'Label','Collect Multiple...', ...
                                       'Callback',@CollectMulti);
    % --------------------------------------------------------------------
    MenuCompareResults = uimenu(MenuFile,'Label','Compare Results', ...
                                         'Separator','on');
    % --------------------------------------------------------------------
    MenuBDFConvert = uimenu(MenuFile,'Label','BDF File Converter...', ...
                                     'Separator','on', ...
                                     'Callback',['run(''', ...
                                     fullfile(AFCDir,'converter.m'), ...
                                                       ''');']);
    % --------------------------------------------------------------------
    MenuPreferences = uimenu(MenuFile,'Label','Preferences...', ...
                                      'Separator','on', ...
                                      'Callback',@Preferences);
    % --------------------------------------------------------------------
    MenuExit = uimenu(MenuFile,'Label','Exit', ...
                               'Accelerator','Q', ...
                               'Separator','on', ...
                               'Callback',@Exit);
                           
% ------------------------------------------------------------------------
MenuInsert = uimenu(MainFig,'Label','Insert');

    MenuWire = uimenu(MenuInsert,'Label','Wire...', ...
                                 'Callback',{@Wire,0}, ...
                                 'Enable','off');
    % --------------------------------------------------------------------
    MenuTriangle = uimenu(MenuInsert,'Label','Triangle...', ...
                                     'Callback',{@Triangle,0}, ...
                                     'Separator','on', ...
                                     'Enable','off');
    MenuRectangle = uimenu(MenuInsert,'Label','Rectangle...', ...
                                      'Callback',{@Rectangle,0}, ...
                                      'Enable','off');
    % --------------------------------------------------------------------
    MenuBrick = uimenu(MenuInsert,'Label','Brick...', ...
                                  'Callback',{@Brick,0}, ...
                                  'Separator','on', ...
                                  'Enable','off');
    % --------------------------------------------------------------------
    MenuSource = uimenu(MenuInsert,'Label','Source...', ...
                                   'Callback',{@Source,0}, ...
                                   'Separator','on', ...
                                   'Enable','off');
    MenuFieldSource = uimenu(MenuInsert,'Label','Field Source...', ...
                                        'Callback',{@FieldSource,0}, ...
                                        'Enable','off');
    % --------------------------------------------------------------------
    MenuResistor = uimenu(MenuInsert,'Label','Resistor...', ...
                                     'Callback',{@Resistor,0}, ...
                                     'Separator','on', ...
                                     'Enable','off');
    MenuCapacitor = uimenu(MenuInsert,'Label','Capacitor...', ...
                                      'Callback',{@Capacitor,0}, ...
                                      'Enable','off');
    % --------------------------------------------------------------------
    MenuProbe = uimenu(MenuInsert,'Label','Probe...', ...
                                  'Callback',{@Probe,0}, ...
                                  'Separator','on', ...
                                  'Enable','off');
    % --------------------------------------------------------------------
    MenuExtObj = uimenu(MenuInsert,'Label','External Object...', ...
                                   'Callback',{@ExtObj,0}, ...
                                   'Separator','on', ...
                                   'Enable','off');
                               
% ------------------------------------------------------------------------
MenuSimulation = uimenu(MainFig,'Label','Simulation');

    MenuDeploy  = uimenu(MenuSimulation,'Label','Deploy...', ...
                                        'Callback',@SimulationDeploy, ...
                                        'Enable','off');
    % --------------------------------------------------------------------
    MenuSubmit  = uimenu(MenuSimulation,'Label','Submit to queue', ...
                                        'Separator','on', ...
                                        'Callback',@SimulationStart, ...
                                        'Enable','off');
    MenuRun     = uimenu(MenuSimulation,'Label','Run', ...
                                        'Callback',@SimulationStart, ...
                                        'Enable','off');
    % --------------------------------------------------------------------
    MenuCollect = uimenu(MenuSimulation,'Label','Collect', ...
                                        'Separator','on', ...
                                        'Callback',@SimulationCollect, ...
                                        'Enable','off');

% ------------------------------------------------------------------------
MenuResults = uimenu(MainFig,'Label','Results');

    MenuSourceVoltage = uimenu(MenuResults, ...
        'Label','Source Voltage', ...
        'Enable','off');
        uimenu(MenuSourceVoltage, ...
            'Label','Time Domain', ...
            'Callback',{@ResponseDisplay,'Source Voltage','Time'});
        uimenu(MenuSourceVoltage, ...
            'Label','Frequency Domain Magnitude',...
            'Callback',{@ResponseDisplay,'Source Voltage','Magnitude'});
        uimenu(MenuSourceVoltage, ...
            'Label','Frequency Domain Phase',...
            'Callback',{@ResponseDisplay,'Source Voltage','Phase'});

    MenuInputVoltage = uimenu(MenuResults, ...
        'Label','Input Voltage', ...
        'Enable','off');
        uimenu(MenuInputVoltage, ...
            'Label','Time Domain', ...
            'Callback',{@ResponseDisplay,'Input Voltage','Time'});
        uimenu(MenuInputVoltage, ...
            'Label','Frequency Domain Magnitude',...
            'Callback',{@ResponseDisplay,'Input Voltage','Magnitude'});
        uimenu(MenuInputVoltage, ...
            'Label','Frequency Domain Phase',...
            'Callback',{@ResponseDisplay,'Input Voltage','Phase'});

    MenuInputCurrent = uimenu(MenuResults, ...
        'Label','Input Current', ...
        'Enable','off');
        uimenu(MenuInputCurrent, ...
            'Label','Time Domain', ...
            'Callback',{@ResponseDisplay,'Input Current','Time'});
        uimenu(MenuInputCurrent, ...
            'Label','Frequency Domain Magnitude',...
            'Callback',{@ResponseDisplay,'Input Current','Magnitude'});
        uimenu(MenuInputCurrent, ...
            'Label','Frequency Domain Phase',...
            'Callback',{@ResponseDisplay,'Input Current','Phase'});
        
    MenuProbeResults = uimenu(MenuResults, ...
        'Label','Probes', ...
        'Enable','off');
        uimenu(MenuProbeResults, ...
            'Label','Time Domain', ...
            'Callback',{@ProbeResults,'Time'});
        uimenu(MenuProbeResults, ...
            'Label','Frequency Domain Magnitude',...
            'Callback',{@ProbeResults,'Magnitude'});
        uimenu(MenuProbeResults, ...
            'Label','Frequency Domain Phase',...
            'Callback',{@ProbeResults,'Phase'});

    % --------------------------------------------------------------------
    MenuImpedance = uimenu(MenuResults,'Label','Impedance', ...
                                 'Separator','on', ...
                                 'Enable','off');
        uimenu(MenuImpedance,'Label','Real Part', ...
                             'Callback',{@Impedance,'real'});
        uimenu(MenuImpedance,'Label','Imaginary Part', ...
                             'Callback',{@Impedance,'imag'});
        uimenu(MenuImpedance,'Label','Magnitude', ...
                             'Callback',{@Impedance,'mag'});
        uimenu(MenuImpedance,'Label','Phase', ...
                             'Callback',{@Impedance,'phase'});
        % ----------------------------------------------------------------
        uimenu(MenuImpedance,'Label','Quality Factor (Q)', ...
                             'Separator','on', ...
                             'Callback',{@Impedance,'Q'});
        % ----------------------------------------------------------------
        uimenu(MenuImpedance,'Label','S11 Parameter', ...
                             'Separator','on', ...
                             'Callback',{@Impedance,'s11'});
        uimenu(MenuImpedance,'Label','VSWR', ...
                             'Callback',{@Impedance,'vswr'});
        % ----------------------------------------------------------------
        uimenu(MenuImpedance,'Label','Smith Chart', ...
                             'Separator','on', ...
                             'Callback',@SmithChart);

    MenuPatt = uimenu(MenuResults,'Label','Radiation Pattern', ...
                                  'Enable','off');
                              
        MenuPattTotal = uimenu(MenuPatt,'Label','Total');
        
            uimenu(MenuPattTotal,'Label','3D Flat', ...
                                 'Callback',{@Pattern,'flat'});
            uimenu(MenuPattTotal,'Label','3D Flat [dB]', ...
                                 'Callback',{@Pattern,'flat_db'});
            uimenu(MenuPattTotal,'Label','3D Sphere', ...
                                 'Callback',{@Pattern,'sphere'});
            uimenu(MenuPattTotal,'Label','3D Sphere [dB]', ...
                                 'Callback',{@Pattern,'sphere_db'});
            % ------------------------------------------------------------
            uimenu(MenuPattTotal,'Label','Polar X-Y Plane', ...
                                 'Separator','on', ...
                                 'Callback',{@Pattern,'xy','polar'});
            uimenu(MenuPattTotal,'Label','Polar X-Z Plane', ...
                                 'Callback',{@Pattern,'xz','polar'});
            uimenu(MenuPattTotal,'Label','Polar Y-Z Plane', ...
                                 'Callback',{@Pattern,'yz','polar'});
            % ------------------------------------------------------------
            uimenu(MenuPattTotal,'Label','Polar X-Y Plane [dB]', ...
                                 'Separator','on', ...
                                 'Callback',{@Pattern,'xy_db','polar'});
            uimenu(MenuPattTotal,'Label','Polar X-Z Plane [dB]', ...
                                 'Callback',{@Pattern,'xz_db','polar'});
            uimenu(MenuPattTotal,'Label','Polar Y-Z Plane [dB]', ...
                                 'Callback',{@Pattern,'yz_db','polar'});
            % ------------------------------------------------------------
            uimenu(MenuPattTotal,'Label','Cartesian X-Y Plane', ...
                                 'Separator','on', ...
                                 'Callback',{@Pattern,'xy','cart'});
            uimenu(MenuPattTotal,'Label','Cartesian X-Z Plane', ...
                                 'Callback',{@Pattern,'xz','cart'});
            uimenu(MenuPattTotal,'Label','Cartesian Y-Z Plane', ...
                                 'Callback',{@Pattern,'yz','cart'});
            % ------------------------------------------------------------
            uimenu(MenuPattTotal,'Label','Cartesian X-Y Plane [dB]', ...
                                 'Separator','on', ...
                                 'Callback',{@Pattern,'xy_db','cart'});
            uimenu(MenuPattTotal,'Label','Cartesian X-Z Plane [dB]', ...
                                 'Callback',{@Pattern,'xz_db','cart'});
            uimenu(MenuPattTotal,'Label','Cartesian Y-Z Plane [dB]', ...
                                 'Callback',{@Pattern,'yz_db','cart'});
        
        MenuPattTheta  = copyobj(MenuPattTotal,MenuPatt);
        MenuPattPhi    = copyobj(MenuPattTotal,MenuPatt);
        
        set(MenuPattTheta,'Label','Theta polarized');
        set(MenuPattPhi  ,'Label','Phi polarized');
                                   
    MenuNearFields = uimenu(MenuResults,'Label','Near Fields', ...
                                        'Enable','off');
        MenuExField = uimenu(MenuNearFields, ...
            'Label','Ex-field', ...
            'Callback',@ArrayDisplay);
        MenuEyField = uimenu(MenuNearFields, ...
            'Label','Ey-field', ...
            'Callback',@ArrayDisplay);
        MenuEzField = uimenu(MenuNearFields, ...
            'Label','Ez-field', ...
            'Callback',@ArrayDisplay);
        % ----------------------------------------------------------------
        MenuHxField = uimenu(MenuNearFields, ...
            'Label','Hx-field', ...
            'Separator','on', ...
            'Callback',@ArrayDisplay);
        MenuHyField = uimenu(MenuNearFields, ...
            'Label','Hy-field', ...
            'Callback',@ArrayDisplay);
        MenuHzField = uimenu(MenuNearFields, ...
            'Label','Hz-field', ...
            'Callback',@ArrayDisplay);
        % ----------------------------------------------------------------
        MenuEmagField = uimenu(MenuNearFields, ...
            'Label','E-field Magnitude', ...
            'Enable','off', ...
            'Separator','on', ...
            'Callback',@ArrayDisplay);
        MenuHmagField = uimenu(MenuNearFields, ...
            'Label','H-field Magnitude', ...
            'Enable','off', ...
            'Callback',@ArrayDisplay);
        MenuJmagField = uimenu(MenuNearFields, ...
            'Label','J-field Magnitude', ...
            'Enable','off', ...
            'Callback',@ArrayDisplay);
        % ----------------------------------------------------------------
        MenuSARDisp = uimenu(MenuNearFields, ...
            'Label','SAR Local', ...
            'Enable','off', ...
            'Separator','on', ...
            'Callback',@ArrayDisplay);
        MenuSAR1gDisp = uimenu(MenuNearFields, ...
            'Label','SAR 1g-averaged', ...
            'Enable','off', ...
            'Callback',@ArrayDisplay);
        MenuSAR10gDisp = uimenu(MenuNearFields, ...
            'Label','SAR 10g-averaged', ...
            'Enable','off', ...
            'Callback',@ArrayDisplay);
        
    MenuCurrents = uimenu(MenuResults,'Label','Wire Currents...', ...
                                      'Enable','off', ...
                                      'Callback',@WireCurrents);
    MenuEnergy = uimenu(MenuResults,'Label','Energy', ...
                                    'Enable','off', ...
                                    'Callback',@EnergyDisplay);
    % --------------------------------------------------------------------
    MenuEfficiency = uimenu(MenuResults,'Label','Radiation Efficiency', ...
                                        'Separator','on', ...
                                        'Enable','off', ...
                                        'Callback',@EfficiencyDisplay);
    MenuQandBW = uimenu(MenuResults,'Label','Q and Bandwidth', ...
                                    'Enable','off', ...
                                    'Callback',@QandBWDisplay);
    MenuCPUTime = uimenu(MenuResults,'Label','CPU Time', ...
                                     'Enable','off', ...
                                     'Callback',@CPUTime);
    % --------------------------------------------------------------------
    MenuPostprocessing = uimenu(MenuResults, ...
        'Label','Postprocessing', ...
        'Separator','on', ...
        'Enable','off');
        MenuSARDist = uimenu(MenuPostprocessing, ...
            'Label','SAR Distribution', ...
            'Callback',@SARDist);
        MenuSARAver1g  = uimenu(MenuPostprocessing, ...
            'Label','SAR Averaging 1g', ...
            'Enable','off', ...
            'Callback',@SAR_Average);
        MenuSARAver10g = uimenu(MenuPostprocessing, ...
            'Label','SAR Averaging 10g', ...
            'Enable','off', ...
            'Callback',@SAR_Average);
        % ----------------------------------------------------------------
        MenuEfield  = uimenu(MenuPostprocessing, ...
            'Label','E-field Magnitude', ...
            'Separator','on', ...
            'Callback',@Efield);
        MenuHfield  = uimenu(MenuPostprocessing, ...
            'Label','H-field Magnitude', ...
            'Callback',@Hfield);
        MenuJfield  = uimenu(MenuPostprocessing, ...
            'Label','Current Density Magnitude', ...
            'Callback',@Jfield);
        
%% copy the Results submenus to CompareResults
ch2 = copyobj(get(MenuResults,'Children'),MenuCompareResults);

%% delete those which cannot be compared
delete(findobj(MenuCompareResults,'Label','Postprocessing','-or', ...
                                  'Label','Wire Currents...','-or', ...
                                  'Label','Near Fields','-or', ...
                                  'Label','3D Flat','-or', ...
                                  'Label','3D Flat [dB]','-or', ...
                                  'Label','3D Sphere','-or', ...
                                  'Label','3D Sphere [dB]','-or', ...
                                  'Label','Polar X-Y Plane','-or', ...
                                  'Label','Polar X-Z Plane','-or', ...
                                  'Label','Polar Y-Z Plane','-or', ...
                                  'Label','Polar X-Y Plane [dB]','-or', ...
                                  'Label','Polar X-Z Plane [dB]','-or', ...
                                  'Label','Polar Y-Z Plane [dB]','-or', ...
                                  'Label','Radiation Efficiency','-or', ...
                                  'Label','Q and Bandwidth','-or', ...
                                  'Label','CPU Time'));

%% find all submenus of Compare Results
CompareSubmenus = findobj(MenuCompareResults);

set(CompareSubmenus,'Enable','on'); %% enable them

%% move callbacks to userdata and put own callback instead
for h = CompareSubmenus',
    cb = get(h,'Callback');
    if ~isempty(cb),
        set(h,'UserData',cb,'Callback',@CompareResults);
    end;
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PANEL DEFINITIONS AND CALLBACKS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%
%% MAIN PANEL %%
%%%%%%%%%%%%%%%%
Panel = uipanel('Parent',MainFig, ...
                'Units','pixels', ...
                'Position',[0 -1 FigWidth+3 FigHeight+4]);

%% panels position parameters            
PanelsWidth = FigWidth - 11;
PanelsLeft = 5;
VerticalPosition = FigHeight + 2;
Separate = 2;



%%%%%%%%%%%%%%%%%%%%
%% GEOMETRY PANEL %%
%%%%%%%%%%%%%%%%%%%%
PanelHeight = 160;
VerticalPosition = VerticalPosition - PanelHeight - Separate;
GeomPanel = uipanel('Parent',Panel, ...
                    'ForegroundColor','b', ...
                    'Title','Geometry', ...
                    'Units','pixels', ...
                    'Position',[PanelsLeft VerticalPosition ...
                                PanelsWidth PanelHeight]);
GeomListBox = uicontrol(GeomPanel,'Style','listbox', ...
                                  'BackgroundColor','w', ...
                                  'String','', ...
                                  'Callback',@GeomCallback, ...
                                  'KeyPressFcn',@GeomKeyPress, ...
                                  'UIContextMenu',[], ...
                                  'Min',0,'Max',2, ...
                                  'Units','pixels', ...
                                  'Position',[8 40 PanelsWidth-18 100]);

%% context menu for listbox
GeomContextMenu = uicontextmenu;
    GeomMenuInsert = uimenu(GeomContextMenu,'Label','Insert');
        h = copyobj(get(MenuInsert,'Children'),GeomMenuInsert);
        set(h,'Enable','on');
    GeomMenuEdit = uimenu(GeomContextMenu,'Label','Edit...', ...
                                          'Callback',@GeomCallback, ...
                                          'Separator','on', ...
                                          'Enable','off');
    GeomMenuDuplicate = uimenu(GeomContextMenu,'Label','Duplicate...', ...
                                               'Callback',@GeomCallback, ...
                                               'Separator','off', ...
                                               'Enable','off');
    GeomMenuTranslate = uimenu(GeomContextMenu,'Label','Translate...', ...
                                               'Callback',@GeomCallback, ...
                                               'Separator','off', ...
                                               'Enable','off');
    GeomMenuDelete = uimenu(GeomContextMenu,'Label','Delete...', ...
                                            'Callback',@GeomCallback, ...
                                            'Separator','on', ...
                                            'Enable','off');
    GeomMenuInvisible = uimenu(GeomContextMenu,'Label','Invisible', ...
                                               'Callback',@GeomCallback, ...
                                               'Separator','on', ...
                                               'Enable','off');
    GeomMenuHighlight = uimenu(GeomContextMenu,'Label','Highlight Off', ...
                                               'Callback',@GeomCallback, ...
                                               'Separator','on', ...
                                               'Enable','on');
    GeomMenuMoveUp   = uimenu(GeomContextMenu,'Label','Move Up', ...
                                              'Callback',@GeomCallback, ...
                                              'Separator','on', ...
                                              'Enable','off');
    GeomMenuMoveDown = uimenu(GeomContextMenu,'Label','Move Down', ...
                                              'Callback',@GeomCallback, ...
                                              'Separator','off', ...
                                              'Enable','off');
                                           
%% Highlighting is disabled by default
Highlight = false;
                                        
%% context menu for Insert pushbutton
GeomInsertContextMenu = uicontextmenu;
    h = copyobj(get(MenuInsert,'Children'),GeomInsertContextMenu);
    set(h,'Enable','on');
    
%% pushbuttons                                        
GeomPushInsert = uipushbutton(GeomPanel,'Insert...',[  8 10 75 23], ...
                              @GeomInsertCallback,'Enable','off');
GeomPushEdit   = uipushbutton(GeomPanel,'Edit...',  [ 88 10 75 23], ...
                              @GeomCallback,'Enable','off');
GeomPushDelete = uipushbutton(GeomPanel,'Delete...',[168 10 75 23], ...
                              @GeomCallback,'Enable','off');

%% set where the Insert context menu should appear                          
set([GeomPanel,GeomPushInsert],'Units','pixels');
ShowContext = get(GeomPanel,'Position') + get(GeomPushInsert,'Position');
    
%% showing the menu under the Insert pushbutton
function GeomInsertCallback(src,eventdata)
    set(GeomInsertContextMenu,'Position',ShowContext(1:2), ...
                              'Visible','on');
end
 
%% callback of the listbox and associated controls
function GeomCallback(src,eventdata)
    
    GeomPointer = get(GeomListBox,'Value');
    
    %% unselect all
    if Highlight && ishandle(GeomAxes), set(findobj(GeomAxes),'Selected','off'); end;
    
    %% if nothing selected, DUPLICATE, TRANSLATE, DELETE, INVISIBLE, 
    %% MOVE UP and MOVE DOWN are not allowed
    if isempty(GeomPointer),
        set([GeomMenuDuplicate,GeomMenuTranslate, ...
             GeomMenuDelete,GeomPushDelete,GeomMenuInvisible, ...
             GeomMenuMoveUp,GeomMenuMoveDown],'Enable','off');
    else
        %% in case file is locked, some are not allowed either
        if v.locked_flag,
            set([GeomMenuInsert,GeomMenuDuplicate,GeomMenuTranslate, ...
                 GeomMenuDelete,GeomMenuMoveUp,GeomMenuMoveDown],'Enable','off');
            set(GeomMenuInvisible,'Enable','on');
        else 
            set([GeomMenuInsert,GeomMenuDuplicate,GeomMenuTranslate, ...
                 GeomMenuDelete,GeomMenuMoveUp,GeomMenuMoveDown, ...
                 GeomMenuInvisible,GeomPushDelete],'Enable','on');
        end;
        
        %% visibility is taken from majority of multiple selected items
        visible = 0;
        for pointer = GeomPointer,
            %% querying visibility and making average
            visible = visible + v.Geometry{pointer}.visible;
            %% selecting objects in Geometry window
            if Highlight && all(ishandle(v.Geometry{pointer}.handle)),
                set(v.Geometry{pointer}.handle,'Selected','on');
            end;
        end;
        visible = round(visible/length(GeomPointer));
        if visible, set(GeomMenuInvisible,'Checked','off');
        else        set(GeomMenuInvisible,'Checked','on');
        end;
        
    end;

    %% if nothing or multiple items selected, EDIT is not allowed
    if isscalar(GeomPointer),
        set([GeomMenuEdit,GeomPushEdit],'Enable','on');
    else
        set([GeomMenuEdit,GeomPushEdit],'Enable','off');
    end;

    %% Editing geometry item
    if ( strcmp(get(MainFig,'SelectionType'),'open') ...
      || (src==GeomMenuEdit) ...
      || (src==GeomPushEdit) ) ...
        && isscalar(GeomPointer),
        Type = v.Geometry{GeomPointer}.type;
        switch Type
            case 'wire',         Wire(src,eventdata,GeomPointer);
            case 'triangle',     Triangle(src,eventdata,GeomPointer);
            case 'rectangle',    Rectangle(src,eventdata,GeomPointer);
            case 'brick',        Brick(src,eventdata,GeomPointer);
            case 'capacitor',    Capacitor(src,eventdata,GeomPointer);
            case 'resistor',     Resistor(src,eventdata,GeomPointer);
            case 'source',       Source(src,eventdata,GeomPointer);
            case 'field source', FieldSource(src,eventdata,GeomPointer);
            case 'probe',        Probe(src,eventdata,GeomPointer);
            case 'extobj',       ExtObj(src,eventdata,GeomPointer);
            case 'volobj',       ExtObj(src,eventdata,GeomPointer); %% obsolete
        end;
        return;
    end;
    
    %% Duplicating geometry items
    if (src==GeomMenuDuplicate),
        Duplicate(src,eventdata,GeomPointer);
        return;
    end;
    
    %% Translating geometry items
    if (src==GeomMenuTranslate),
        Translate(src,eventdata,GeomPointer);
        return;
    end;
    
    %% Deleting geometry items
    if (src==GeomMenuDelete) || (src==GeomPushDelete),
        Button = questdlg('Delete selected item(s)?','Delete Item(s)','No');
        if strcmp(Button,'Yes'),
            for pointer = GeomPointer,
                delete(v.Geometry{pointer}.handle);
            end;
            v.Geometry(GeomPointer) = [];
            v.GeometryString(GeomPointer,:) = [];
            
            %% ListboxTop should point to a valid component
            ListboxTopValue = get(GeomListBox,'ListboxTop');
            ListLength = size(v.GeometryString,1);
            if ListboxTopValue > ListLength,
                ListboxTopValue = ListLength;
                if ListboxTopValue == 0; ListboxTopValue = 1; end;
            end;
            set(GeomListBox,'String',v.GeometryString, ...
                            'Value',[], ...
                            'ListboxTop',ListboxTopValue);
                        
            set([GeomMenuEdit,GeomMenuDelete, ...
                 GeomPushEdit,GeomPushDelete],'Enable','off');
        end;
        return;
    end;

    %% Making geometry items invisible
    if (src==GeomMenuInvisible),
        if visible, %% if majority is visible, make invisible
            for pointer = GeomPointer,
                if v.Geometry{pointer}.visible, %% working only on visible
                    v.Geometry{pointer}.visible = false;
                    set(v.Geometry{pointer}.handle,'Visible','off');
                    temp = v.GeometryString(pointer,:);
                    temp = [InvisibilityChar temp];
                    %% including line
                    v.GeometryString = strvcat( ...
                        v.GeometryString(1:pointer-1,:), ...
                        temp, ...
                        v.GeometryString(pointer+1:end,:) );
                end;
            end;
            set(GeomMenuInvisible,'Checked','on');
        else %% if majority is invisible, make visible
            for pointer = GeomPointer,
                if ~v.Geometry{pointer}.visible, %% working only on invisible
                    v.Geometry{pointer}.visible = true;
                    set(v.Geometry{pointer}.handle,'Visible','on');
                    temp = v.GeometryString(pointer,:);
                    temp(1:length(InvisibilityChar)) = [];
                    %% including line
                    v.GeometryString = strvcat( ...
                        v.GeometryString(1:pointer-1,:), ...
                        temp, ...
                        v.GeometryString(pointer+1:end,:) );
                end;
            end;
            set(GeomMenuInvisible,'Checked','off');
        end;
        %% geometry string update
        set(GeomListBox,'String',v.GeometryString);
    end;

    %% Toggle highlighting of objects in the Geometry window
    if (src==GeomMenuHighlight),
        Highlight = ~Highlight;
        if Highlight,
            set(GeomMenuHighlight,'Label','Highlight On','Checked','on');
            if ishandle(GeomAxes),
                set(findobj(GeomAxes),'SelectionHighlight','on', ...
                                      'Selected','off');
            end;
            %% selecting objects in Geometry window
            for pointer = GeomPointer,
                if ishandle(v.Geometry{pointer}.handle),
                    set(v.Geometry{pointer}.handle,'Selected','on');
                end;
            end;
        else
            set(GeomMenuHighlight,'Label','Highlight Off','Checked','off');
            if ishandle(GeomAxes),
                set(findobj(GeomAxes),'SelectionHighlight','off');
            end;
        end;
    end;
    
    %% Move geometry items up
    if (src==GeomMenuMoveUp),
        %% skipping items at the beginning
        start = find(diff([0 GeomPointer])>1,1);
        for i = GeomPointer(start:end), %% for selected items...
            %% ...swap the preceding and actual items
            v.Geometry([i-1,i]) = v.Geometry([i,i-1]);
            v.GeometryString([i-1,i],:) = v.GeometryString([i,i-1],:);
        end;
        %% geometry string and highlighted item update
        GeomPointer(start:end) = GeomPointer(start:end) - 1;
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);
    end;
    
    %% Move geometry items down
    if (src==GeomMenuMoveDown),
        %% skipping items at the end
        stop = find(diff([GeomPointer length(v.Geometry)+1])>1,1,'last');
        for i = GeomPointer(stop:-1:1), %% for selected items...
            %% ...swap the succeeding and actual items
            v.Geometry([i+1,i]) = v.Geometry([i,i+1]);
            v.GeometryString([i+1,i],:) = v.GeometryString([i,i+1],:);
        end;
        %% geometry string and highlighted item update
        GeomPointer(1:stop) = GeomPointer(1:stop) + 1;
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);
    end;

end

%% listbox keypress callback
function GeomKeyPress(src,eventdata)
    if isequal(eventdata.Modifier,{'alt'}),
        if     isequal(eventdata.Key,'uparrow'),
            GeomCallback(GeomMenuMoveUp,eventdata);
        elseif isequal(eventdata.Key,'downarrow'),
            GeomCallback(GeomMenuMoveDown,eventdata);
        end;
    end;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD PARAMETERS PANEL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
PanelHeight = 110;
VerticalPosition = VerticalPosition - PanelHeight - Separate;
ParamPanel = uipanel('Parent',Panel, ...
                     'ForegroundColor','b', ...
                     'Title','FDTD parameters', ...
                     'Units','pixels', ...
                     'Position',[PanelsLeft VerticalPosition ...
                                 PanelsWidth PanelHeight]);

uitext(ParamPanel,'Cell size:'           ,[10 70 130 17]);
uitext(ParamPanel,'Frequency resolution:',[10 45 130 17]);

Dx_Edit = uiedit(ParamPanel,'',[140 70 80 21],'Tag','Dx');
Df_Edit = uiedit(ParamPanel,'',[140 45 80 21],'Tag','Df');

uitext(ParamPanel,'  m', [220 70 20 17]);
uitext(ParamPanel,'  Hz',[220 45 20 17]);

set([Dx_Edit Df_Edit],'Callback',@ParamCallback);

function ParamCallback(src,eventdata)
    
    if ~CheckEditboxValue(src,[]), return; end;
    
    v.(get(src,'Tag')) = str2num( get(src,'String') );

    if src==Dx_Edit && isequal(v.Dt_mode,'auto'),
        v.Dt = v.Dx/(v.c*sqrt(3))*v.CFL_coef;
        
        %% refresh the excitation pulse -- Dt has changed
        if isequal(v.ExcitationPulse,'gauss_sine'),
            switch v.GaussSineShape
                case 'auto'
                    [fc,fs,offset] = ComputeGaussSinePars( ...
                        v.fmin,v.fmax,v.td_thres,v.fd_thres);
                    v.Vs = GaussSinePulse(v.mag,fc,fs,offset);
                case 'manual'
                    v.Vs = GaussSinePulse(v.mag,v.fc,v.fs,v.offset);
            end;
        end;

    end;
    
end

ExcitePush = uipushbutton(ParamPanel,'Excitation...',[10 10 75 23], ...
                          @ExciteCallback);

FDTDAdvancedPush = uipushbutton(ParamPanel,'Advanced...',[165 10 75 23], ...
    @FDTDAdvancedCallback);

function ExciteCallback(src,eventdata)
    
    %% dialog window
    ExciteDialog = dialog('Visible','off', ...
        'Name','Excitation Pulse', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[450 250] 900 600], ...
        'KeyPressFcn',{@RetEscPressed,@ExciteStore,@CloseFigure});

    
    %% Gauss weighted sine pulse panel
    GaussSinePulsePanel = uipanel( ...
        'Parent',ExciteDialog, ...
        'ForegroundColor','b', ...
        'Units','pixels', ...
        'Position',[7 45 590 170]);

        GaussSinePulseRadio = uiradiobutton( ...
            GaussSinePulsePanel,0,[10 161 147 15], ...
            'ForegroundColor','b', ...
            'String','Gauss weighted sine pulse', ...
            'Callback',@ExciteRadiobuttonCallback);
        
        %% magnitude
        uitext(GaussSinePulsePanel,'Magnitude:',[10 130 70 17]);
        uiedit(GaussSinePulsePanel,'',[80 130 80 21],'Tag','mag');
        uitext(GaussSinePulsePanel,'  V',[160 130 25 17]);
        uipushbutton(GaussSinePulsePanel,'>',[318 53 23 23],@GaussSineRecompute);
        
        FormulaAxes = axes( ...
            'Parent',GaussSinePulsePanel, ...
            'Units','pixels', ...
            'Position',[270 142 1 1], ...
            'Visible','off');
        
        FormulaText = text( ...
            'Parent',FormulaAxes, ...
            'Position',[0 0], ...
            'Interpreter','latex', ...
            'String',['$f(t) = M\cdot\sin[2{\pi}{f_C}(t-\tau)]\cdot ', ...
                      'e^{-[f_S(t-\tau)]^2}$'], ...
            'FontUnits','pixels', ...
            'FontSize',13);
        
        %% automatic shape
        GaussSineAutoPanel = uipanel( ...
            'Parent',GaussSinePulsePanel, ...
            'Units','pixels', ...
            'Position',[10 20 295 90], ...
            'BorderType','none');
        
            uipanel('Parent',GaussSineAutoPanel, ...
                    'Units','pixels', ...
                    'Position',[0 89 295 2]); %% horizontal line
        
            GaussSineAutoRadio = uiradiobutton( ...
                GaussSineAutoPanel,0,[10 83 100 15], ...
                'String','Automatic shape', ...
                'Callback',@GaussSineRadioCallback);
            
            uitext(GaussSineAutoPanel,'Frequency range:',[0 55 100 17]);
            uitext(GaussSineAutoPanel,'--',[180 55 15 17], ...
                'HorizontalAlignment','center');
            uitext(GaussSineAutoPanel,'Frequency domain threshold:',[0 25 195 17]);
            uitext(GaussSineAutoPanel,'Time domain threshold:',[0 0 195 17]);
    
            uiedit(GaussSineAutoPanel,'',[100 55 80 21],'Tag','fmin');
            uiedit(GaussSineAutoPanel,'',[195 55 80 21],'Tag','fmax');
            uiedit(GaussSineAutoPanel,'',[195 25 80 21],'Tag','fd_thres');
            uiedit(GaussSineAutoPanel,'',[195  0 80 21],'Tag','td_thres');

            uitext(GaussSineAutoPanel,'  Hz',[275 55 20 17]);
            uitext(GaussSineAutoPanel,'  dB',[275 25 20 17]);
            uitext(GaussSineAutoPanel,'  dB',[275  0 20 17]);
            

        %% manual shape
        GaussSineManualPanel = uipanel( ...
            'Parent',GaussSinePulsePanel, ...
            'Units','pixels', ...
            'Position',[355 20 220 90], ...
            'BorderType','none');

            uipanel('Parent',GaussSineManualPanel, ...
                    'Units','pixels', ...
                    'Position',[0 89 220 2]); %% horizontal line
    
            GaussSineManualRadio = uiradiobutton( ...
                GaussSineManualPanel,0,[10 83 87 15], ...
                'String','Manual shape', ...
                'Callback',@GaussSineRadioCallback);

            uitext(GaussSineManualPanel,'Center frequency:'     ,[0 55 120 17]);
            uitext(GaussSineManualPanel,'Frequency spread:'     ,[0 30 120 17]);
            uitext(GaussSineManualPanel,'Impulse center offset:',[0  5 120 17]);

            uiedit(GaussSineManualPanel,'',[120 55 80 21],'Tag','fc');
            uiedit(GaussSineManualPanel,'',[120 30 80 21],'Tag','fs');
            uiedit(GaussSineManualPanel,'',[120  5 80 21],'Tag','offset');

            uitext(GaussSineManualPanel,'  Hz',[200 55 20 17]);
            uitext(GaussSineManualPanel,'  Hz',[200 30 20 17]);
            uitext(GaussSineManualPanel,'  s' ,[200  5 20 17]);
            

    %% Custom pulse panel
    CustomPulsePanel = uipanel( ...
        'Parent',ExciteDialog, ...
        'ForegroundColor','b', ...
        'Units','pixels', ...
        'Position',[605 45 288 170]);
                      
        CustomPulseRadio = uiradiobutton( ...
            CustomPulsePanel,0,[10 161 86 15], ...
            'ForegroundColor','b', ...
            'String','Custom pulse', ...
            'Callback',@ExciteRadiobuttonCallback);
        
        CustomPulseLoadButton = uipushbutton( ...
            CustomPulsePanel,'Load...',[10 130 75 23],@CustomPulseLoad);
        
            function CustomPulseLoad(src,eventdata)
                
                while 1,
                
                    %% Show dialog window for opening a file
                    [FileName,PathName] = uigetfile( ...
                        {'*.mat','MAT-files (*.mat)'; ...
                         '*','All Files (*.*)'}, ...
                        'Load Custom Pulse',[WorkDir,filesep]);

                    %% If nothing selected or Cancel, return
                    if isequal(FileName,0) || isequal(PathName,0), return; end;
                    WorkDir = PathName;

                    %% If opening fails, go back the dialog window 
                    try
                        VarList = who('-file',fullfile(PathName,FileName));
                    catch
                        errordlg(['Error when opening file ',FileName], ...
                            'Open Error','modal');
                        uiwait;
                        continue;
                    end;
                    
                    %% otherwise everything is ok, and leave the loop
                    break;
                
                end;
                
                while 1,
                    
                    %% Variable selection dialog window
                    [VarPtr,ok] = listdlg( ...
                        'ListString',VarList, ...
                        'SelectionMode','single', ...
                        'ListSize',[160 160], ...
                        'Name','Load Variable', ...
                        'PromptString','Select variable to load:');
                    if ~ok, return; end; %% If Cancel pressed, return
                    VarName = VarList{VarPtr};
                    
                    %% If opening fails, go back to the variable selection
                    try
                        s = load('-mat',fullfile(PathName,FileName),VarName);
                    catch
                        errordlg(['Error when opening file ',FileName], ...
                            'Open Error','modal');
                        uiwait;
                        continue;
                    end;
                    
                    LoadedVar = s.(VarName);

                    %% If the variable is not usable, go back to the
                    %% variable selection window
                    if     ~isvector(LoadedVar),      ErrorString = 'a vector.';
                    elseif ~isnumeric(LoadedVar),     ErrorString = 'a numeric array.';
                    elseif ~all(isfinite(LoadedVar)), ErrorString = 'finite.';
                    elseif ~isreal(LoadedVar),        ErrorString = 'real.';
                    else break; %% otherwise leave the loop
                    end;
                    errordlg({['Error when loading variable ''',VarName,''':'], ...
                        ['Variable is not ',ErrorString]},'Load Error','modal');
                    uiwait;
                    
                end;
                
                temp.CustomPulse = LoadedVar(:);
                
            end
    
        
    %% control handles for enabling process
    GaussSinePulseControls = [ ...
        findobj(GaussSinePulsePanel,'-depth',1,'Style','text','-or', ...
            'Style','edit','-or','Style','pushbutton'); ...
        GaussSineAutoRadio; GaussSineManualRadio ];
        GaussSineAutoControls = findobj(GaussSineAutoPanel, ...
            'Style','text','-or','Style','edit');
        GaussSineManualControls = findobj(GaussSineManualPanel, ...
            'Style','text','-or','Style','edit');
    CustomPulseControls = CustomPulseLoadButton;
       
    %% pulse selection radiobutton callbacks
    function ExciteRadiobuttonCallback(src,eventdata)
        
        switch src
            
            case GaussSinePulseRadio
                temp.ExcitationPulse = 'gauss_sine';
                set(GaussSinePulseRadio,'Value',1);
                set(CustomPulseRadio   ,'Value',0);
                set(GaussSinePulseControls,'Enable','on' );
                switch temp.GaussSineShape
                    case 'auto',   set(GaussSineAutoControls  ,'Enable','on');
                    case 'manual', set(GaussSineManualControls,'Enable','on');
                end;
                set(CustomPulseControls,'Enable','off');
                set(FormulaText,'Color','k');
                
            case CustomPulseRadio
                temp.ExcitationPulse = 'custom';
                set(GaussSinePulseRadio,'Value',0);
                set(CustomPulseRadio   ,'Value',1);
                set([ GaussSinePulseControls  ; ...
                      GaussSineAutoControls   ; ...
                      GaussSineManualControls ],'Enable','off');
                set(CustomPulseControls,'Enable','on' );
                set(FormulaText,'Color',[.6 .6 .6]);
                
        end;
        
    end
        
    %% automatic - manual pulse shape radiobutton callback
    function GaussSineRadioCallback(src,eventdata)
        
        switch src
            
            case GaussSineAutoRadio
                temp.GaussSineShape = 'auto';
                set(GaussSineAutoRadio  ,'Value',1);
                set(GaussSineManualRadio,'Value',0);
                set(GaussSineAutoControls  ,'Enable','on' );
                set(GaussSineManualControls,'Enable','off');
                
            case GaussSineManualRadio
                temp.GaussSineShape = 'manual';
                set(GaussSineAutoRadio  ,'Value',0);
                set(GaussSineManualRadio,'Value',1);
                set(GaussSineAutoControls  ,'Enable','off');
                set(GaussSineManualControls,'Enable','on' );
                
        end;
        
    end
    
    %% filling the temporary variables
    VarNames = {'mag','fmin','fmax','td_thres','fd_thres', ...
                'fc','fs','offset','GaussSineShape','ExcitationPulse','Vs'};
    for i = 1:length(VarNames),
        temp.(VarNames{i}) = v.(VarNames{i});
    end;
    
    
    %% filling the states of the radiobuttons
    switch temp.GaussSineShape,
        case 'auto',   GaussSineRadioCallback(GaussSineAutoRadio);
        case 'manual', GaussSineRadioCallback(GaussSineManualRadio);
    end;
    switch temp.ExcitationPulse,
        case 'gauss_sine'
            ExciteRadiobuttonCallback(GaussSinePulseRadio);
            temp.CustomPulse = [];
        case 'custom'
            ExciteRadiobuttonCallback(CustomPulseRadio);
            temp.CustomPulse = temp.Vs;
    end;
    
   
    AxesPanel = uipanel( ...
        'Parent',ExciteDialog, ...
        'Units','pixels', ...
        'Position',[7 230 886 363]);
    
    %% time-domain plot
    TD_Axes = axes( ...
        'Parent',AxesPanel, ...
        'Units','pixels', ...
        'OuterPosition',[0 50 443 300], ...
        'ActivePositionProperty','outerposition', ...
        'Units','normalized', ...
        'Box','on');
    
    %% frequency-domain plot
    FD_Axes = axes( ...
        'Parent',AxesPanel, ...
        'Units','pixels', ...
        'OuterPosition',[443 50 443 300], ...
        'ActivePositionProperty','outerposition', ...
        'Units','normalized', ...
        'Box','on');
    
    uipushbutton(AxesPanel,'Refresh',[10 12 75 23],@ExciteDisplayRefresh);
    uitogglebutton(AxesPanel,'Zoom' ,0,[ 94 12 75 23],@ExciteZoom);
    
    %% finding handles of the editboxes in the order as in VarNames
    for i = 1:8,
        Editboxes(i) = findobj(GaussSinePulsePanel,'Tag',VarNames{i});
    end;
    
    FillEditboxes(Editboxes);
    
    ExciteDisplayRefresh;
    
    %% default behavior of Enter and Escape on controls
    set(Editboxes,'KeyPressFcn',{@RetEscPressed,@ExciteStore,@CloseFigure});
    
    %% buttons
    OKButton     = Button(ExciteDialog,'OK'    ,2,@ExciteStore);
    CancelButton = Button(ExciteDialog,'Cancel',1,@CloseFigure);

    %% locking the uicontrols
    if v.locked_flag, LockDialog(ExciteDialog,OKButton,CancelButton); end;
    
    uicontrol(GaussSineAutoRadio); %% give focus to the first enabled control
    set(ExciteDialog,'Visible','on'); %% display figure

    %% callback in case of pressing Enter on editboxes or pushing OK -
    %% storing and exit
    function ExciteStore(src,eventdata);
    
        %% filling the temporary variables from the editboxes
        if ~FillTempVariables(Editboxes), return; end;
        
        RefreshExcitationPulse;
        
        %% filling the global variables from the temporary ones
        for i = 1:length(VarNames),
            v.(VarNames{i}) = temp.(VarNames{i});
        end;
        
        close(ExciteDialog);
        
    end
        
    %% recomputing the Gauss-sine manual values from the auto values
    function GaussSineRecompute(src,eventdata)

        %% filling the temporary variables from the editboxes
        if ~FillTempVariables(Editboxes(2:5)), return; end;
        
        %% recomputing
        [temp.fc,temp.fs,temp.offset] = ...
            ComputeGaussSinePars(temp.fmin,temp.fmax,temp.td_thres,temp.fd_thres);
        
        %% filling the editboxes from the temporary variables
        FillEditboxes(Editboxes(6:8));
        
    end
        
    %% refreshing the plots
    function ExciteDisplayRefresh(src,eventdata)
        
        if isequal(temp.ExcitationPulse,'gauss_sine'),
            
            %% choose the actual editboxes
            switch temp.GaussSineShape
                case 'auto',   ActualEditboxes = Editboxes(1:5);
                case 'manual', ActualEditboxes = Editboxes([1 6:8]);
            end;
            
            %% filling the temporary variables from the editboxes
            if ~FillTempVariables(ActualEditboxes), return; end;
            
        end;
            
        RefreshExcitationPulse;
        
        %% if Vs empty, make empty plots
        if isempty(temp.Vs),
            plot(TD_Axes,0,0);
            plot(FD_Axes,0,0);
            return;
        end;

        %% FFT

        %% time domain length with zero padding for particular v.Df
        Nt_fft = round( 1/(v.Dt*v.Df) );

        %% for too low Nt_fft the pulse would be cropped - avoid this
        if Nt_fft < length(temp.Vs), Nt_fft = length(temp.Vs); end;
        
        %% "real" Df taken from already rounded Nt_fft
        Df = 1/(v.Dt*Nt_fft);
        
        temp.Vs_fft = 2.*fft(temp.Vs,Nt_fft)./Nt_fft;
        temp.Vs_fft = temp.Vs_fft(1:floor(Nt_fft/2)+1);
        
        %% display
        t = (0:length(temp.Vs)-1).*v.Dt;
        plot(  TD_Axes,t,temp.Vs);
        title( TD_Axes,'Time domain');
        xlabel(TD_Axes,'time [s]');
        ylabel(TD_Axes,'voltage [V]');

        f = (0:length(temp.Vs_fft)-1).*Df;
        plot(  FD_Axes,f,abs(temp.Vs_fft));
        title( FD_Axes,'Frequency domain');
        xlabel(FD_Axes,'frequency [Hz]');
        ylabel(FD_Axes,'voltage [V]');
        
    end
    
    %% zooming button callback
    function ExciteZoom(src,eventdata)
        switch get(src,'Value')
            case 0, zoom off;
            case 1, zoom on;
        end;
    end

    %% check editboxes for valid entries and fill the temp. variables
    function status = FillTempVariables(handles)
        
        for h = handles,
            
            if any(strcmp(get(h,'Tag'),{'td_thres','fd_thres'})),
                 sign = 'negative';
            else sign = []; end;
            
            %% check the value if it is correct
            status = CheckEditboxValue(h,sign);
            
            %% if not, return with status=0
            if ~status, return; end;
            
            %% if yes, fill the temporary variable
            temp.(get(h,'Tag')) = str2num( get(h,'String') );
            
        end;
        
    end
    
    %% fill editboxes from the temp. variables
    function FillEditboxes(handles)
        for h = handles
            var_name = get(h,'Tag');
            switch var_name
                case {'td_thres','fd_thres'}, num_format = '%d';
                otherwise                     num_format = '%.4g';
            end;
            set(h,'String',num2str(temp.(var_name),num_format));
        end;
    end
    
    %% no comment...
    function RefreshExcitationPulse
        
        switch temp.ExcitationPulse
            case 'gauss_sine'
                switch temp.GaussSineShape
                    case 'auto'
                        [fc,fs,offset] = ComputeGaussSinePars( ...
                            temp.fmin,temp.fmax,temp.td_thres,temp.fd_thres);
                        temp.Vs = GaussSinePulse(temp.mag,fc,fs,offset);
                    case 'manual'
                        temp.Vs = GaussSinePulse(temp.mag,temp.fc,temp.fs,temp.offset);
                end;
            case 'custom'
                temp.Vs = temp.CustomPulse;
        end;
        
    end
    
end

function FDTDAdvancedCallback(src,eventdata)
    
    %% dialog window
    FDTDAdvancedDialog = dialog('Visible','off', ...
        'Name','Advanced FDTD Settings', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[200 200] 340 337], ...
        'KeyPressFcn',{@RetEscPressed,@FDTDAdvancedStore,@CloseFigure});
   
    %% Time step setting panel
    TimeStepPanel = uipanel( ...
        'Parent',FDTDAdvancedDialog, ...
        'ForegroundColor','b', ...
        'Title','Time step', ...
        'Units','pixels', ...
        'Position',[7 219 328 115]);
        
        TimeStepRadios = uiradiogroup(TimeStepPanel,[10 15 310 85], ...
            @TimeStepCallback);
        TimeOptions(1) = uiradiobutton(TimeStepRadios,0,[0 50 130 15], ...
            'String','From CFL condition:','Tag','auto');
        TimeOptions(2) = uiradiobutton(TimeStepRadios,0,[0 5 75 15], ...
            'String','Manual:','Tag','manual');
        
        FormulaAxes = axes( ...
            'Parent',TimeStepRadios, ...
            'Units','pixels', ...
            'Position',[130 52 1 1], ...
            'Visible','off');
        
        FormulaText = text( ...
            'Parent',FormulaAxes, ...
            'Position',[0 0], ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','baseline', ...
            'Interpreter','latex', ...
            'String',['$$\Delta t = \frac{\Delta x}{c\sqrt{3}} \times$$'], ...
            'FontUnits','pixels', ...
            'FontSize',13);
        
        CFL_Edit = uiedit(TimeStepRadios,'',[225 47 80 21],'Tag','CFL_coef');
        Dt_Edit  = uiedit(TimeStepRadios,'',[ 75  2 80 21],'Tag','Dt');

        TimeStepEdit = [CFL_Edit Dt_Edit];
        
        uitext(TimeStepRadios,'  s',[155 2 30 17]);
        
        ComputePush = uipushbutton(TimeStepRadios,'Compute from CFL', ...
            [185 2 120 23],@ComputeFromCFLCallback);
        
        function ComputeFromCFLCallback(src,eventdata)
            
            val = str2num( get(CFL_Edit,'String') );
            if isempty(val) || ~isfinite(val) || imag(val) ~= 0 || val <= 0,
                return;
            end;
            
            new_val = v.Dx/v.c/sqrt(3) * val;
            set(Dt_Edit,'String',num2str(new_val,'%.4g'));
            
        end
        
        function TimeStepCallback(src,eventdata)
            
            TimeStepSwitch = TimeOptions == get(src,'SelectedObject');
            set(TimeStepEdit( TimeStepSwitch),'enable','on' );
            set(TimeStepEdit(~TimeStepSwitch),'enable','off');
            if TimeStepSwitch(2),
                set(ComputePush,'Enable','on');
                set(FormulaText,'Color',[.6 .6 .6]);
            else
                set(ComputePush,'Enable','off');
                set(FormulaText,'Color','k');
            end;
            
        end
        
    %% Free-space properties panel
    SpacePanel = uipanel( ...
        'Parent',FDTDAdvancedDialog, ...
        'ForegroundColor','b', ...
        'Title','Free-space properties', ...
        'Units','pixels', ...
        'Position',[7 97 328 120]);

        uitext(SpacePanel,'Permittivity:'  ,[10 78 80 17]);
        uitext(SpacePanel,'Permeability:'  ,[10 53 80 17]);
        uitext(SpacePanel,'Speed of light:',[10 28 80 17]);

        Eps_Edit = uiedit(SpacePanel,'',[90 78 80 21],'Tag','eps0');
        Mu_Edit  = uiedit(SpacePanel,'',[90 53 80 21],'Tag','mu0');
        C_Edit   = uiedit(SpacePanel,'',[90 28 80 21],'Tag','c');

        uitext(SpacePanel,'  F/m',[170 78 30 17]);
        uitext(SpacePanel,'  H/m',[170 53 30 17]);
        uitext(SpacePanel,'  m/s',[170 28 30 17]);

        SpaceEdit = [Eps_Edit Mu_Edit C_Edit];
        
        SpaceDep = uiradiogroup(SpacePanel,[205 22 24 80],@SpaceDepSelect);
        SpaceDepButton(1) = uiradiobutton(SpaceDep,0,[5 60 15 15]);
        SpaceDepButton(2) = uiradiobutton(SpaceDep,0,[5 35 15 15]);
        SpaceDepButton(3) = uiradiobutton(SpaceDep,0,[5 10 15 15]);

        uiline(SpaceDep,[.6 .6 .6],[12 0 1 61]);
        uitext(SpacePanel,'Dependent',[187 4 60 17],'HorizontalAlignment','center');

        function SpaceDepSelect(src,eventdata)
            SpaceDepSwitch = SpaceDepButton == eventdata.NewValue;
            set(SpaceEdit( SpaceDepSwitch),'enable','off');
            set(SpaceEdit(~SpaceDepSwitch),'enable','on' );
        end
        
    %% Free-space properties panel
    GridPanel = uipanel( ...
        'Parent',FDTDAdvancedDialog, ...
        'ForegroundColor','b', ...
        'Title','Grid properties', ...
        'Units','pixels', ...
        'Position',[7 45 328 50]);

        MuAvgCheckbox = uicheckbox(GridPanel,'Permeability averaging', ...
            v.MuAveraging,[10 10 305 21],'');

    %% Activation
    Editboxes = [TimeStepEdit SpaceEdit];
    
    %% fill editboxes
    for h = Editboxes,
        set(h,'String',num2str(v.(get(h,'Tag')),'%.4g'));
    end;

    %% setting the free space radiobuttons and editboxes
    set(SpaceEdit(v.SpaceDepSwitch==1),'Enable','off');
    set(SpaceDep,'SelectedObject',SpaceDepButton(~~v.SpaceDepSwitch))

    %% setting the time step radiobuttons and editboxes
    set(TimeStepRadios,'SelectedObject',findobj(TimeStepRadios,'Tag',v.Dt_mode));
    TimeStepCallback(TimeStepRadios);

    %% default behavior of Enter and Escape on controls
    set(Editboxes,'KeyPressFcn',{@RetEscPressed,@FDTDAdvancedStore,@CloseFigure});
    
    %% buttons
    OKButton     = Button(FDTDAdvancedDialog,'OK'    ,2,@FDTDAdvancedStore);
    CancelButton = Button(FDTDAdvancedDialog,'Cancel',1,@CloseFigure);

    %% locking the uicontrols
    if v.locked_flag, LockDialog(FDTDAdvancedDialog,OKButton,CancelButton); end;
    
    uicontrol(TimeOptions(1)); %% give focus to the first enabled control
    set(FDTDAdvancedDialog,'Visible','on'); %% display figure
    
    %% callback in case of pressing Enter on editboxes or pushing OK -
    %% storing and exit
    function FDTDAdvancedStore(src,eventdata);

        %% checking for correct values
        for h = Editboxes,
            if ~CheckEditboxValue(h,[]), return; end;
        end;

        %% writing values
        for h = Editboxes,
            v.(get(h,'Tag')) = str2num( get(h,'String') );
        end;
        
        %% writing radiobutton switches
        v.Dt_mode = get(get(TimeStepRadios,'SelectedObject'),'Tag');
        v.SpaceDepSwitch = SpaceDepButton == get(SpaceDep,'SelectedObject');

        %% writing permeability averaging
        v.MuAveraging = get(MuAvgCheckbox,'Value');

        %% computing dependent variables
        switch find(v.SpaceDepSwitch),
            case 1, v.eps0 = 1/(v.mu0*v.c^2);      %% eps0 dependent
            case 2, v.mu0  = 1/(v.eps0*v.c^2);     %% mu0  dependent
            case 3, v.c    = 1/sqrt(v.eps0*v.mu0); %% c    dependent
        end;
        
        if isequal(v.Dt_mode,'auto'),
            v.Dt = v.Dx/(v.c*sqrt(3))*v.CFL_coef;
        end;

        %% refresh the excitation pulse
        if isequal(v.ExcitationPulse,'gauss_sine'),
            switch v.GaussSineShape
                case 'auto'
                    [fc,fs,offset] = ComputeGaussSinePars( ...
                        v.fmin,v.fmax,v.td_thres,v.fd_thres);
                    v.Vs = GaussSinePulse(v.mag,fc,fs,offset);
                case 'manual'
                    v.Vs = GaussSinePulse(v.mag,v.fc,v.fs,v.offset);
            end;
        end;

        close(FDTDAdvancedDialog);
        
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TERMINATION CONDITION PANEL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PanelHeight = 100;
VerticalPosition = VerticalPosition - PanelHeight - Separate;
TermPanel = uipanel('Parent',Panel, ...
                    'ForegroundColor','b', ...
                    'Title','Termination condition', ...
                    'Units','pixels', ...
                    'Position',[PanelsLeft VerticalPosition ...
                                PanelsWidth PanelHeight]);

uitext(TermPanel,'Maximum number of time steps:',[10 60 165 17]);
Nt_Edit = uiedit(TermPanel,'',[175  60 50 21]);

DropCheckbox = uicheckbox(TermPanel,'Total energy below:',0,[10 10 165 21], ...
           @DropCheckCallback);
Drop_Edit = uiedit(TermPanel,'',[175 10 50 21]);
uitext(TermPanel,'  dB',[225 10 20 17]);

Drop2Checkbox = uicheckbox(TermPanel,'Source energy below:',0,[10 35 165 21], ...
           @Drop2CheckCallback);
Drop2_Edit = uiedit(TermPanel,'',[175 35 50 21]);
uitext(TermPanel,'  dB',[225 35 20 17]);

set([Nt_Edit Drop_Edit Drop2_Edit],'Callback',@TermCallback);

function DropCheckCallback(src,eventdata)
    v.drop_flag = get(src,'Value')==1;
    if v.drop_flag, set(Drop_Edit,'Enable','on');
    else            set(Drop_Edit,'Enable','off');
    end;
end

function Drop2CheckCallback(src,eventdata)
    v.drop2_flag = get(src,'Value')==1;
    if v.drop2_flag, set(Drop2_Edit,'Enable','on');
    else             set(Drop2_Edit,'Enable','off');
    end;
end

function TermCallback(src,eventdata)
    switch src,
        case Nt_Edit,
            if CheckEditboxValue(src,[]),
                v.Nt   = str2num( get(src,'String') );
            end;
        case Drop_Edit,
            if CheckEditboxValue(src,'negative'),
                v.Drop_dB = str2num( get(src,'String') );
            end;
        case Drop2_Edit,
            if CheckEditboxValue(src,'negative'),
                v.Drop2_dB = str2num( get(src,'String') );
            end;
    end;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOMAIN BOUNDARIES PANEL %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PanelHeight = 120;
VerticalPosition = VerticalPosition - PanelHeight - Separate;
DomainPanel = uipanel('Parent',Panel, ...
                      'ForegroundColor','b', ...
                      'Title','Domain boundaries', ...
                      'Units','pixels', ...
                      'Position',[PanelsLeft VerticalPosition ...
                                  PanelsWidth PanelHeight]);

uitext(DomainPanel,'x',[ 50 86 60 16],'HorizontalAlignment','center');
uitext(DomainPanel,'y',[115 86 60 16],'HorizontalAlignment','center');
uitext(DomainPanel,'z',[180 86 60 16],'HorizontalAlignment','center');

uitext(DomainPanel,'lower:',[10 65 40 17]);
uitext(DomainPanel,'upper:',[10 40 40 17]);

Domain_xlim1 = uiedit(DomainPanel,'',[ 50 65 60 21]);
Domain_ylim1 = uiedit(DomainPanel,'',[115 65 60 21]);
Domain_zlim1 = uiedit(DomainPanel,'',[180 65 60 21]);
Domain_xlim2 = uiedit(DomainPanel,'',[ 50 40 60 21]);
Domain_ylim2 = uiedit(DomainPanel,'',[115 40 60 21]);
Domain_zlim2 = uiedit(DomainPanel,'',[180 40 60 21]);

set([ Domain_xlim1 Domain_ylim1 Domain_zlim1 ...
      Domain_xlim2 Domain_ylim2 Domain_zlim2 ],'Callback',@DomainCallback);

function DomainCallback(src,eventdata)
    
    value = str2num( get(src,'String') );
    if isempty(value) || ~isfinite(value) || ...
       imag(value)~=0 || value~=round(value),
        errordlg('Domain limits must be integer numeric values', ...
                 'Bad input','modal');
        pause(1);
        uicontrol(src);
        return;
    end;

    temp_limits = v.limits;
    switch src,
        case Domain_xlim1, temp_limits(1,1) = value;
        case Domain_ylim1, temp_limits(1,2) = value;
        case Domain_zlim1, temp_limits(1,3) = value;
        case Domain_xlim2, temp_limits(2,1) = value;
        case Domain_ylim2, temp_limits(2,2) = value;
        case Domain_zlim2, temp_limits(2,3) = value;
    end;
    
    if any( temp_limits(1,:) >= temp_limits(2,:) ),
        errordlg('Upper limit must be higher than lower limit', ...
                 'Bad input','modal');
        pause(1);
        uicontrol(src);
        return;
    end;
        
    v.limits = temp_limits;

    %% setting the new limits ('axis' didn't work as it resets the x-y view)
    xlim(GeomAxes,v.limits(:,1).');
    ylim(GeomAxes,v.limits(:,2).');
    zlim(GeomAxes,v.limits(:,3).');

    %% changing the axis limits even for rotate reset view
    %% (otherwise 'Reset to Original View' will revoke the original limits)
    if isappdata(GeomAxes,'matlab_graphics_resetplotview')
        ResetViewStruct = getappdata(GeomAxes,'matlab_graphics_resetplotview');
        ResetViewStruct.XLim = v.limits(:,1).';
        ResetViewStruct.YLim = v.limits(:,2).';
        ResetViewStruct.ZLim = v.limits(:,3).';
        setappdata(GeomAxes,'matlab_graphics_resetplotview',ResetViewStruct);
    end;

end

%% PML setup
PML_Checkbox = uicheckbox(DomainPanel,'Perfectly matched layers',0,[10 10 150 21], ...
           @PML_CheckCallback);
function PML_CheckCallback(src,eventdata)
    v.pml = get(src,'Value');
    if v.pml==1, set(PML_SettingsPush,'Enable','on');
    else         set(PML_SettingsPush,'Enable','off');
    end;
end

PML_SettingsPush = uipushbutton(DomainPanel,'Settings...',[165 10 75 23], ...
                                @PML_SettingsCallback);
function PML_SettingsCallback(src,eventdata)
    PML_Dialog = dialog('Visible','off', ...
        'Name','PML Settings', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 70] 375 210], ...
        'KeyPressFcn',{@RetEscPressed,@PML_Store,@CloseFigure});

    uitext(PML_Dialog,'PML depth:',[10 180 70 17]);
    pmax_Edit = uiedit(PML_Dialog,num2str(v.pmax,'%d'),[80 180 80 21]);
    uitext(PML_Dialog,'  cells',[160 180 40 17]);
    
    PML_ParamPanel = uipanel('Parent',PML_Dialog, ...
                             'ForegroundColor','b', ...
                             'Title','PML parameters', ...
                             'Units','pixels', ...
                             'Position',[10 45 355 120]);
    
    uitext(PML_ParamPanel,'sigma',[ 90 86 80 16],'HorizontalAlignment','center');
    uitext(PML_ParamPanel,'kappa',[175 86 80 16],'HorizontalAlignment','center');
    uitext(PML_ParamPanel,'a',    [260 86 80 16],'HorizontalAlignment','center');
                         
    uitext(PML_ParamPanel,'Grading order:',[10 65 80 17]);
    uitext(PML_ParamPanel,'Maximum:',      [10 40 80 17]);

    ms_Edit       = uiedit(PML_ParamPanel,num2str(v.m_s,'%.4g'),[ 90 65 80 21]);
    mk_Edit       = uiedit(PML_ParamPanel,num2str(v.m_k,'%.4g'),[175 65 80 21]);
    ma_Edit       = uiedit(PML_ParamPanel,num2str(v.m_a,'%.4g'),[260 65 80 21]);
    sigmax_Edit   = uiedit(PML_ParamPanel,num2str(v.sig_max,  '%.4g'),[ 90 40 80 21]);
    kappamax_Edit = uiedit(PML_ParamPanel,num2str(v.kappa_max,'%.4g'),[175 40 80 21]);
    amax_Edit     = uiedit(PML_ParamPanel,num2str(v.a_max,    '%.4g'),[260 40 80 21]);

    PML_SigmaAutoPush = uipushbutton(PML_ParamPanel,'Auto', ...
                                     [90 10 75 23],@PML_SigmaAutoCallback);
    function PML_SigmaAutoCallback(src,eventdata)
        value = str2num( get(ms_Edit,'String') );
        if isempty(value) || ~isfinite(value) || imag(value)~=0 || value<=0,
            errordlg('Sigma grading order must be real numeric value', ...
                'Bad input','modal');
            pause(1);
            uicontrol(ms_Edit);
            return;
        end;
        value = (value + 1)/150/pi/v.Dx;
        set(sigmax_Edit,'String',num2str(value,'%.4g'));
        %% suggesting the a parameter (see p. 220, 245):
        value = 0.2 * 0.001 / v.Dx;
        set(amax_Edit,'String',num2str(value,'%.4g'));
    end
    
    PML_Edit = [pmax_Edit ms_Edit mk_Edit ma_Edit ...
                sigmax_Edit kappamax_Edit amax_Edit];

    %% default behavior of Enter and Escape on controls
    set(PML_Edit,'KeyPressFcn',{@RetEscPressed,@PML_Store,@CloseFigure});

    %% buttons
    OKButton = Button(PML_Dialog,'OK',2,@PML_Store);
    CancelButton = Button(PML_Dialog,'Cancel',1,@CloseFigure);

    %% locking the uicontrols
    if v.locked_flag, LockDialog(PML_Dialog,OKButton,CancelButton); end;
    
    uicontrol(pmax_Edit); %% give focus to the first field
    set(PML_Dialog,'Visible','on'); %% display figure

    function PML_Store(src,eventdata);
        for h = PML_Edit,
            value = str2num( get(h,'String') );
            if isempty(value) || ~isfinite(value) || imag(value)~=0 || value<0,
                errordlg('PML parameters must be real numeric values', ...
                    'Bad input','modal');
                pause(1);
                uicontrol(h);
                return;
            end;
            switch h,
                case pmax_Edit, v.pmax = value;
                case ms_Edit, v.m_s = value;
                case mk_Edit, v.m_k = value;
                case ma_Edit, v.m_a = value;
                case sigmax_Edit,   v.sig_max   = value;
                case kappamax_Edit, v.kappa_max = value;
                case amax_Edit,     v.a_max     = value;
            end;
        end;
        close(PML_Dialog);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISCELLANEOUS PANEL %%
%%%%%%%%%%%%%%%%%%%%%%%%%
PanelHeight = 100;
VerticalPosition = VerticalPosition - PanelHeight - Separate;
MiscPanel = uipanel('Parent',Panel, ...
                    'ForegroundColor','b', ...
                    'Title','Miscellaneous', ...
                    'Units','pixels', ...
                    'Position',[PanelsLeft VerticalPosition ...
                                PanelsWidth PanelHeight]);

NF_Checkbox = uicheckbox(MiscPanel,'Compute near fields',0,[10 60 130 21], ...
           @MiscCallback);
FF_Checkbox = uicheckbox(MiscPanel,'Compute far fields',0,[10 35 130 21], ...
           @MiscCallback);
EnergyCheckbox = uicheckbox(MiscPanel,'Compute energy every',0,[10 10 135 21], ...
           @MiscCallback);

NF_SettingsPush = uipushbutton(MiscPanel,'Settings...',[140 60 75 23], ...
                               @NF_SettingsCallback);
FF_SettingsPush = uipushbutton(MiscPanel,'Settings...',[140 35 75 23], ...
                               @FF_SettingsCallback);
EnergyEdit = uiedit(MiscPanel,'',[145 10 60 21],'Callback',@MiscCallback);
uitext(MiscPanel,'  step',[205 10 30 17]);


function MiscCallback(src,eventdata)
    
    switch src,
        case NF_Checkbox
            v.field_flag = get(src,'Value');
            if v.field_flag==1, set(NF_SettingsPush,'Enable','on');
            else                set(NF_SettingsPush,'Enable','off');
            end;
        case FF_Checkbox
            v.nfff_flag = get(src,'Value');
            if v.nfff_flag==1, set(FF_SettingsPush,'Enable','on');
            else               set(FF_SettingsPush,'Enable','off');
            end;
        case EnergyCheckbox
            v.energy_flag = get(src,'Value');
            if v.energy_flag==1,
                set([EnergyEdit,DropCheckbox],'Enable','on');
            else
                set([EnergyEdit,DropCheckbox,Drop_Edit],'Enable','off');
                set(DropCheckbox,'Value',0);
                v.drop_flag = false;
            end;
        case EnergyEdit
            value = str2num( get(src,'String') );
            if isempty(value) || ~isfinite(value) || imag(value)~=0 || value<=0,
                errordlg('Energy computation interval must be real and positive numeric value', ...
                    'Bad input','modal');
                pause(1);
                uicontrol(src);
                return;
            end;
            v.e_step = value;
    end;
    
end


function NF_SettingsCallback(src,eventdata)
    
    NF_Dialog = dialog('Visible','off', ...
        'Name','Near Field Settings', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 70] 260 130], ...
        'KeyPressFcn',{@RetEscPressed,@NF_Store,@CloseFigure});

    uitext(NF_Dialog,'Frequencies (comma separated):',[10 100 240 17]);
    str = num2str(v.freq_nf,' %.4g,'); str(end) = [];
    freqnf_Edit = uiedit(NF_Dialog,str,[10 75 220 21]);
    uitext(NF_Dialog,'  Hz',[230 75 20 17]);

    %% Field value normalizing
    FieldCheckbox = uicheckbox(NF_Dialog, ...
        'Near field normalizing to 1 W of input power', ...
        v.FieldNormalizing,[10 45 240 21],'');

    %% default behavior of Enter and Escape on controls
    set(freqnf_Edit,'KeyPressFcn',{@RetEscPressed,@NF_Store,@CloseFigure});

    %% buttons
    OKButton = Button(NF_Dialog,'OK',2,@NF_Store);
    CancelButton = Button(NF_Dialog,'Cancel',1,@CloseFigure);

    %% locking the uicontrols
    if v.locked_flag, LockDialog(NF_Dialog,OKButton,CancelButton); end;
    
    uicontrol(freqnf_Edit); %% give focus to the first field
    set(NF_Dialog,'Visible','on'); %% display figure

    function NF_Store(src,eventdata);
        
        %% writing the NF frequency
        value = str2num( get(freqnf_Edit,'String') );
        if isempty(value) || ~all(isfinite(value)) || ...
           any(imag(value)~=0) || any(value<=0),
            errordlg('near field parameters must be real numeric values', ...
                'Bad input','modal');
            pause(1);
            uicontrol(freqnf_Edit);
            return;
        end;
        v.freq_nf = value;
        
        %% writing field normalizing flag
        v.FieldNormalizing = get(FieldCheckbox,'Value');
        
        close(NF_Dialog);
        
    end
    
end


function FF_SettingsCallback(src,eventdata)
    FF_Dialog = dialog('Visible','off', ...
        'Name','Far Field Settings', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 70] 510 265], ...
        'KeyPressFcn',{@RetEscPressed,@FF_Store,@CloseFigure});
    
    uitext(FF_Dialog,'Frequencies (comma separated):',[10 235 240 17]);
    str = num2str(v.freq_ff,' %.4g,'); str(end) = [];
    freqff_Edit = uiedit(FF_Dialog,str,[10 210 250 21]);
    uitext(FF_Dialog,'  Hz',[260 210 20 17]);

    %% Integration surfaces
    FF_IntSurfPanel = uipanel('Parent',FF_Dialog, ...
                              'ForegroundColor','b', ...
                              'Title','Integration surfaces', ...
                              'Units','pixels', ...
                              'Position',[10 45 270 150]);
    
    uitext(FF_IntSurfPanel,'lower',[100 111 60 16],'HorizontalAlignment','center');
    uitext(FF_IntSurfPanel,'upper',[195 111 60 16],'HorizontalAlignment','center');

    uitext(FF_IntSurfPanel,'x-distance',[10 90 70 17]);
    uitext(FF_IntSurfPanel,'y-distance',[10 65 70 17]);
    uitext(FF_IntSurfPanel,'z-distance',[10 40 70 17]);
    uitext(FF_IntSurfPanel,'Units:',    [10 10 40 17]);
                          
    FF_xdist1 = uiedit(FF_IntSurfPanel,num2str(v.x_nfff(1),'%d'),[100 90 60 21]);
    FF_xdist2 = uiedit(FF_IntSurfPanel,num2str(v.x_nfff(2),'%d'),[195 90 60 21]);
    FF_ydist1 = uiedit(FF_IntSurfPanel,num2str(v.y_nfff(1),'%d'),[100 65 60 21]);
    FF_ydist2 = uiedit(FF_IntSurfPanel,num2str(v.y_nfff(2),'%d'),[195 65 60 21]);
    FF_zdist1 = uiedit(FF_IntSurfPanel,num2str(v.z_nfff(1),'%d'),[100 40 60 21]);
    FF_zdist2 = uiedit(FF_IntSurfPanel,num2str(v.z_nfff(2),'%d'),[195 40 60 21]);
    
    FF_face(1) = uicheckbox(FF_IntSurfPanel,[],v.nfff_faces(1),[ 80 90 20 21],'');
    FF_face(2) = uicheckbox(FF_IntSurfPanel,[],v.nfff_faces(2),[175 90 20 21],'');
    FF_face(3) = uicheckbox(FF_IntSurfPanel,[],v.nfff_faces(3),[ 80 65 20 21],'');
    FF_face(4) = uicheckbox(FF_IntSurfPanel,[],v.nfff_faces(4),[175 65 20 21],'');
    FF_face(5) = uicheckbox(FF_IntSurfPanel,[],v.nfff_faces(5),[ 80 40 20 21],'');
    FF_face(6) = uicheckbox(FF_IntSurfPanel,[],v.nfff_faces(6),[175 40 20 21],'');
    
    FF_IntSurfUnits = uipopupmenu(FF_IntSurfPanel,Units,1,[50 10 55 21]);
    
    %% Phase center
    FF_CenterPanel = uipanel('Parent',FF_Dialog, ...
                             'ForegroundColor','b', ...
                             'Title','Phase center', ...
                             'Units','pixels', ...
                             'Position',[290 125 215 100]);
    
    uitext(FF_CenterPanel,'x',[ 10 61 60 16],'HorizontalAlignment','center');
    uitext(FF_CenterPanel,'y',[ 75 61 60 16],'HorizontalAlignment','center');
    uitext(FF_CenterPanel,'z',[140 61 60 16],'HorizontalAlignment','center');

    uitext(FF_CenterPanel,'Units:',[10 10 40 17]);
                          
    FF_center1 = uiedit(FF_CenterPanel,num2str(v.phcenter(1),'%g'),[ 10 40 60 21]);
    FF_center2 = uiedit(FF_CenterPanel,num2str(v.phcenter(2),'%g'),[ 75 40 60 21]);
    FF_center3 = uiedit(FF_CenterPanel,num2str(v.phcenter(3),'%g'),[140 40 60 21]);
    
    FF_CenterUnits = uipopupmenu(FF_CenterPanel,Units,1,[50 10 55 21]);
    
    %% Radiation pattern resolution
    FF_PatternPanel = uipanel('Parent',FF_Dialog, ...
                              'ForegroundColor','b', ...
                              'Title','Radiation pattern resolution', ...
                              'Units','pixels', ...
                              'Position',[290 45 215 70]);
    
    uitext(FF_PatternPanel,'theta',[ 10 31 60 16],'HorizontalAlignment','center');
    uitext(FF_PatternPanel,'phi',  [100 31 60 16],'HorizontalAlignment','center');

    ThetaResEdit = uiedit(FF_PatternPanel,num2str(v.resolution(1),'%d'),[ 10 10 60 21]);
    uitext(FF_PatternPanel,char(215),[70 10 30 17],'HorizontalAlignment','center');
    PhiResEdit   = uiedit(FF_PatternPanel,num2str(v.resolution(2),'%d'),[100 10 60 21]);

    FF_Edit = [freqff_Edit ...
               FF_xdist1 FF_xdist2 FF_ydist1 FF_ydist2 FF_zdist1 FF_zdist2 ...
               FF_center1 FF_center2 FF_center3 ...
               ThetaResEdit PhiResEdit];

    %% default behavior of Enter and Escape on controls
    set(FF_Edit,'KeyPressFcn',{@RetEscPressed,@FF_Store,@CloseFigure});

    %% buttons
    OKButton = Button(FF_Dialog,'OK',2,@FF_Store);
    CancelButton = Button(FF_Dialog,'Cancel',1,@CloseFigure);

    %% locking the uicontrols
    if v.locked_flag, LockDialog(FF_Dialog,OKButton,CancelButton); end;
    
    uicontrol(freqff_Edit); %% give focus to the first field
    set(FF_Dialog,'Visible','on'); %% display figure

    function FF_Store(src,eventdata);
        
        %% checking and writing the FF frequency
        value = str2num( get(freqff_Edit,'String') );
        if isempty(value) || ~all(isfinite(value)) || ...
           any(imag(value)~=0) || any(value<=0),
            errordlg('far field parameters must be real numeric values', ...
                'Bad input','modal');
            pause(1);
            uicontrol(freqff_Edit);
            return;
        end;
        v.freq_ff = value;

        %% checking and writing the rest
        for h = FF_Edit(2:end),
            value = str2num( get(h,'String') );
            if isempty(value) || ~isfinite(value) || imag(value)~=0,
                errordlg('far field parameters must be real numeric values', ...
                    'Bad input','modal');
                pause(1);
                uicontrol(h);
                return;
            end;
            switch h,
                case FF_xdist1, v.x_nfff(1) = value;
                case FF_xdist2, v.x_nfff(2) = value;
                case FF_ydist1, v.y_nfff(1) = value;
                case FF_ydist2, v.y_nfff(2) = value;
                case FF_zdist1, v.z_nfff(1) = value;
                case FF_zdist2, v.z_nfff(2) = value;
                case FF_center1, v.phcenter(1) = value;
                case FF_center2, v.phcenter(2) = value;
                case FF_center3, v.phcenter(3) = value;
                case ThetaResEdit, v.resolution(1) = value;
                case PhiResEdit,   v.resolution(2) = value;
            end;
        end;
        for i = 1:6,
            v.nfff_faces(i) = get(FF_face(i),'Value');
        end;
        close(FF_Dialog);
    end
end


%% finding all controls for deactivation
UIActive = findobj(Panel,'Type','uicontrol');
set(UIActive,'Enable','off');

%% activation
set(MainFig,'Visible','on');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%% MENU CALLBACKS %%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Open dialog
function Open(src,eventdata);

    if isequal(src,MenuNew), %% new file is being opened

        [FileName,PathName] = uiputfile( ...
            {'*.aaf','AAU3 FDTD Files (*.aaf)'; ...
             '*','All Files (*.*)'}, ...
            'New File...',[WorkDir,filesep]);

        if isequal(FileName,0) || isequal(PathName,0), return; end;
        WorkDir = PathName;
        try
            save(fullfile(PathName,FileName),'-struct','DefaultVars',SavingVersion);
        catch
            errordlg('Create operation failed!','Create Error');
            return;
        end;
        v = DefaultVars;
        
    elseif isequal(src,MenuOpen), %% existing file being opened

        [FileName,PathName] = uigetfile( ...
            {'*.aaf','AAU3 FDTD Files (*.aaf)'; ...
             '*.bdf','Binary Definition Files (*.bdf)'; ...
             '*','All Files (*.*)'}, ...
            'Open',[WorkDir,filesep]);

        if isequal(FileName,0) || isequal(PathName,0), return; end;
        WorkDir = PathName;
        try
            v = load('-mat',fullfile(PathName,FileName));
        catch
            errordlg(['Error when opening file ',FileName],'Open Error');
            return;
        end;
        
    elseif isequal(src,MenuDeployMulti) || ...
           isequal(src,MenuCollectMulti) || ...
           isequal(src,MenuCompareResults), %% other cases
        
        FileName = eventdata;
        PathName = WorkDir;
        try
            v = load('-mat',fullfile(PathName,FileName));
        catch
            errordlg(['Error when opening file ',FileName],'Open Error');
            return;
        end;
        
    else
        error('Unspecified OPEN callback source');
        
    end;

    ActualFileName = FileName;
    ActualPathName = PathName;

    %% Compatibility issues (again!)
    DefaultFieldNames = fieldnames(DefaultVars);
    LoadedFieldNames  = fieldnames(v);
    MissingFields = setdiff(DefaultFieldNames,LoadedFieldNames);
    if length(MissingFields) > 0,
        for FieldName = MissingFields.',
            v = setfield(v,FieldName{1},getfield(DefaultVars,FieldName{1}));
        end;
    end;
    v = orderfields(v,DefaultVars);
    
    for GeomPointer = 1:length(v.Geometry),
        
        switch v.Geometry{GeomPointer}.type
            
            case 'source'
                if ~isfield(v.Geometry{GeomPointer},'resistance'),
                    v.Geometry{GeomPointer}.resistance = DefaultSourceResistance;
                end;
                
            case 'brick'
                if ~isfield(v.Geometry{GeomPointer},'massdensity'),
                    v.Geometry{GeomPointer}.massdensity = DefaultMassDensity;
                end;
                
        end;

    end;
    
    %% enabling uicontrols and menus
    set(UIActive,'Enable','on');
    set([GeomPushEdit,GeomPushDelete],'Enable','off');
    set(GeomListBox,'UIContextMenu',GeomContextMenu);
    set(get(MenuFile,      'Children'),'Enable','on');
    set(get(MenuInsert,    'Children'),'Enable','on');
    set(get(MenuSimulation,'Children'),'Enable','on');
    set(get(MenuResults,   'Children'),'Enable','off');
    
    set(MenuEmagField, 'Enable','off');
    set(MenuHmagField, 'Enable','off');
    set(MenuJmagField, 'Enable','off');
    set(MenuSARDisp,   'Enable','off');
    set(MenuSAR1gDisp, 'Enable','off');
    set(MenuSAR10gDisp,'Enable','off');
    
    set(MenuSARAver1g, 'Enable','off');
    set(MenuSARAver10g,'Enable','off');

    %% only some menus should be enabled depending on the engine
    if pref.ActiveEngine == 1, % Matlab engine
        set([MenuDeployMulti,MenuCollectMulti, ...
             MenuDeploy,MenuCollect,MenuSubmit],'Enable','off');
    else % other engines
        switch Engine.Processing
            case 'Run immediately', set(MenuSubmit,'Enable','off');
            case 'Submit to queue', set(MenuRun   ,'Enable','off');
        end;
    end;
    
    %% filling the editboxes and setting controls
    set(GeomListBox,'ListboxTop',1,'String',v.GeometryString,'Value',[]);
    set(Dx_Edit,'String',num2str(v.Dx,'%.4g'));
    set(Df_Edit,'String',num2str(v.Df,'%.4g'));
    set(Nt_Edit,'String',num2str(v.Nt,'%d'));
    set(Drop_Edit ,'String',num2str(v.Drop_dB ,'%.4g'));
    set(Drop2_Edit,'String',num2str(v.Drop2_dB,'%.4g'));
    set(DropCheckbox ,'Value',v.drop_flag);
    set(Drop2Checkbox,'Value',v.drop2_flag);
    if v.drop_flag, set(Drop_Edit,'Enable','on');
    else            set(Drop_Edit,'Enable','off');
    end;
    if v.drop2_flag, set(Drop2_Edit,'Enable','on');
    else             set(Drop2_Edit,'Enable','off');
    end;
    set(Domain_xlim1,'String',num2str(v.limits(1,1)));
    set(Domain_xlim2,'String',num2str(v.limits(2,1)));
    set(Domain_ylim1,'String',num2str(v.limits(1,2)));
    set(Domain_ylim2,'String',num2str(v.limits(2,2)));
    set(Domain_zlim1,'String',num2str(v.limits(1,3)));
    set(Domain_zlim2,'String',num2str(v.limits(2,3)));
    set(PML_Checkbox,'Value',v.pml);
    if v.pml==1, set(PML_SettingsPush,'Enable','on');
    else         set(PML_SettingsPush,'Enable','off');
    end;
    set(NF_Checkbox,'Value',v.field_flag);
    set(FF_Checkbox,'Value',v.nfff_flag);
    set(EnergyCheckbox,'Value',v.energy_flag);
    set(EnergyEdit,'String',num2str(v.e_step,'%d'));
    if v.field_flag==1, set(NF_SettingsPush,'Enable','on');
    else                set(NF_SettingsPush,'Enable','off');
    end;
    if v.nfff_flag==1, set(FF_SettingsPush,'Enable','on');
    else               set(FF_SettingsPush,'Enable','off');
    end;
    if v.energy_flag==1, set([EnergyEdit,DropCheckbox],'Enable','on');
    else                 set([EnergyEdit,DropCheckbox],'Enable','off');
    end;
    set(MainFig,'Name',[FigTitle,' - ',ActualFileName]);
    
    %% creating the geometry figure if it does not exist
    if isempty(GeometryFig) || ~ishandle(GeometryFig),
        GeometryFig = figure( ...
            'Visible','off', ...
            'Color','w', ...
            'NumberTitle','off', ...
            'Units','pixels', ...
            'Position',[30+FigWidth ScreenSize(4)-FigHeight-60 700 FigHeight-30], ...
            'Renderer',pref.Renderer, ...
            'Resize','on', ...
            'Toolbar','figure', ...
            'DockControls','on', ...
            'DeleteFcn',@Close );
    end;

    %% putting the filename in titlebar
    set(GeometryFig,'Name',['Geometry - ',ActualFileName], ...
                    'Visible','on');

    %% when there's old axes, delete it
    if ishandle(GeomAxes), delete(GeomAxes); end;

    %% displaying axes
    GeomAxes = axes( ...
        'Visible','off', ...
        'Parent',GeometryFig, ...
        'ActivePositionProperty','outerposition', ...
        'Box','on', ...
        'DataAspectRatio',[1 1 1], ...
        'Projection','orthographic' );

    axes(GeomAxes);
    axis(GeomAxes,v.limits(:));
    view(142.5,30);
    light;
    axis(GeomAxes,'vis3d');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    for GeomPointer = 1:length(v.Geometry),
        Type = v.Geometry{GeomPointer}.type;
        switch Type
            case 'wire',         DisplayWire(GeomPointer,'k');
            case 'triangle',     DisplayPlate(GeomPointer);
            case 'rectangle',    DisplayPlate(GeomPointer);
            case 'brick',        DisplayBrick(GeomPointer);
            case 'capacitor',    DisplayWire(GeomPointer,'g');
            case 'resistor',     DisplayWire(GeomPointer,'b');
            case 'source',       DisplaySource(GeomPointer);
            case 'field source', DisplayFieldSource(GeomPointer);
            case 'probe',        DisplayProbe(GeomPointer);
            case 'extobj',       DisplayExtObj(GeomPointer);
        end;
    end;

    if ~isempty(v.V), %% results are present

        set(MenuSourceVoltage,'Enable','on');
        set(MenuInputVoltage, 'Enable','on');
        set(MenuInputCurrent, 'Enable','on');
        set(MenuProbeResults, 'Enable','on');
        set(MenuImpedance,    'Enable','on');
        set(MenuQandBW,       'Enable','on');
        set(MenuCPUTime,      'Enable','on');
        if v.nfff_flag, set(MenuPatt,'Enable','on'); end;
        if v.nfff_flag, set(MenuEfficiency,'Enable','on'); end;
        if ~isempty(v.energy), set(MenuEnergy,'Enable','on'); end;
        if v.field_flag, set([MenuNearFields,MenuCurrents,MenuPostprocessing],'Enable','on'); end;
        if ~isempty(v.SAR), set([MenuSARAver1g,MenuSARAver10g,MenuSARDisp],'Enable','on'); end;
        if ~isempty(v.SAR_Aver_1g),  set(MenuSAR1gDisp ,'Enable','on'); end;
        if ~isempty(v.SAR_Aver_10g), set(MenuSAR10gDisp,'Enable','on'); end;
        if ~isempty(v.E_field),      set(MenuEmagField ,'Enable','on'); end;
        if ~isempty(v.H_field),      set(MenuHmagField ,'Enable','on'); end;
        if ~isempty(v.J_field),      set(MenuJmagField ,'Enable','on'); end;
        
    end;
        
    %% Locking the controls
    if ishandle(LockImageAxes), delete(LockImageAxes); end;
    if v.locked_flag, LockMainFig; end;
    
    %% curtain up!
    set(GeomAxes,'Visible','on');
    
end;
    
    
%% Close file and window
function Close(src,eventdata)
        
    v = [];
    delete(GeometryFig);
    
    if ishandle(LockImageAxes), delete(LockImageAxes); end;

    %% setting uicontrols and menus inactive
    set(UIActive,'Enable','off');
    set(GeomListBox,'String','','UIContextMenu',[]);

    set(get(MenuFile,      'Children'),'Enable','off');
    set(get(MenuInsert,    'Children'),'Enable','off');
    set(get(MenuSimulation,'Children'),'Enable','off');
    set(get(MenuResults,   'Children'),'Enable','off');

    set(MenuNew,           'Enable','on');
    set(MenuOpen,          'Enable','on');
    set(MenuDeployMulti,   'Enable','on');
    set(MenuCollectMulti,  'Enable','on');
    set(MenuCompareResults,'Enable','on');
    set(MenuBDFConvert,    'Enable','on');
    set(MenuPreferences,   'Enable','on');
    set(MenuExit,          'Enable','on');

    set(MenuEmagField, 'Enable','off');
    set(MenuHmagField, 'Enable','off');
    set(MenuJmagField, 'Enable','off');
    set(MenuSARDisp,   'Enable','off');
    set(MenuSAR1gDisp, 'Enable','off');
    set(MenuSAR10gDisp,'Enable','off');
    
    set(MenuSARAver1g, 'Enable','off');
    set(MenuSARAver10g,'Enable','off');

    set(findobj(Panel,'Type','uicontrol','Style','edit'),       'String','');
    set(findobj(Panel,'Type','uicontrol','Style','radiobutton'),'Value',0);
    set(findobj(Panel,'Type','uicontrol','Style','checkbox'),   'Value',0);
    
    set(MainFig,'Name',FigTitle);
    
end
   

%% Save
function Save(src,eventdata)
        
    try
        save(fullfile(ActualPathName,ActualFileName),'-struct','v',SavingVersion);
    catch
        errordlg('Save operation failed!','Save Error');
        return;
    end;
        
end


%% Save Ass ;-) dialog
function SaveAs(src,eventdata)
        
    [FileName,PathName] = uiputfile( ...
        {'*.aaf','AAU3 FDTD Files (*.aaf)'; ...
         '*','All Files (*.*)'}, ...
        'Save As',fullfile(ActualPathName,ActualFileName));

    if isequal(FileName,0) || isequal(PathName,0), return; end;
    WorkDir = PathName;
    try
        save(fullfile(PathName,FileName),'-struct','v',SavingVersion);
    catch
        errordlg('Save operation failed!','Save Error');
        return;
    end;
    
    ActualFileName = FileName;
    ActualPathName = PathName;
    set(GeometryFig,'Name',[  'Geometry - ',ActualFileName]);
    set(MainFig,    'Name',[ FigTitle,' - ',ActualFileName]);
        
end


%% Export dialog
function Export(src,eventdata)

    %% indices into a list of variables (provided that the order is kept)
    VariableNames = fieldnames(v);
    Index_of_V = strmatch('V',VariableNames,'exact');
    SimParVariables = 1:Index_of_V-1;
    ResultVariables = Index_of_V:length(VariableNames);

    %% main dialog window
    ExportDialog = dialog( ...
        'Visible','off', ...
        'Name','Export...', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[150 135] 250 250], ...
        'KeyPressFcn',{@RetEscPressed,@ExportProceed,@CloseFigure} );
    
    %% radiobuttons group
    ExportOptionsGroup = uiradiogroup(ExportDialog,[10 165 230 75], ...
        @ExportOptionsCallback);
    ExportOption1 = uiradiobutton(ExportOptionsGroup,1,[5 55 220 15], ...
        'String','Simulation Parameters - Unlocked');
    ExportOption2 = uiradiobutton(ExportOptionsGroup,0,[5 30 220 15], ...
        'String','Results');
    ExportOption3 = uiradiobutton(ExportOptionsGroup,0,[5  5 220 15], ...
        'String','Selected Variables:');
    
    function ExportOptionsCallback(src,eventdata)
        if eventdata.NewValue==ExportOption3,
             set(ExportListbox,'Enable','on');
        else set(ExportListbox,'Enable','off');
        end;
    end
    
    %% list of variables
    ExportListbox = uicontrol(ExportDialog, ...
        'Style','listbox', ...
        'BackgroundColor','w', ...
        'String',VariableNames, ...
        'Min',0,'Max',2, ...
        'Units','pixels', ...
        'Position',[40 55 200 100], ...
        'Enable','off');

    %% default behavior of Enter and Escape on controls
    set([ExportOption1,ExportOption2,ExportOption3,ExportListbox], ...
        'KeyPressFcn',{@RetEscPressed,@ExportProceed,@CloseFigure});

    %% buttons
    OKButton = Button(ExportDialog,'OK',2,@ExportProceed);
    CancelButton = Button(ExportDialog,'Cancel',1,@CloseFigure);

    uicontrol(ExportOption1); %% give focus to the first field
    set(ExportDialog,'Visible','on'); %% display figure
    
    function ExportProceed(src,eventdata)
        
        %% first store the important things from the dialog...
        OptionSelected = get(ExportOptionsGroup,'SelectedObject');
        if OptionSelected==ExportOption3,
            VariablesSelected = get(ExportListbox,'Value');
        end
        
        %% ...before we close it
        close(ExportDialog);

        %% filename retrieval
        [null,FileNameBody,null] = fileparts(ActualFileName);
        switch OptionSelected
            case ExportOption1
                FileNameSuggest = [FileNameBody '_unlocked.aaf'];
                [FileName,PathName] = uiputfile( ...
                    {'*.aaf','AAU3 FDTD Files (*.aaf)'; ...
                    '*.mat','Matlab MAT-files (*.mat)'; ...
                    '*','All Files (*.*)'}, ...
                    'Export to File',fullfile(ActualPathName,FileNameSuggest));
            case ExportOption2
                FileNameSuggest = [FileNameBody '_results.mat'];
                [FileName,PathName] = uiputfile( ...
                    {'*.mat','Matlab MAT-files (*.mat)'; ...
                    '*','All Files (*.*)'}, ...
                    'Export to File',fullfile(ActualPathName,FileNameSuggest));
            case ExportOption3
                FileNameSuggest = [FileNameBody '_variables.mat'];
                [FileName,PathName] = uiputfile( ...
                    {'*.mat','Matlab MAT-files (*.mat)'; ...
                    '*','All Files (*.*)'}, ...
                    'Export to File',fullfile(ActualPathName,FileNameSuggest));
        end;

        if isequal(FileName,0) || isequal(PathName,0), return; end;
        WorkDir = PathName;

        switch OptionSelected
            case ExportOption1, VariableNumbers = SimParVariables;
            case ExportOption2, VariableNumbers = ResultVariables;
            case ExportOption3, VariableNumbers = VariablesSelected;
        end

        %% copying only the desired variables
        for i = VariableNumbers,
            v_temp.(VariableNames{i}) = v.(VariableNames{i});
        end;
        
        %% if option 1, unlocking and adding the results variables as empty
        if OptionSelected==ExportOption1,
            v_temp.locked_flag = false;
            for i = ResultVariables,
                v_temp.(VariableNames{i}) = [];
            end;
        end;
        
        try
            save(fullfile(PathName,FileName),'-struct','v_temp',SavingVersion);
        catch
            errordlg('Save operation failed!','Save Error');
            return;
        end;
        
    end
    
end


%% Deploy multiple files for processing by Fortran kernel
function DeployMulti(src,eventdata)
    
    %% which files to take
    [FileNames,PathName] = uigetfile( ...
        {'*.aaf','AAU3 FDTD Files'; ...
         '*','All Files (*.*)'}, ...
        'Deploy Multiple Files for Fortran Processing',[WorkDir,filesep], ...
        'MultiSelect','on');

    if isequal(FileNames,0) || isequal(PathName,0), return; end;
    WorkDir = PathName;
    
    if ~iscell(FileNames), FileNames = cellstr(FileNames); end;
    
    SuccessCount = 0;
    
    for n = 1:length(FileNames),
        
        Open(src,FileNames{n});
        
        if v.locked_flag,
            h = msgbox({['File ',FileNames{n},' already contains results.'], ...
                        'Not deploying.'},'Error','error','modal');
            uiwait(h);
            continue;
        end;
        
        SuccessFlag = false;
        
        SimulationDeploy;
        
        if SuccessFlag, SuccessCount = SuccessCount + 1; end;
        
    end;
    
    Close;

    msgbox([num2str(SuccessCount),' file(s) successfully deployed.'],'Info','modal');
    
end


%% Collect multiple results from Fortran
function CollectMulti(src,eventdata)
    
    %% which files to take
    [FileNames,PathName] = uigetfile( ...
        {'*.aaf','AAU3 FDTD Files'; ...
         '*','All Files (*.*)'}, ...
        'Collect Multiple Fortran Results',[WorkDir,filesep],'MultiSelect','on');

    if isequal(FileNames,0) || isequal(PathName,0), return; end;
    
    if ~iscell(FileNames), FileNames = {FileNames}; end;
    WorkDir = PathName;
    
    %% where to put the results
    PathName = uigetdir(WorkDir,'Save Files to Folder...');
    if isequal(PathName,0), return; end;
    SaveDir = PathName;
    
    %% here should be Warning if the WorkDir and SaveDir are the same
    
    for n = 1:length(FileNames),
        
        Open(src,FileNames{n});
        
        SimulationCollect;
        
        try
            save(fullfile(SaveDir,ActualFileName),'-struct','v',SavingVersion);
        catch
            errordlg('Save operation failed!','Save Error');
            return;
        end;

    end;
    
    Close;
    
end



function CompareResults(src,eventdata)

    %% which variable to plot
    func_args = get(src,'UserData');
    if iscell(func_args), func_handle = func_args{1};
    else                  func_handle = func_args;
    end;
    
    %% for patterns (dirty!)
    h = get(get(src,'Parent'),'Parent');
    if isequal(get(h,'Type'),'uimenu') && ...
       isequal(get(h,'Label'),'Radiation Pattern'),
        pol_label = get(get(src,'Parent'),'Label');
    else
        pol_label = [];
    end;

    %% which files to take
    [FileNames,PathName] = uigetfile( ...
        {'*.aaf','AAU3 FDTD Files'; ...
         '*','All Files (*.*)'}, ...
        'Compare Multiple Results',[WorkDir,filesep],'MultiSelect','on');

    if isequal(FileNames,0) || isequal(PathName,0), return; end;
    
    if ~iscell(FileNames), FileNames = {FileNames}; end;
    WorkDir = PathName;

    for n = 1:length(FileNames),
        
        Open(MenuCompareResults,FileNames{n});
        
        if ~v.locked_flag,
            h = msgbox({['File ',FileNames{n},' does not contain results.'], ...
                        'Aborting.'},'Error','error','modal');
            uiwait(h);
            continue;
        end;

        %% adding the next plot to the axes
        if n > 1, eventdata = 'add'; else eventdata = []; end;
        
        %% running the selected display function
        if iscell(func_args),
            switch length(func_args)
                case 2, func_handle(MenuCompareResults,eventdata,func_args{2});
                case 3, func_handle(MenuCompareResults,eventdata,func_args{2},func_args{3});
            end;
        else
            func_handle(MenuCompareResults,eventdata);
        end;
        
    end;
    
    box(plot_axes,'on');
    legend(plot_axes,FileNames,'Interpreter','none');
    
    Close;

end



%%%%%%%%%%%%%%%%%%%%%%%%
%% PREFERENCES DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%%%
function Preferences(src,eventdata)
    
    PrefDialog = dialog('Visible','off', ...
        'Name','Preferences', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[30+FigWidth ScreenSize(4)-420-40 267 420], ...
        'KeyPressFcn',{@RetEscPressed,@PrefSave,@CloseFigure});

    %% Engines panel with drop-down menu and control buttons
    EnginePanel = uipanel('Parent',PrefDialog, ...
                          'ForegroundColor','b', ...
                          'Title','Simulation engine', ...
                          'Units','pixels', ...
                          'Position',[7 336 253 80]);

    %% making the cell array of strings for the drop-down menu
    for i = 1:length(pref.Engines),
        EngineString{i} = pref.Engines(i).Name;
    end;
                      
    %% drop-down menu
    EngineDropMenu = uipopupmenu(EnginePanel,EngineString, ...
        pref.ActiveEngine,[8 40 235 21],'UserData',pref.Engines);
    
    %% pushbuttons                                        
    EngineAdd    = uipushbutton(EnginePanel,'Add...',[  8 10 75 23], ...
                               @EngineCallback,'Enable','on');
    EngineModify = uipushbutton(EnginePanel,'Modify...',  [ 88 10 75 23], ...
                               @EngineCallback,'Enable','on');
    EngineDelete = uipushbutton(EnginePanel,'Delete...',[168 10 75 23], ...
                               @EngineCallback,'Enable','on');

    %% callback of the listbox and associated controls
    function EngineCallback(src,eventdata)

        Engines       = get(EngineDropMenu,'UserData');
        EnginePointer = get(EngineDropMenu,'Value'); %% which item is active

        switch src
            
        %% ADD : Adding new engine
        case EngineAdd,
            
            %% open window for editing properties of engines
            EnginesDefine('add',EngineDropMenu);

        %% MODIFY : Modifying existing engine
        case EngineModify,

            %% Matlab engine cannot be modified
            if strcmp(Engines(EnginePointer).Name,'Matlab')
                errordlg('Matlab engine is built-in and cannot be modified!','Modify Error');
                return;
            end;

            %% open window for editing properties of engines
            EnginesDefine('modify',EngineDropMenu);

        %% DELETE : Deleting engine list items
        case EngineDelete,

            %% Matlab engine cannot be deleted
            if strcmp(Engines(EnginePointer).Name,'Matlab')
                errordlg('Matlab engine is built-in and cannot be deleted!','Delete Error');
                return;
            end;

            Button = questdlg('Delete selected engine?','Delete Engine','No');

            if strcmp(Button,'Yes'),

                %% delete the engine
                Engines(EnginePointer) = [];

                %% update the string in the drop-down menu
                EngineString = get(EngineDropMenu,'String');
                EngineString(EnginePointer) = [];
                set(EngineDropMenu,'String',EngineString, ...
                                   'Value',EnginePointer-1, ...
                                   'UserData',Engines);
                
            end;
            
        end;

    end

    %% Geometry display panel
    GeomDispPanel = uipanel('Parent',PrefDialog, ...
                            'ForegroundColor','b', ...
                            'Title','Geometry display', ...
                            'Units','pixels', ...
                            'Position',[7 254 253 80]);
    
    %% Renderer
    uitext(GeomDispPanel,'Renderer:',[10 40 60 17]);
    RenderPopup = uipopupmenu(GeomDispPanel,Renderers, ...
        strmatch(pref.Renderer,Renderers),[70 40 70 21]);

    %% Grid
    GridCheckbox = uicheckbox(GeomDispPanel,'Display grid on objects', ...
        pref.DisplayGrid,[10 10 230 21],'');

    %% User interface panel
    UIPanel = uipanel('Parent',PrefDialog, ...
                      'ForegroundColor','b', ...
                      'Title','User Interface', ...
                      'Units','pixels', ...
                      'Position',[7 172 253 80]);

    %% Restart warning
    uitext(UIPanel,'(These settings need restart to take effect)', ...
        [10 40 230 17]);

    %% Color and Font pushbuttons                                        
    ColorSetup = uipushbutton(UIPanel,'Color...',[  8 10 75 23], ...
                              @UISetupCallback,'Enable','on', ...
                                               'UserData',pref.Color);
    FontSetup  = uipushbutton(UIPanel,'Font...',  [ 88 10 75 23], ...
                              @UISetupCallback,'Enable','on', ...
                                               'UserData',pref.Font);
                          
    function UISetupCallback(src,eventdata)
        
        switch src
            
            case ColorSetup
                Color = uisetcolor(get(src,'UserData'),'GUI Background Color');
                if ~isequal(Color,0), set(src,'UserData',Color); end;
                
            case FontSetup
                Font = uisetfont(get(src,'UserData'),'GUI Font Properties');
                if ~isequal(Font,0), set(src,'UserData',Font); end;
                
        end;
        
    end

    %% Miscellaneous panel
    MiscPrefPanel = uipanel('Parent',PrefDialog, ...
                            'ForegroundColor','b', ...
                            'Title','Miscellaneous', ...
                            'Units','pixels', ...
                            'Position',[7 45 253 125]);
    
    %% Check brick limits
    BrickCheckbox = uicheckbox(MiscPrefPanel,'Check brick limits', ...
        pref.CheckBrickLimits,[10 85 230 21],'');
    
    %% New saving
    SavingCheckbox = uicheckbox(MiscPrefPanel,'Saving in v7.3 format', ...
        pref.NewSaving,[10 60 230 21],'');
    
    %% New .in and .out files
    InOutCheckbox = uicheckbox(MiscPrefPanel,'New .in and .out files', ...
        pref.NewInOutFiles,[10 35 230 21],'');
    
    %% Cache optimization
    CacheCheckbox = uicheckbox(MiscPrefPanel,'Cache optimization for Opteron', ...
        pref.CacheOptim,[10 10 230 21],'');
    
    %% buttons
    SaveButton   = Button(PrefDialog,'Save'  ,2,@PrefSave);
    CancelButton = Button(PrefDialog,'Cancel',1,@CloseFigure);

    set(PrefDialog,'Visible','on'); %% display figure

    function PrefSave(src,eventdata);

        %% Engines
        pref.Engines      = get(EngineDropMenu,'UserData');
        pref.ActiveEngine = get(EngineDropMenu,'Value');
        
        %% shorter name for the active engine
        Engine = pref.Engines(pref.ActiveEngine);

        %% only some menus should be enabled depending on the engine
        if ~isempty(v) && ~v.locked_flag, % only if file open and not locked
            if pref.ActiveEngine == 1, % Matlab engine
                set([MenuDeployMulti;MenuCollectMulti; ...
                     get(MenuSimulation,'Children')],'Enable','off');
                set(MenuRun,'Enable','on');
            else % other engines
                set([MenuDeployMulti;MenuCollectMulti; ...
                     get(MenuSimulation,'Children')],'Enable','on');
                switch Engine.Processing
                    case 'Run immediately', set(MenuSubmit,'Enable','off');
                    case 'Submit to queue', set(MenuRun   ,'Enable','off');
                end;
            end;
        end;

        %% Renderer
        pref.Renderer = Renderers{get(RenderPopup,'Value')};
        if ishandle(GeometryFig),
            set(GeometryFig,'Renderer',pref.Renderer);
        end;
        
        %% Grid display
        if pref.DisplayGrid ~= get(GridCheckbox,'Value'), % save redraw
            pref.DisplayGrid = get(GridCheckbox,'Value');
            if pref.DisplayGrid, PatchLineStyle = '-';
            else                 PatchLineStyle = 'none';
            end;
            if ishandle(GeomAxes),
                set(findobj(GeomAxes,'Type','patch'),'LineStyle',PatchLineStyle);
            end;
        end;

        %% Font
        pref.Font = get(FontSetup,'UserData');
        
        %% Color
        pref.Color = get(ColorSetup,'UserData');

        %% Check brick limits
        pref.CheckBrickLimits = get(BrickCheckbox,'Value');
        
        %% New saving
        pref.NewSaving = get(SavingCheckbox,'Value');
        if pref.NewSaving, SavingVersion = '-v7.3';
        else               SavingVersion = '-v7';
        end;
        
        %% New .in and .out files
        pref.NewInOutFiles = get(InOutCheckbox,'Value');
        
        %% Cache optimization
        pref.CacheOptim = get(CacheCheckbox,'Value');

        %% Save the preferences and close the dialog window
        save(fullfile(AFCDir,'preferences.mat'),'-struct','pref');
        close(PrefDialog);
        
    end

end



%%%%%%%%%%%%%%%%%%%%
%% ENGINES DIALOG %%
%%%%%%%%%%%%%%%%%%%%
function EnginesDefine(mode,EngineDropMenu)
    
    %% Active engine in the drop-down menu
    Engines       = get(EngineDropMenu,'UserData');
    EnginePointer = get(EngineDropMenu,'Value');

    %% Open window with engine parameters
    EngineDialog = dialog('Visible','off', ...
        'Name','Engine Definition', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[50+FigWidth+267 ScreenSize(4)-405-40 448 405], ...
        'KeyPressFcn',{@RetEscPressed,@EngineSave,@CloseFigure});

    %% Engine name editbox
    uitext(EngineDialog,'Engine name:',[10 376 75 17]);
    EngineNameEdit = uiedit(EngineDialog,'',[85 376 353 21], ...
        'HorizontalAlignment','left','Tag','Name');
    
    %% Location panel with radiobuttons
    LocationPanel = uipanel('Parent',EngineDialog, ...
                            'ForegroundColor','b', ...
                            'Title','Location', ...
                            'Units','pixels', ...
                            'Position',[7 299 140 70]);
    LocationGroup = uiradiogroup(LocationPanel,[10 10 120 45], ...
        @ChangeCallback,'Tag','Location');
    LocationOptions(1) = uiradiobutton(LocationGroup,0,[0 25 120 20], ...
        'String','Local');
    LocationOptions(2) = uiradiobutton(LocationGroup,0,[0 0 120 20], ...
        'String','Remote');

    %% Job processing panel with radiobuttons
    ProcessingPanel = uipanel('Parent',EngineDialog, ...
                              'ForegroundColor','b', ...
                              'Title','Job processing', ...
                              'Units','pixels', ...
                              'Position',[154 299 140 70]);
    ProcessingGroup = uiradiogroup(ProcessingPanel,[10 10 120 45], ...
        @ChangeCallback,'Tag','Processing');
    ProcessingOptions(1) = uiradiobutton(ProcessingGroup,0,[0 25 120 20], ...
        'String','Run immediately');
    ProcessingOptions(2) = uiradiobutton(ProcessingGroup,0,[0 0 120 20], ...
        'String','Submit to queue');

    %% Script panel with radiobuttons
    ScriptPanel = uipanel('Parent',EngineDialog, ...
                          'ForegroundColor','b', ...
                          'Title','Run command/script', ...
                          'Units','pixels', ...
                          'Position',[301 299 140 70]);
    ScriptGroup = uiradiogroup(ScriptPanel,[10 10 120 45], ...
        @ChangeCallback,'Tag','Script');
    ScriptOptions(1) = uiradiobutton(ScriptGroup,0,[0 25 120 20], ...
        'String','Ready to use');
    ScriptOptions(2) = uiradiobutton(ScriptGroup,0,[0 0 120 20], ...
        'String','Edit from template');

    %% Host panel
    HostPanel = uipanel('Parent',EngineDialog, ...
                        'ForegroundColor','b', ...
                        'Title','Host', ...
                        'Units','pixels', ...
                        'Position',[7 222 434 75]);

    %% Host name editbox
    uitext(HostPanel,'Host name:',[10 35 65 17]);
    HostNameEdit = uiedit(HostPanel,'',[75 35 345 21], ...
        'HorizontalAlignment','left','Tag','HostName');

    %% Username editbox
    uitext(HostPanel,'Username:',[10 10 65 17]);
    UsernameEdit = uiedit(HostPanel,'',[75 10 125 21], ...
        'HorizontalAlignment','left','Tag','Username');

    %% Password editbox
    uitext(HostPanel,'Password:',[230 10 65 17]);
    PasswordEdit = uiedit(HostPanel,'',[295 10 125 21], ...
        'HorizontalAlignment','left','Tag','Password');

    %% Paths panel
    PathsPanel = uipanel('Parent',EngineDialog, ...
                         'ForegroundColor','b', ...
                         'Title','Paths', ...
                         'Units','pixels', ...
                         'Position',[7 45 434 175]);

    %% PuTTy path editbox
    uitext(PathsPanel,'PuTTy path:',[10 135 110 17]);
    PuTTyPathEdit = uiedit(PathsPanel,'',[120 135 275 21], ...
        'HorizontalAlignment','left','Tag','PuTTyPath');
    
    %% Local deploy folder editbox
    uitext(PathsPanel,'Local deploy folder:',[10 110 110 17]);
    LocalDeployEdit = uiedit(PathsPanel,'',[120 110 275 21], ...
        'HorizontalAlignment','left','Tag','LocalDeployFolder');

    %% Default job editbox
    uitext(PathsPanel,'Default job folder:',[10 85 110 17]);
    DefJobFolderEdit = uiedit(PathsPanel,'',[120 85 275 21], ...
        'HorizontalAlignment','left','Tag','DefJobFolder');
    
    %% Submit command editbox
    uitext(PathsPanel,'Submit command:',[10 60 110 17]);
    SubmitCommEdit = uiedit(PathsPanel,'',[120 60 275 21], ...
        'HorizontalAlignment','left','Tag','SubmitCommand');

    %% Run script editbox
    uitext(PathsPanel,'Run command/script:',[10 35 110 17]);
    RunScriptEdit = uiedit(PathsPanel,'',[120 35 275 21], ...
        'HorizontalAlignment','left','Tag','RunScript');

    %% Pushbuttons for interactive folder retrieval
    GetDirButton(1) = uipushbutton(PathsPanel,'...',[400 135 21 21], ...
        {@GetDir,PuTTyPathEdit,'PuTTy Path...'});
    GetDirButton(2) = uipushbutton(PathsPanel,'...',[400 110 21 21], ...
        {@GetDir,LocalDeployEdit,'Local Deploy Folder...'});
    GetDirButton(3) = uipushbutton(PathsPanel,'...',[400 85 21 21], ...
        {@GetDir,DefJobFolderEdit,'Default Job Folder...'});
    GetDirButton(4) = uipushbutton(PathsPanel,'...',[400 60 21 21], ...
        {@GetFile,SubmitCommEdit,'Queue Submit Command...'});
    GetDirButton(5) = uipushbutton(PathsPanel,'...',[400 35 21 21], ...
        {@GetFile,RunScriptEdit,'Command/Script File...'});
    
    %% checkbox for existing script
    RunScriptCheckbox = uicheckbox(PathsPanel, ...
        'Check the job folder for existing script first', ...
        0,[120 10 275 21],'');
        
    %% Callback when settings in radiobuttons are changed
    %% tells which things should be greyed out with which settings
    function ChangeCallback(src,eventdata)
        switch get(get(src,'SelectedObject'),'String')
            case 'Local'
                set(get(HostPanel,'Children'),'Enable','off');
                set([PuTTyPathEdit,LocalDeployEdit, ...
                     GetDirButton(1:2)],'Enable','off');
                set(GetDirButton([3,5]),'Enable','on');
                if get(ProcessingOptions(2),'Value') == 1,
                    set(GetDirButton(4),'Enable','on');
                end;
            case 'Remote'
                set(get(HostPanel,'Children'),'Enable','on');
                set([PuTTyPathEdit,LocalDeployEdit, ...
                     GetDirButton(1:2)],'Enable','on');
                set(GetDirButton(3:5),'Enable','off');
            case 'Run immediately'
                set([SubmitCommEdit,GetDirButton(4)],'Enable','off');
            case 'Submit to queue'
                set(SubmitCommEdit,'Enable','on');
                if get(LocationOptions(1),'Value') == 1,
                    set(GetDirButton(4),'Enable','on');
                end;
            case 'Ready to use'
                set(RunScriptCheckbox,'Enable','off');
            case 'Edit from template'
                set(RunScriptCheckbox,'Enable','on');
        end;
    end

    %% handles of the controls
    Editboxes = findobj(EngineDialog,'Style','edit')';
    RadioGroups = [LocationGroup,ProcessingGroup,ScriptGroup];

    %% control buttons
    SaveButton   = Button(EngineDialog,'Save'  ,2,@EngineSave );
    CancelButton = Button(EngineDialog,'Cancel',1,@CloseFigure);

    %% default behavior of Enter and Escape on edit controls
    set(Editboxes,'KeyPressFcn',{@RetEscPressed,@EngineSave,@CloseFigure});

    %% filling controls if editing and not new
    if strcmp(mode,'modify'),
        
        %% editboxes
        for h = Editboxes,
            set(h,'String',Engines(EnginePointer).(get(h,'Tag')));
        end;
        
        %% radiobuttons
        for group = RadioGroups,
            for button = findobj(group,'Style','radiobutton')',
                if strcmp( Engines(EnginePointer).(get(group,'Tag')), ...
                   get(button,'String') ),
                    set(group,'SelectedObject',button);
                end;
            end;
        end;
        
        %% checkbox
        set(RunScriptCheckbox,'Value',Engines(EnginePointer).CheckScript);
        
    end;

    %% set the correct greying
    for group = RadioGroups,
        ChangeCallback(group,[]);
    end;

    uicontrol(EngineNameEdit); %% give focus to the first field
    set(EngineDialog,'Visible','on'); %% display figure

    %% Saving the parameters
    function EngineSave(src,eventdata);

        %% active strings except username and password cannot be empty
        for h = [EngineNameEdit,HostNameEdit,PuTTyPathEdit,LocalDeployEdit, ...
                 DefJobFolderEdit,SubmitCommEdit,RunScriptEdit],
             if strcmp(get(h,'Enable'),'on') && isempty(get(h,'String')),
                 errordlg('Input cannot be empty','Bad input','modal');
                 pause(1);
                 uicontrol(h);
                 return;
             end;
        end;
        
        %% If the appropriate getdir pushbutton is active (which means the 
        %% folder is local)), then check if the folder is valid
        DirEdits = [PuTTyPathEdit,LocalDeployEdit,DefJobFolderEdit];
        for i = 1:2,
            if strcmp(get(GetDirButton(i),'Enable'),'on') && ...
                ~isdir(get(DirEdits(i),'String')),
                errordlg('Path must point to a valid local folder', ...
                         'Bad input','modal');
                pause(1);
                uicontrol(DirEdits(i));
                return;
            end;
        end;
        
        %% check if the script file exists
        if strcmp(get(GetDirButton(5),'Enable'),'on') && ...
           ~exist(get(RunScriptEdit,'String'),'file'),
            errordlg('Script file not found', ...
                     'Bad input','modal');
            pause(1);
            uicontrol(RunScriptEdit);
            return;
        end;
        
        %% emptying inactive controls
        for h = Editboxes,
            if strcmp(get(h,'Enable'),'off'), set(h,'String',''); end;
        end;

        %% if new engine, set the pointer after the last engine
        if strcmp(mode,'add'), EnginePointer = length(Engines) + 1; end;

        %% updating (or adding) the engine definition
        for h = Editboxes,
            Engines(EnginePointer).(get(h,'Tag')) = get(h,'String');
        end;
        for h = RadioGroups,
            Engines(EnginePointer).(get(h,'Tag')) ...
                = get(get(h,'SelectedObject'),'String');
        end;
        Engines(EnginePointer).CheckScript = get(RunScriptCheckbox,'Value');

        %% including the new string in the drop-down menu
        EngineString = get(EngineDropMenu,'String'); % retrieve
        EngineString = { EngineString{1:EnginePointer-1}, ...
                         Engines(EnginePointer).Name, ...
                         EngineString{EnginePointer+1:end} }; % update
        set(EngineDropMenu,'String',EngineString, ...
                           'Value',EnginePointer, ...
                           'UserData',Engines); % save

        close(EngineDialog);
    
    end;
    
end



%% Exit menu
function Exit(src,eventdata)
        
    if ishandle(GeometryFig), delete(GeometryFig); end;
    delete(MainFig);
        
end



%%%%%%%%%%%%%%%%%
%% WIRE DIALOG %%
%%%%%%%%%%%%%%%%%
function Wire(src,eventdata,GeomPointer);
        
    %% main dialog window
    WireDialog = dialog( ...
        'Visible','off', ...
        'Name','Wire', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 260], ...
        'KeyPressFcn',{@RetEscPressed,@WireStore,@CloseFigure} );

    %% local geometry panel
    WireGeometry = uipanel( ...
        'Parent',WireDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 120], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    %% texts
    uitext(WireGeometry,'x',[ 80 91 60 16],'HorizontalAlignment','center');
    uitext(WireGeometry,'y',[145 91 60 16],'HorizontalAlignment','center');
    uitext(WireGeometry,'z',[210 91 60 16],'HorizontalAlignment','center');
    uitext(WireGeometry,'Endpoint 1:',[10 70 70 17]);
    uitext(WireGeometry,'Endpoint 2:',[10 45 70 17]);
    uitext(WireGeometry,'Units:'     ,[10 10 40 17]);

    %% uicontrols
    x1 = uiedit(WireGeometry,'',[ 80 70 60 21]);
    y1 = uiedit(WireGeometry,'',[145 70 60 21]);
    z1 = uiedit(WireGeometry,'',[210 70 60 21]);
    x2 = uiedit(WireGeometry,'',[ 80 45 60 21]);
    y2 = uiedit(WireGeometry,'',[145 45 60 21]);
    z2 = uiedit(WireGeometry,'',[210 45 60 21]);
    UnitsUI = uipopupmenu(WireGeometry,Units,1,[50 10 55 21]);

    data_handles = [x1 y1 z1 x2 y2 z2 UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@WireStore,@CloseFigure});

    %% Comment
    uitext(WireDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(WireDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(WireDialog,'OK',2,@WireStore);
    CancelButton = Button(WireDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(WireDialog,OKButton,CancelButton); end;
    
    uicontrol(x1); %% give focus to the first field
    set(WireDialog,'Visible','on'); %% display figure

    function WireStore(src,eventdata);

        %% get endpoints of the wire
        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];

        check = ~isfinite([point1 point2]) | imag([point1 point2])~=0;
        if any(check),
            errordlg('Endpoint coordinates must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;
        
        %% check for zero-length wire (after staircasing)
        if round(point1)==round(point2),
            errordlg('Length of the wire is too small', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(x1);
            return;
        end;

        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','wire', ...
            'point1',point1, ...
            'point2',point2, ...
            'units',get(UnitsUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character in comment
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'wire' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayWire(GeomPointer,'k');

        close(WireDialog);

    end %% end WireStore

end %% end Wire



%%%%%%%%%%%%%%%%%%%%%
%% TRIANGLE DIALOG %%
%%%%%%%%%%%%%%%%%%%%%
function Triangle(src,eventdata,GeomPointer);
        
    %% main dialog window
    TriangleDialog = dialog( ...
        'Visible','off', ...
        'Name','Triangle', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 285], ...
        'KeyPressFcn',{@RetEscPressed,@TriangleStore,@CloseFigure} );

    %% local geometry panel
    TriangleGeometry = uipanel( ...
        'Parent',TriangleDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 145], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    %% texts
    uitext(TriangleGeometry,'x',[ 80 116 60 16],'HorizontalAlignment','center');
    uitext(TriangleGeometry,'y',[145 116 60 16],'HorizontalAlignment','center');
    uitext(TriangleGeometry,'z',[210 116 60 16],'HorizontalAlignment','center');
    uitext(TriangleGeometry,'Vertex 1:',[10 95 70 17]);
    uitext(TriangleGeometry,'Vertex 2:',[10 70 70 17]);
    uitext(TriangleGeometry,'Vertex 3:',[10 45 70 17]);
    uitext(TriangleGeometry,'Units:'   ,[10 10 40 17]);

    %% uicontrols
    x1 = uiedit(TriangleGeometry,'',[ 80 95 60 21]);
    y1 = uiedit(TriangleGeometry,'',[145 95 60 21]);
    z1 = uiedit(TriangleGeometry,'',[210 95 60 21]);
    x2 = uiedit(TriangleGeometry,'',[ 80 70 60 21]);
    y2 = uiedit(TriangleGeometry,'',[145 70 60 21]);
    z2 = uiedit(TriangleGeometry,'',[210 70 60 21]);
    x3 = uiedit(TriangleGeometry,'',[ 80 45 60 21]);
    y3 = uiedit(TriangleGeometry,'',[145 45 60 21]);
    z3 = uiedit(TriangleGeometry,'',[210 45 60 21]);
    UnitsUI = uipopupmenu(TriangleGeometry,Units,1,[50 10 55 21]);

    data_handles = [x1 y1 z1 x2 y2 z2 x3 y3 z3 UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@TriangleStore,@CloseFigure});

    %% Comment
    uitext(TriangleDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(TriangleDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(TriangleDialog,'OK',2,@TriangleStore);
    CancelButton = Button(TriangleDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(x3,'String',num2str(v.Geometry{GeomPointer}.point3(1)));
        set(y3,'String',num2str(v.Geometry{GeomPointer}.point3(2)));
        set(z3,'String',num2str(v.Geometry{GeomPointer}.point3(3)));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(TriangleDialog,OKButton,CancelButton); end;
    
    uicontrol(x1); %% give focus to the first field
    set(TriangleDialog,'Visible','on'); %% display figure

    function TriangleStore(src,eventdata);

        %% get vertex points of the triangle
        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];
        point3 = [str2double(get(x3,'String')), ...
                  str2double(get(y3,'String')), ...
                  str2double(get(z3,'String'))];

        check = ~isfinite([point1 point2 point3]) | imag([point1 point2 point3])~=0;
        if any(check),
            errordlg('Vertex coordinates must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;
        
        %% compute the faces and vertices and check if renderable
        [faces,vertices,patch_count] = ...
            StaircasePlane(point1,point2,point3,'triangle');
        if isempty(faces),
            errordlg('Unable to render the plane - too small size', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(x1);
            return;
        end;

        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','triangle', ...
            'point1',point1, ...
            'point2',point2, ...
            'point3',point3, ...
            'units',get(UnitsUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character in comment
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'triangle' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayPlate(GeomPointer,faces,vertices,patch_count);

        close(TriangleDialog);

    end %% end TriangleStore

end %% end Triangle



%%%%%%%%%%%%%%%%%%%%%%
%% RECTANGLE DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%
function Rectangle(src,eventdata,GeomPointer);
        
    %% main dialog window
    RectangleDialog = dialog( ...
        'Visible','off', ...
        'Name','Rectangle', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 285], ...
        'KeyPressFcn',{@RetEscPressed,@RectangleStore,@CloseFigure} );

    %% local geometry panel
    RectangleGeometry = uipanel( ...
        'Parent',RectangleDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 145], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    %% texts
    uitext(RectangleGeometry,'x',[ 80 116 60 16],'HorizontalAlignment','center');
    uitext(RectangleGeometry,'y',[145 116 60 16],'HorizontalAlignment','center');
    uitext(RectangleGeometry,'z',[210 116 60 16],'HorizontalAlignment','center');
    uitext(RectangleGeometry,'Basepoint 1:' ,[10 95 70 17]);
    uitext(RectangleGeometry,'Basepoint 2:' ,[10 70 70 17]);
    uitext(RectangleGeometry,'Extent point:',[10 45 70 17]);
    uitext(RectangleGeometry,'Units:'       ,[10 10 40 17]);

    %% uicontrols
    x1 = uiedit(RectangleGeometry,'',[ 80 95 60 21]);
    y1 = uiedit(RectangleGeometry,'',[145 95 60 21]);
    z1 = uiedit(RectangleGeometry,'',[210 95 60 21]);
    x2 = uiedit(RectangleGeometry,'',[ 80 70 60 21]);
    y2 = uiedit(RectangleGeometry,'',[145 70 60 21]);
    z2 = uiedit(RectangleGeometry,'',[210 70 60 21]);
    x3 = uiedit(RectangleGeometry,'',[ 80 45 60 21]);
    y3 = uiedit(RectangleGeometry,'',[145 45 60 21]);
    z3 = uiedit(RectangleGeometry,'',[210 45 60 21]);
    UnitsUI = uipopupmenu(RectangleGeometry,Units,1,[50 10 55 21]);

    data_handles = [x1 y1 z1 x2 y2 z2 x3 y3 z3 UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@RectangleStore,@CloseFigure});

    %% Comment
    uitext(RectangleDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(RectangleDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(RectangleDialog,'OK',2,@RectangleStore);
    CancelButton = Button(RectangleDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(x3,'String',num2str(v.Geometry{GeomPointer}.point3(1)));
        set(y3,'String',num2str(v.Geometry{GeomPointer}.point3(2)));
        set(z3,'String',num2str(v.Geometry{GeomPointer}.point3(3)));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(RectangleDialog,OKButton,CancelButton); end;

    uicontrol(x1); %% give focus to the first field
    set(RectangleDialog,'Visible','on'); %% display figure

    function RectangleStore(src,eventdata);

        %% get points of the rectangle
        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];
        point3 = [str2double(get(x3,'String')), ...
                  str2double(get(y3,'String')), ...
                  str2double(get(z3,'String'))];

        check = ~isfinite([point1 point2 point3]) | imag([point1 point2 point3])~=0;
        if any(check),
            errordlg('Point coordinates must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;
        
        %% compute the faces and vertices and check if renderable
        [faces,vertices,patch_count] = ...
            StaircasePlane(point1,point2,point3,'rectangle');
        if isempty(faces),
            errordlg('Unable to render the plane - too small size', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(x1);
            return;
        end;

        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','rectangle', ...
            'point1',point1, ...
            'point2',point2, ...
            'point3',point3, ...
            'units',get(UnitsUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character in comment
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'rectangle' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayPlate(GeomPointer,faces,vertices,patch_count);

        close(RectangleDialog);

    end %% end RectangleStore

end %% end Rectangle



%%%%%%%%%%%%%%%%%%
%% BRICK DIALOG %%
%%%%%%%%%%%%%%%%%%
function Brick(src,eventdata,GeomPointer);
        
    %% main dialog window
    BrickDialog = dialog( ...
        'Visible','off', ...
        'Name','Brick', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 400], ...
        'KeyPressFcn',{@RetEscPressed,@BrickStore,@CloseFigure} );

    %% local geometry panel
    BrickGeometry = uipanel( ...
        'Parent',BrickDialog, ...
        'Units','pixels', ...
        'Position',[10 275 285 120], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );
    
    %% properties panel
    BrickProperties = uipanel( ...
        'Parent',BrickDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 135], ...
        'ForegroundColor','b', ...
        'Title','Properties' );

    %% texts
    uitext(BrickGeometry,'x',[ 80 91 60 16],'HorizontalAlignment','center');
    uitext(BrickGeometry,'y',[145 91 60 16],'HorizontalAlignment','center');
    uitext(BrickGeometry,'z',[210 91 60 16],'HorizontalAlignment','center');
    uitext(BrickGeometry,'Corner 1:',[10 70 70 17]);
    uitext(BrickGeometry,'Corner 2:',[10 45 70 17]);
    uitext(BrickGeometry,'Units:'   ,[10 10 40 17]);

    uitext(BrickProperties,'Permittivity:',[10 65 70 17]);
    uitext(BrickProperties,'Conductivity:',[10 40 70 17]);
    uitext(BrickProperties,'Mass density:',[10 15 70 17]);
    uitext(BrickProperties,'(relative)',[145 65 70 17]);
    uitext(BrickProperties,'S/m'       ,[145 40 70 17]);
    uitext(BrickProperties,'kg/m^3'    ,[145 15 70 17]);

    %% uicontrols
    x1 = uiedit(BrickGeometry,'',[ 80 70 60 21]);
    y1 = uiedit(BrickGeometry,'',[145 70 60 21]);
    z1 = uiedit(BrickGeometry,'',[210 70 60 21]);
    x2 = uiedit(BrickGeometry,'',[ 80 45 60 21]);
    y2 = uiedit(BrickGeometry,'',[145 45 60 21]);
    z2 = uiedit(BrickGeometry,'',[210 45 60 21]);
    UnitsUI = uipopupmenu(BrickGeometry,Units,1,[50 10 55 21]);

    PEC_Checkbox = uicheckbox(BrickProperties, ...
        'Perfect Electric Conductor',0,[10 95 270 21],@PEC_Callback);
    function PEC_Callback(src,eventdata)
        PEC = get(src,'Value')==1;
        if PEC, set([EpsUI SigUI RhoUI],'Enable','off');
        else    set([EpsUI SigUI RhoUI],'Enable','on');
        end;
    end
    EpsUI = uiedit(BrickProperties,'1'               ,[80 65 60 21]);
    SigUI = uiedit(BrickProperties,'0'               ,[80 40 60 21]);
    RhoUI = uiedit(BrickProperties,DefaultMassDensity,[80 15 60 21]);
    
    data_handles = [x1 y1 z1 x2 y2 z2 UnitsUI PEC_Checkbox EpsUI SigUI RhoUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@BrickStore,@CloseFigure});

    %% Comment
    uitext(BrickDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(BrickDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(BrickDialog,'OK',2,@BrickStore);
    CancelButton = Button(BrickDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(PEC_Checkbox,'Value',v.Geometry{GeomPointer}.PEC);
        set(EpsUI,'String',num2str(v.Geometry{GeomPointer}.permittivity));
        set(SigUI,'String',num2str(v.Geometry{GeomPointer}.conductivity));
        set(RhoUI,'String',num2str(v.Geometry{GeomPointer}.massdensity));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
        if v.Geometry{GeomPointer}.PEC, set([EpsUI SigUI RhoUI],'Enable','off'); end;
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(BrickDialog,OKButton,CancelButton); end;

    uicontrol(x1); %% give focus to the first field
    set(BrickDialog,'Visible','on'); %% display figure

    function BrickStore(src,eventdata);

        %% get corners of the brick
        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];
        
        PEC = get(PEC_Checkbox,'Value')==1;
        if PEC, permittivity = 1;
                conductivity = 0;
                massdensity  = DefaultMassDensity;
        else permittivity = str2double(get(EpsUI,'String'));
             conductivity = str2double(get(SigUI,'String'));
             massdensity  = str2double(get(RhoUI,'String'));
        end;

        check = ~isfinite([point1 point2 0 0 permittivity conductivity massdensity]) | ...
                     imag([point1 point2 0 0 permittivity conductivity massdensity])~=0;
        if any(check),
            errordlg('Brick parameters must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% check for zero-size brick (after rounding)
        if any( round(point1) == round(point2) ),
            errordlg('Size of the brick is too small', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(x1);
            return;
        end;
        
        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;
        
        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;
 
        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','brick', ...
            'point1',point1, ...
            'point2',point2, ...
            'PEC',PEC, ...
            'permittivity',permittivity, ...
            'conductivity',conductivity, ...
            'massdensity',massdensity, ...
            'units',get(UnitsUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'brick' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayBrick(GeomPointer);

        close(BrickDialog);

    end %% end BrickStore

end %% end Brick



%%%%%%%%%%%%%%%%%%%
%% SOURCE DIALOG %%
%%%%%%%%%%%%%%%%%%%
function Source(src,eventdata,GeomPointer);
        
    %% main dialog window
    SourceDialog = dialog( ...
        'Visible','off', ...
        'Name','Source', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 330], ...
        'KeyPressFcn',{@RetEscPressed,@SourceStore,@CloseFigure} );

    %% local geometry panel
    SourceGeometry = uipanel( ...
        'Parent',SourceDialog, ...
        'Units','pixels', ...
        'Position',[10 195 285 130], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    SourceProperties = uipanel( ...
        'Parent',SourceDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 55], ...
        'ForegroundColor','b', ...
        'Title','Properties');

    %% texts
    uitext(SourceGeometry,'x',[ 80 101 60 16],'HorizontalAlignment','center');
    uitext(SourceGeometry,'y',[145 101 60 16],'HorizontalAlignment','center');
    uitext(SourceGeometry,'z',[210 101 60 16],'HorizontalAlignment','center');
    uitext(SourceGeometry,'Ground point:',[10 80 70 17]);
    uitext(SourceGeometry,'Direction:',[10 50 70 17]);
    uitext(SourceGeometry,'Units:'    ,[10 10 40 17]);

    uitext(SourceProperties,'Internal resistance:',[10 15 100 17]);
    uitext(SourceProperties,'Ohm',[175 15 70 17]);
    
    %% uicontrols
    x = uiedit(SourceGeometry,'',[ 80 80 60 21]);
    y = uiedit(SourceGeometry,'',[145 80 60 21]);
    z = uiedit(SourceGeometry,'',[210 80 60 21]);
    DirectionUI = uipopupmenu(SourceGeometry,Directions,1,[80 50 50 21]);
    UnitsUI = uipopupmenu(SourceGeometry,Units,1,[50 10 55 21]);

    ResistUI = uiedit(SourceProperties,num2str(DefaultSourceResistance),[110 15 60 21]);
    
    data_handles = [x y z DirectionUI UnitsUI ResistUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@SourceStore,@CloseFigure});

    %% Comment
    uitext(SourceDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(SourceDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(SourceDialog,'OK',2,@SourceStore);
    CancelButton = Button(SourceDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x,'String',num2str(v.Geometry{GeomPointer}.gndpoint(1)));
        set(y,'String',num2str(v.Geometry{GeomPointer}.gndpoint(2)));
        set(z,'String',num2str(v.Geometry{GeomPointer}.gndpoint(3)));
        set(DirectionUI,'Value',v.Geometry{GeomPointer}.direction);
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(ResistUI,'String',num2str(v.Geometry{GeomPointer}.resistance));
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(SourceDialog,OKButton,CancelButton); end;

    uicontrol(x); %% give focus to the first field
    set(SourceDialog,'Visible','on'); %% display figure

    function SourceStore(src,eventdata);

        %% get groundpoint of the source
        gndpoint = [str2double(get(x,'String')), ...
                    str2double(get(y,'String')), ...
                    str2double(get(z,'String'))];
        resistance = str2double(get(ResistUI,'String'));

        check = ~isfinite([gndpoint 0 0 resistance]) | ...
                imag([gndpoint 0 0 resistance])~=0;
        if any(check),
            errordlg('Inserted values must be real numeric values', ...
                'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% rounding in case that the coordinates are not integers
        gndpoint = round(gndpoint);
        
        %% second vertex of the source
        point2 = gndpoint;
        switch get(DirectionUI,'Value')
            case 1, point2(1) = point2(1) + 1;
            case 2, point2(1) = point2(1) - 1;
            case 3, point2(2) = point2(2) + 1;
            case 4, point2(2) = point2(2) - 1;
            case 5, point2(3) = point2(3) + 1;
            case 6, point2(3) = point2(3) - 1;
        end;

        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','source', ...
            'gndpoint',gndpoint, ...
            'point2',point2, ...
            'units',get(UnitsUI,'Value'), ...
            'direction',get(DirectionUI,'Value'), ...
            'resistance',resistance, ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'source' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplaySource(GeomPointer);

        close(SourceDialog);

    end %% end SourceStore

end %% end source



%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIELD SOURCE DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function FieldSource(src,eventdata,GeomPointer);
        
    %% main dialog window
    FieldSourceDialog = dialog( ...
        'Visible','off', ...
        'Name','Field Source', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 340], ...
        'KeyPressFcn',{@RetEscPressed,@FieldSourceStore,@CloseFigure} );

    %% local geometry panel
    FieldSourceGeometry = uipanel( ...
        'Parent',FieldSourceDialog, ...
        'Units','pixels', ...
        'Position',[10 220 285 115], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    FieldSourceProperties = uipanel( ...
        'Parent',FieldSourceDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 80], ...
        'ForegroundColor','b', ...
        'Title','Properties');

    %% texts
    uitext(FieldSourceGeometry,'x',[ 80 86 60 16],'HorizontalAlignment','center');
    uitext(FieldSourceGeometry,'y',[145 86 60 16],'HorizontalAlignment','center');
    uitext(FieldSourceGeometry,'z',[210 86 60 16],'HorizontalAlignment','center');
    uitext(FieldSourceGeometry,'Corner 1:' ,[ 10 65 70 17]);
    uitext(FieldSourceGeometry,'Corner 2:' ,[ 10 40 70 17]);
    uitext(FieldSourceGeometry,'Component:',[ 10 10 70 17]);
    uitext(FieldSourceGeometry,'Units:'    ,[180 10 40 17]);
    
    uitext(FieldSourceProperties,'Field update type:',[10 40 95 17]);
    uitext(FieldSourceProperties,'Data:',[10 10 40 17]);
    
    %% button "Load"
    LoadFieldPush = uipushbutton(FieldSourceProperties,'Load...',[195 10 75 23], ...
                                 @LoadFieldCallback);

    function LoadFieldCallback(src,eventdata)

        %% choose file - uigetfile
        [FileName,PathName] = uigetfile( ...
            {'*.mat','Matlab Data Files (*.mat)'; ...
             '*','All Files (*.*)'}, ...
             'Load Field Source',[WorkDir,filesep]);

        if isequal(FileName,0) || isequal(PathName,0), return; end;
        WorkDir = PathName;

        try
            ListOfVariables = who('-file',fullfile(PathName,FileName));
        catch
            errordlg(['Error when opening file ',FileName],'Load Error');
            return;
        end;

        %% Variable selection
        while 1,
            
            [pointer,ok] = listdlg('ListString',ListOfVariables, ...
                                   'ListSize',[160 160], ...
                                   'Name','Field Source', ...
                                   'PromptString','Select Variable:');
            if not(ok), return; end;

            try
                TempData = load('-mat',fullfile(PathName,FileName), ...
                                       ListOfVariables{pointer});
            catch
                errordlg(['Error when opening file ',FileName],'Load Error');
                return;
            end;
            TempData = TempData.(ListOfVariables{pointer});
            
            %% check if numeric
            if isnumeric(TempData), break; end;

            %% if not, give error
            errordlg('Selected variable is not numeric!','Variable Error');
            uiwait;
            
        end;
        
        FieldSourceData = TempData;

        %% show info before the pushbutton - class and size
        set(DataText,'String',VarString(FieldSourceData));

    end

    %% function generating formatted string describing variable
    function var_string = VarString(var)
        
        size_var = size(var);
        if size_var < 3, size_var(3) = 1; end;
        var_string = num2str(size_var(1));
        for i = 2:length(size_var),
            var_string = [var_string,' x ',num2str(size_var(i))];
        end;
        var_string = [class(var),' (',var_string,')'];
        
    end

    %% uicontrols
    x1 = uiedit(FieldSourceGeometry,'',[ 80 65 60 21]);
    y1 = uiedit(FieldSourceGeometry,'',[145 65 60 21]);
    z1 = uiedit(FieldSourceGeometry,'',[210 65 60 21]);
    x2 = uiedit(FieldSourceGeometry,'',[ 80 40 60 21]);
    y2 = uiedit(FieldSourceGeometry,'',[145 40 60 21]);
    z2 = uiedit(FieldSourceGeometry,'',[210 40 60 21]);
    UnitsUI     = uipopupmenu(FieldSourceGeometry,Units     ,1,[220 10 50 21]);
    ComponentUI = uipopupmenu(FieldSourceGeometry,Components,1,[ 80 10 50 21]);
    UpdateUI    = uipopupmenu(FieldSourceProperties,{'replace','add'},1,[105 40 80 21]);
    DataText = uitext(FieldSourceProperties,'<empty>',[45 10 140 17],'FontWeight','bold');

    data_handles = [x1 y1 z1 x2 y2 z2 UnitsUI ComponentUI];

    FieldSourceData = [];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@FieldSourceStore,@CloseFigure});

    %% Comment
    uitext(FieldSourceDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(FieldSourceDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(FieldSourceDialog,'OK',2,@FieldSourceStore);
    CancelButton = Button(FieldSourceDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,

        FieldSourceData = v.Geometry{GeomPointer}.data;

        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(ComponentUI,'Value',v.Geometry{GeomPointer}.component);
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(UpdateUI,'Value',v.Geometry{GeomPointer}.update);
        set(DataText,'String',VarString(FieldSourceData));
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
        
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(FieldSourceDialog,OKButton,CancelButton); end;

    uicontrol(x1); %% give focus to the first field
    set(FieldSourceDialog,'Visible','on'); %% display figure

    function FieldSourceStore(src,eventdata);

        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];

        check = ~isfinite([point1 point2]) | imag([point1 point2])~=0;
        if any(check),
            errordlg('Corner coordinates must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;
        
        %% rounding in case that the coordinates are not integers
        point1 = round(point1);
        point2 = round(point2);
        
        %% check the size of the data array
        data_size = size(FieldSourceData);
        if length(data_size) < 3, data_size(3) = 1; end;
        if not(isequal(data_size,abs(point2 - point1) + 1)),
            errordlg('Data size does not match', ...
                     'Bad input','modal');
            uiwait;
            uicontrol(x1);
            return;
        end;
        
        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','field source', ...
            'point1',point1, ...
            'point2',point2, ...
            'component',get(ComponentUI,'Value'), ...
            'units',get(UnitsUI,'Value'), ...
            'update',get(UpdateUI,'Value'), ...
            'data',FieldSourceData, ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character in comment
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'field source' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString,'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayFieldSource(GeomPointer);

        close(FieldSourceDialog);

    end %% end FieldSourceStore

end %% end field source



%%%%%%%%%%%%%%%%%%%%%
%% RESISTOR DIALOG %%
%%%%%%%%%%%%%%%%%%%%%
function Resistor(src,eventdata,GeomPointer);
        
    %% main dialog window
    ResistorDialog = dialog( ...
        'Visible','off', ...
        'Name','Resistor', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 320], ...
        'KeyPressFcn',{@RetEscPressed,@ResistorStore,@CloseFigure} );

    %% local geometry panel
    ResistorGeometry = uipanel( ...
        'Parent',ResistorDialog, ...
        'Units','pixels', ...
        'Position',[10 195 285 120], ...
        'ForegroundColor','b', ...
        'Title','Geometry');
    
    ResistorProperties = uipanel( ...
        'Parent',ResistorDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 55], ...
        'ForegroundColor','b', ...
        'Title','Properties');

    %% texts
    uitext(ResistorGeometry,'x',[ 80 91 60 16],'HorizontalAlignment','center');
    uitext(ResistorGeometry,'y',[145 91 60 16],'HorizontalAlignment','center');
    uitext(ResistorGeometry,'z',[210 91 60 16],'HorizontalAlignment','center');
    uitext(ResistorGeometry,'Endpoint 1:',[10 70 70 17]);
    uitext(ResistorGeometry,'Endpoint 2:',[10 45 70 17]);
    uitext(ResistorGeometry,'Units:'     ,[10 10 40 17]);

    uitext(ResistorProperties,'Resistance:',[10 15 70 17]);
    uitext(ResistorProperties,'Ohm',[145 15 70 17]);
    
    %% uicontrols
    x1 = uiedit(ResistorGeometry,'',[ 80 70 60 21]);
    y1 = uiedit(ResistorGeometry,'',[145 70 60 21]);
    z1 = uiedit(ResistorGeometry,'',[210 70 60 21]);
    x2 = uiedit(ResistorGeometry,'',[ 80 45 60 21]);
    y2 = uiedit(ResistorGeometry,'',[145 45 60 21]);
    z2 = uiedit(ResistorGeometry,'',[210 45 60 21]);
    UnitsUI = uipopupmenu(ResistorGeometry,Units,1,[50 10 55 21]);

    ResistUI = uiedit(ResistorProperties,'',[80 15 60 21]);
    
    data_handles = [x1 y1 z1 x2 y2 z2 UnitsUI ResistUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@ResistorStore,@CloseFigure});

    %% Comment
    uitext(ResistorDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(ResistorDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(ResistorDialog,'OK',2,@ResistorStore);
    CancelButton = Button(ResistorDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(ResistUI,'String',num2str(v.Geometry{GeomPointer}.resistance));
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(ResistorDialog,OKButton,CancelButton); end;

    uicontrol(x1); %% give focus to the first field
    set(ResistorDialog,'Visible','on'); %% display figure

    function ResistorStore(src,eventdata);

        %% get endpoints and resistance
        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];
        resistance = str2double(get(ResistUI,'String'));

        check = ~isfinite([point1 point2 0 resistance]) | ...
                imag([point1 point2 0 resistance])~=0;
        if any(check),
            errordlg('Resistor parameters must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% rounding in case that the coordinates are not integers
        point1 = round(point1);
        point2 = round(point2);
        
        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','resistor', ...
            'point1',point1, ...
            'point2',point2, ...
            'resistance',resistance, ...
            'units',get(UnitsUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle);

        %% adding a space for every newline character
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'resistor' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
            ListLine, ...
            v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayWire(GeomPointer,'b');
        
        close(ResistorDialog);

    end %% end ResistorStore

end %% end resistor



%%%%%%%%%%%%%%%%%%%%%%
%% CAPACITOR DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%
function Capacitor(src,eventdata,GeomPointer);
        
    %% main dialog window
    CapacitorDialog = dialog( ...
        'Visible','off', ...
        'Name','Capacitor', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 320], ...
        'KeyPressFcn',{@RetEscPressed,@CapacitorStore,@CloseFigure} );

    %% local geometry panel
    CapacitorGeometry = uipanel( ...
        'Parent',CapacitorDialog, ...
        'Units','pixels', ...
        'Position',[10 195 285 120], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    CapacitorProperties = uipanel( ...
        'Parent',CapacitorDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 55], ...
        'ForegroundColor','b', ...
        'Title','Properties');

    %% texts
    uitext(CapacitorGeometry,'x',[ 80 91 60 16],'HorizontalAlignment','center');
    uitext(CapacitorGeometry,'y',[145 91 60 16],'HorizontalAlignment','center');
    uitext(CapacitorGeometry,'z',[210 91 60 16],'HorizontalAlignment','center');
    uitext(CapacitorGeometry,'Endpoint 1:',[10 70 70 17]);
    uitext(CapacitorGeometry,'Endpoint 2:',[10 45 70 17]);
    uitext(CapacitorGeometry,'Units:'     ,[10 10 40 17]);

    uitext(CapacitorProperties,'Capacitance:',[10 15 70 17]);
    uitext(CapacitorProperties,'nF',[145 15 70 17]);

    %% uicontrols
    x1 = uiedit(CapacitorGeometry,'',[ 80 70 60 21]);
    y1 = uiedit(CapacitorGeometry,'',[145 70 60 21]);
    z1 = uiedit(CapacitorGeometry,'',[210 70 60 21]);
    x2 = uiedit(CapacitorGeometry,'',[ 80 45 60 21]);
    y2 = uiedit(CapacitorGeometry,'',[145 45 60 21]);
    z2 = uiedit(CapacitorGeometry,'',[210 45 60 21]);
    UnitsUI = uipopupmenu(CapacitorGeometry,Units,1,[50 10 55 21]);

    CapacitUI = uiedit(CapacitorProperties,'',[80 15 60 21]);
    
    data_handles = [x1 y1 z1 x2 y2 z2 UnitsUI CapacitUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@CapacitorStore,@CloseFigure});

    %% Comment
    uitext(CapacitorDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(CapacitorDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(CapacitorDialog,'OK',2,@CapacitorStore);
    CancelButton = Button(CapacitorDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(x2,'String',num2str(v.Geometry{GeomPointer}.point2(1)));
        set(y2,'String',num2str(v.Geometry{GeomPointer}.point2(2)));
        set(z2,'String',num2str(v.Geometry{GeomPointer}.point2(3)));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CapacitUI,'String',num2str(v.Geometry{GeomPointer}.capacitance));
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(CapacitorDialog,OKButton,CancelButton); end;

    uicontrol(x1); %% give focus to the first field
    set(CapacitorDialog,'Visible','on'); %% display figure

    function CapacitorStore(src,eventdata);

        %% get endpoints and capacitance
        point1 = [str2double(get(x1,'String')), ...
                  str2double(get(y1,'String')), ...
                  str2double(get(z1,'String'))];
        point2 = [str2double(get(x2,'String')), ...
                  str2double(get(y2,'String')), ...
                  str2double(get(z2,'String'))];
        capacitance = str2double(get(CapacitUI,'String'));

        check = ~isfinite([point1 point2 0 capacitance]) | ...
                imag([point1 point2 0 capacitance])~=0;
        if any(check),
            errordlg('Capacitor parameters must be real numeric values', ...
                'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;
        
        %% rounding in case that the coordinates are not integers
        point1 = round(point1);
        point2 = round(point2);
        
        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','capacitor', ...
            'point1',point1, ...
            'point2',point2, ...
            'capacitance',capacitance, ...
            'units',get(UnitsUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle);

        %% adding a space for every newline character
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'capacitor' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayWire(GeomPointer,'g');
        
        close(CapacitorDialog);

    end %% end CapacitorStore

end %% end capacitor



%%%%%%%%%%%%%%%%%%%
%% PROBE DIALOG %%
%%%%%%%%%%%%%%%%%%%
function Probe(src,eventdata,GeomPointer);
        
    %% main dialog window
    ProbeDialog = dialog( ...
        'Visible','off', ...
        'Name','Probe', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 305 270], ...
        'KeyPressFcn',{@RetEscPressed,@ProbeStore,@CloseFigure} );

    %% local geometry panel
    ProbeGeometry = uipanel( ...
        'Parent',ProbeDialog, ...
        'Units','pixels', ...
        'Position',[10 135 285 130], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    %% texts
    uitext(ProbeGeometry,'x',[ 80 101 60 16],'HorizontalAlignment','center');
    uitext(ProbeGeometry,'y',[145 101 60 16],'HorizontalAlignment','center');
    uitext(ProbeGeometry,'z',[210 101 60 16],'HorizontalAlignment','center');
    uitext(ProbeGeometry,'Coordinates:',[10 80 70 17]);
    uitext(ProbeGeometry,'Component:'  ,[10 50 70 17]);
    uitext(ProbeGeometry,'Units:'      ,[10 10 40 17]);

    %% uicontrols
    x = uiedit(ProbeGeometry,'',[ 80 80 60 21]);
    y = uiedit(ProbeGeometry,'',[145 80 60 21]);
    z = uiedit(ProbeGeometry,'',[210 80 60 21]);
    ComponentUI = uipopupmenu(ProbeGeometry,Components,1,[80 50 50 21]);
    UnitsUI = uipopupmenu(ProbeGeometry,Units,1,[50 10 55 21]);
    
    data_handles = [x y z ComponentUI UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@ProbeStore,@CloseFigure});

    %% Comment
    uitext(ProbeDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(ProbeDialog,'',[10 55 285 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(ProbeDialog,'OK',2,@ProbeStore);
    CancelButton = Button(ProbeDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(ComponentUI,'Value',v.Geometry{GeomPointer}.component);
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(ProbeDialog,OKButton,CancelButton); end;

    uicontrol(x); %% give focus to the first field
    set(ProbeDialog,'Visible','on'); %% display figure

    function ProbeStore(src,eventdata);

        %% get coordinates of the probe
        point1 = [str2double(get(x,'String')), ...
                  str2double(get(y,'String')), ...
                  str2double(get(z,'String'))];

        check = ~isfinite(point1) | imag(point1)~=0;
        if any(check),
            errordlg('Inserted values must be real numeric values', ...
                'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% rounding in case that the coordinates are not integers
        point1 = round(point1);
        
        %% handle and visibility of the displayed entity, if exists
        if GeomPointer~=0,
            handle = v.Geometry{GeomPointer}.handle;
            visible = v.Geometry{GeomPointer}.visible;
        else
            handle = [];
            visible = true;
        end;

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = struct( ...
            'type','probe', ...
            'point1',point1, ...
            'units',get(UnitsUI,'Value'), ...
            'component',get(ComponentUI,'Value'), ...
            'comment',{get(CommentUI,'String')}, ...
            'visible',visible, ...
            'handle',handle );

        %% adding a space for every newline character
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'probe' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(handle), delete(handle); end;
        DisplayProbe(GeomPointer);

        close(ProbeDialog);

    end %% end ProbeStore

end %% end probe



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTERNAL OBJECT DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExtObj(src,eventdata,GeomPointer);
        
    if GeomPointer==0,

        while 1,
            [FileName,PathName] = uigetfile( ...
                {'*.mat','External Object MAT Files (*.mat)'; ...
                 '*','All Files (*.*)'}, ...
                'Open External Object',[ExtObjDir,filesep]);
            if isequal(FileName,0) || isequal(PathName,0), return; end;
            ExtObjFile = fullfile(PathName,FileName);
            ExtObjDir = PathName;
            try
                e = load('-mat',ExtObjFile);
            catch
                h = errordlg('File cannot be opened!','Open Error');
                waitfor(h);
                continue; %% open dialog again
            end;
            if ~isfield(e,{'VolArray','PlxArray','PlyArray','PlzArray'}),
                h = errordlg('Object not defined in file!','Error');
                waitfor(h);
                continue; %% open dialog again
            end;
            break;
        end;

        e.type = 'extobj';
        e.visible = true;
        e.handle = [];
        
        %% what to do, if some arrays do not exist
        if ~isfield(e,'VolArray'), e.VolArray = []; end;
        if ~isfield(e,'PlxArray'), e.PlxArray = []; end;
        if ~isfield(e,'PlyArray'), e.PlyArray = []; end;
        if ~isfield(e,'PlzArray'), e.PlzArray = []; end;
            
        if isempty(e.VolArray), VolMaxValue = 0; else VolMaxValue = max(e.VolArray(:)); end;
        if isempty(e.PlxArray), PlxMaxValue = 0; else PlxMaxValue = max(e.PlxArray(:)); end;
        if isempty(e.PlyArray), PlyMaxValue = 0; else PlyMaxValue = max(e.PlyArray(:)); end;
        if isempty(e.PlzArray), PlzMaxValue = 0; else PlzMaxValue = max(e.PlzArray(:)); end;

        PlMaxValue = max([ PlxMaxValue , PlyMaxValue , PlzMaxValue ]);
        
        if ~isfield(e,'VolEps'), e.VolEps = ones(VolMaxValue,1);  end;
        if ~isfield(e,'VolMu'),  e.VolMu  = ones(VolMaxValue,1);  end;
        if ~isfield(e,'VolSig'), e.VolSig = zeros(VolMaxValue,1); end;
        if ~isfield(e,'VolFaceColor'), e.VolFaceColor = repmat(.5,VolMaxValue,3); end;
        if ~isfield(e,'VolEdgeColor'), e.VolEdgeColor = repmat( 0,VolMaxValue,3); end;
        if ~isfield(e,'VolAlpha'), e.VolAlpha = ones(VolMaxValue,1); end;
        if ~isfield(e,'PlMetal'),  e.PlMetal  = ones(PlMaxValue, 1); end;
        if ~isfield(e,'PlFaceColor'), e.PlFaceColor = repmat(.5,PlMaxValue,3); end;
        if ~isfield(e,'PlEdgeColor'), e.PlEdgeColor = repmat( 0,PlMaxValue,3); end;
        if ~isfield(e,'PlAlpha'), e.PlAlpha = ones(PlMaxValue,1); end;
        if ~isfield(e,'GNDPoint'),  e.GNDPoint  = []; end;
        if ~isfield(e,'Direction'), e.Direction = []; end;
        if ~isfield(e,'Tweaks'), e.Tweaks = struct([]); end;

        %% compatibility issues
        if ~iscell(e.VolEdgeColor),
            e.VolEdgeColor = num2cell(e.VolEdgeColor,2);
        end;
        if ~iscell(e.PlEdgeColor),
            e.PlEdgeColor = num2cell(e.PlEdgeColor,2);
        end;

    else

        e = v.Geometry{GeomPointer};

    end;

    %% main dialog window
    ExtObjDialog = dialog( ...
        'Visible','off', ...
        'Name','External object', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 295 235], ...
        'KeyPressFcn',{@RetEscPressed,@ExtObjStore,@CloseFigure} );

    %% local geometry panel
    ExtObjGeometry = uipanel( ...
        'Parent',ExtObjDialog, ...
        'Units','pixels', ...
        'Position',[10 135 275 95], ...
        'ForegroundColor','b', ...
        'Title','Geometry' );

    %% texts
    uitext(ExtObjGeometry,'x',[ 70 66 60 16],'HorizontalAlignment','center');
    uitext(ExtObjGeometry,'y',[135 66 60 16],'HorizontalAlignment','center');
    uitext(ExtObjGeometry,'z',[200 66 60 16],'HorizontalAlignment','center');
    uitext(ExtObjGeometry,'Position:',[10 45 60 17]);
    uitext(ExtObjGeometry,'Units:'   ,[10 10 40 17]);

    %% uicontrols
    x1 = uiedit(ExtObjGeometry,'',[ 70 45 60 21]);
    y1 = uiedit(ExtObjGeometry,'',[135 45 60 21]);
    z1 = uiedit(ExtObjGeometry,'',[200 45 60 21]);
    UnitsUI = uipopupmenu(ExtObjGeometry,Units,1,[50 10 55 21]);

    data_handles = [x1 y1 z1 UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@ExtObjStore,@CloseFigure});

    %% Comment
    uitext(ExtObjDialog,'Comment:',[10 105 80 16]);
    CommentUI = uimultiedit(ExtObjDialog,'',[10 55 275 50], ...
        'KeyPressFcn',@TabFocus);

    %% buttons
    OKButton = Button(ExtObjDialog,'OK',2,@ExtObjStore);
    CancelButton = Button(ExtObjDialog,'Cancel',1,@CloseFigure);

    %% casting data if editing and not new
    if GeomPointer~=0,
        set(x1,'String',num2str(v.Geometry{GeomPointer}.point1(1)));
        set(y1,'String',num2str(v.Geometry{GeomPointer}.point1(2)));
        set(z1,'String',num2str(v.Geometry{GeomPointer}.point1(3)));
        set(UnitsUI,'Value',v.Geometry{GeomPointer}.units);
        set(CommentUI,'String',v.Geometry{GeomPointer}.comment);
    end;

    %% locking the uicontrols
    if v.locked_flag, LockDialog(ExtObjDialog,OKButton,CancelButton); end;

    uicontrol(x1); %% give focus to the first field
    set(ExtObjDialog,'Visible','on'); %% display figure

    function ExtObjStore(src,eventdata);

        %% get entry point
        e.point1 = [str2double(get(x1,'String')), ...
                    str2double(get(y1,'String')), ...
                    str2double(get(z1,'String'))];

        check = ~isfinite(e.point1) | imag(e.point1)~=0;
        if any(check),
            errordlg('Position coordinates must be real numeric values', ...
                'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% rounding in case that the coordinates are not integers
        e.point1 = round(e.point1);
        
        %% checking the limits of the object and the domain boundaries
        ExtObjSize = zeros(1,3);
        for i = 1:3,
            ExtObjSize(i) = max( [ size(e.VolArray,i) ; ...
                                   size(e.PlxArray,i)-(i==1) ; ...
                                   size(e.PlyArray,i)-(i==2) ; ...
                                   size(e.PlzArray,i)-(i==3) ] );
        end;
        if any(                e.point1<v.limits(1,:) ) | ...
           any( (e.point1 + ExtObjSize)>v.limits(2,:) ),
            h = warndlg({'Object frame extends beyond the domain.'; ...
                         'The object will be clipped.'}, ...
                         'Domain excessed','modal');
            waitfor(h);
        end;
        
        e.units = get(UnitsUI,'Value');
        e.comment = get(CommentUI,'String');

        %% if new record, set the pointer after the last record
        if GeomPointer==0, GeomPointer = length(v.Geometry) + 1; end;

        %% updating (or adding) the geometry record
        v.Geometry{GeomPointer} = e;

        %% adding a space for every newline character
        CommentString = strcat(v.Geometry{GeomPointer}.comment,{' '});
        CommentString = cat(2,CommentString{:}); %% in one line

        %% invisibility flag
        if e.visible, VisFlag = ''; else VisFlag = InvisibilityChar; end;
        
        %% line to include in the geometry listbox
        ListLine = [VisFlag 'extobj' SeparationString CommentString];

        %% including line
        v.GeometryString = strvcat(v.GeometryString(1:GeomPointer-1,:), ...
                                   ListLine, ...
                                   v.GeometryString(GeomPointer+1:end,:));

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        %% display
        axes(GeomAxes);
        if ~isempty(e.handle), delete(e.handle); end;
        DisplayExtObj(GeomPointer);

        close(ExtObjDialog);

    end %% end ExtObjStore

end %% end extobj



%%%%%%%%%%%%%%%%%%%%%%
%% DUPLICATE DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%
function Duplicate(src,eventdata,GeomPointer);
        
    %% main dialog window
    DuplicateDialog = dialog( ...
        'Visible','off', ...
        'Name','Duplicate', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 285 140], ...
        'KeyPressFcn',{@RetEscPressed,@DuplicateStore,@CloseFigure} );

    %% local geometry panel
    DuplicateGeometry = uipanel( ...
        'Parent',DuplicateDialog, ...
        'Units','pixels', ...
        'Position',[10 45 265 90], ...
        'ForegroundColor','b', ...
        'Title','Geometry shift' );

    %% texts
    uitext(DuplicateGeometry,'x',[ 60 61 60 16],'HorizontalAlignment','center');
    uitext(DuplicateGeometry,'y',[125 61 60 16],'HorizontalAlignment','center');
    uitext(DuplicateGeometry,'z',[190 61 60 16],'HorizontalAlignment','center');
    uitext(DuplicateGeometry,'Shift:',[10 40 50 17]);
    uitext(DuplicateGeometry,'Units:'   ,[10 10 40 17]);

    %% uicontrols
    x1 = uiedit(DuplicateGeometry,'',[ 60 40 60 21]);
    y1 = uiedit(DuplicateGeometry,'',[125 40 60 21]);
    z1 = uiedit(DuplicateGeometry,'',[190 40 60 21]);
    UnitsUI = uipopupmenu(DuplicateGeometry,Units,1,[50 10 55 21]);

    data_handles = [x1 y1 z1 UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@DuplicateStore,@CloseFigure});

    %% buttons
    OKButton = Button(DuplicateDialog,'OK',2,@DuplicateStore);
    CancelButton = Button(DuplicateDialog,'Cancel',1,@CloseFigure);

    uicontrol(x1); %% give focus to the first field
    set(DuplicateDialog,'Visible','on'); %% display figure

    function DuplicateStore(src,eventdata);

        %% get shift
        shift = [str2double(get(x1,'String')), ...
                 str2double(get(y1,'String')), ...
                 str2double(get(z1,'String'))];

        check = ~isfinite(shift) | imag(shift)~=0;
        if any(check),
            errordlg('Shift coordinates must be real numeric values', ...
                     'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% adding lines
        v.GeometryString = strvcat(v.GeometryString, ...
                                   v.GeometryString(GeomPointer,:));

        %% adding duplicated elements to the end
        v.Geometry = { v.Geometry{:} v.Geometry{GeomPointer} };

        %% deselecting the old objects
        if Highlight, set(findobj(GeomAxes),'Selected','off'); end;
        
        %% pointing at the duplicated elements
        GeomPointer = length(v.Geometry)-(length(GeomPointer)-1:-1:0);
        
        %% translating and displaying
        axes(GeomAxes);
        for i = GeomPointer,
            Type = v.Geometry{i}.type;
            
            switch Type
                case {'wire','brick','capacitor','resistor','field source'}
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
                    v.Geometry{i}.point2 = v.Geometry{i}.point2 + shift;
                case {'triangle','rectangle'}
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
                    v.Geometry{i}.point2 = v.Geometry{i}.point2 + shift;
                    v.Geometry{i}.point3 = v.Geometry{i}.point3 + shift;
                case 'extobj'
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
                case 'source'
                    v.Geometry{i}.gndpoint = v.Geometry{i}.gndpoint + shift;
                    v.Geometry{i}.point2 = v.Geometry{i}.point2 + shift;
                case 'probe'
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
            end;
            
            switch Type
                case 'wire',         DisplayWire(i,'k');
                case 'triangle',     DisplayPlate(i);
                case 'rectangle',    DisplayPlate(i);
                case 'brick',        DisplayBrick(i);
                case 'capacitor',    DisplayWire(i,'g');
                case 'resistor',     DisplayWire(i,'b');
                case 'source',       DisplaySource(i);
                case 'field source', DisplayFieldSource(i);
                case 'probe',        DisplayProbe(i);
                case 'extobj',       DisplayExtObj(i);
            end;
            
            %% selecting object in Geometry window
            if Highlight, set(v.Geometry{i}.handle,'Selected','on'); end;
            
        end;

        %% string and highlighted item update
        set(GeomListBox,'String',v.GeometryString, ...
                        'Value',GeomPointer);

        close(DuplicateDialog);

    end %% end DuplicateStore

end %% end Duplicate



%%%%%%%%%%%%%%%%%%%%%%
%% TRANSLATE DIALOG %%
%%%%%%%%%%%%%%%%%%%%%%
function Translate(src,eventdata,GeomPointer);
        
    %% main dialog window
    TranslateDialog = dialog( ...
        'Visible','off', ...
        'Name','Translate', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[155 135] 285 140], ...
        'KeyPressFcn',{@RetEscPressed,@TranslateStore,@CloseFigure} );

    %% local geometry panel
    TranslateGeometry = uipanel( ...
        'Parent',TranslateDialog, ...
        'Units','pixels', ...
        'Position',[10 45 265 90], ...
        'ForegroundColor','b', ...
        'Title','Geometry shift' );

    %% texts
    uitext(TranslateGeometry,'x',[ 60 61 60 16],'HorizontalAlignment','center');
    uitext(TranslateGeometry,'y',[125 61 60 16],'HorizontalAlignment','center');
    uitext(TranslateGeometry,'z',[190 61 60 16],'HorizontalAlignment','center');
    uitext(TranslateGeometry,'Shift:',[10 40 50 17]);
    uitext(TranslateGeometry,'Units:'   ,[10 10 40 17]);

    %% uicontrols
    x1 = uiedit(TranslateGeometry,'',[ 60 40 60 21]);
    y1 = uiedit(TranslateGeometry,'',[125 40 60 21]);
    z1 = uiedit(TranslateGeometry,'',[190 40 60 21]);
    UnitsUI = uipopupmenu(TranslateGeometry,Units,1,[50 10 55 21]);

    data_handles = [x1 y1 z1 UnitsUI];

    %% default behavior of Enter and Escape on controls
    set(data_handles,'KeyPressFcn',{@RetEscPressed,@TranslateStore,@CloseFigure});

    %% buttons
    OKButton = Button(TranslateDialog,'OK',2,@TranslateStore);
    CancelButton = Button(TranslateDialog,'Cancel',1,@CloseFigure);

    uicontrol(x1); %% give focus to the first field
    set(TranslateDialog,'Visible','on'); %% display figure

    function TranslateStore(src,eventdata);

        %% get shift
        shift = [str2double(get(x1,'String')), ...
                 str2double(get(y1,'String')), ...
                 str2double(get(z1,'String'))];

        check = ~isfinite(shift) | imag(shift)~=0;
        if any(check),
            errordlg('Shift coordinates must be real numeric values', ...
                'Bad input','modal');
            pause(1); %% otherwise the focus is not given to the handle!
            uicontrol(data_handles(find(check,1)));
            return;
        end;

        %% translating and displaying
        axes(GeomAxes);
        for i = GeomPointer,
            Type = v.Geometry{i}.type;
            if ~isempty(v.Geometry{i}.handle),
                delete(v.Geometry{i}.handle);
            end;
            
            switch Type
                case {'wire','brick','capacitor','resistor','field source'}
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
                    v.Geometry{i}.point2 = v.Geometry{i}.point2 + shift;
                case {'triangle','rectangle'}
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
                    v.Geometry{i}.point2 = v.Geometry{i}.point2 + shift;
                    v.Geometry{i}.point3 = v.Geometry{i}.point3 + shift;
                case 'extobj'
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
                case 'source'
                    v.Geometry{i}.gndpoint = v.Geometry{i}.gndpoint + shift;
                    v.Geometry{i}.point2 = v.Geometry{i}.point2 + shift;
                case 'probe'
                    v.Geometry{i}.point1 = v.Geometry{i}.point1 + shift;
            end;
            
            switch Type
                case 'wire',         DisplayWire(i,'k');
                case 'triangle',     DisplayPlate(i);
                case 'rectangle',    DisplayPlate(i);
                case 'brick',        DisplayBrick(i);
                case 'capacitor',    DisplayWire(i,'g');
                case 'resistor',     DisplayWire(i,'b');
                case 'source',       DisplaySource(i);
                case 'field source', DisplayFieldSource(i);
                case 'probe',        DisplayProbe(i);
                case 'extobj',       DisplayExtObj(i);
            end;
            
            %% selecting object in Geometry window
            if Highlight, set(v.Geometry{i}.handle,'Selected','on'); end;

        end;

        close(TranslateDialog);

    end %% end TranslateStore

end %% end Translate



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%% SIMULATION %%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Simulate(src,eventdata,arg,DeployDir)

    %% checking the stability limit
    if v.Dt > (v.Dx/(v.c*sqrt(3))),
        Button = questdlg( ...
            {'Stability condition violated.', ...
             '', ...
             'Proceed?'}, ...
            'Warning','No');
        if ~strcmp(Button,'Yes'), return; end;
    end;

    %% checking the FDTD dispersion limit in free space
    if strcmp(v.ExcitationPulse,'gauss_sine') && ...
       strcmp(v.GaussSineShape,'auto') && ...
       v.fmax > (v.c/(10*v.Dx)),
        Button = questdlg( ...
            {'Maximum frequency of the excitation pulse results in', ...
             'less than ten cells per wavelength in free space.', ...
             '', ...
             'Proceed?'}, ...
            'Warning','No');
        if ~strcmp(Button,'Yes'), return; end;
    end;

    %% flags
    pml_flag    = v.pml;
    nfff_flag   = v.nfff_flag;
    field_flag  = v.field_flag;
    energy_flag = v.energy_flag;
    drop_flag   = v.drop_flag;
    drop2_flag  = v.drop2_flag;
    
    %% translation to 'old' 'non-structured' variables
    f    = v.fc;
    Dx   = v.Dx;
    Dt   = v.Dt;
    Df   = v.Df;
    Nt   = v.Nt;
    eps0 = v.eps0;
    mu0  = v.mu0;
    c    = v.c;

    %% origin of the FDTD coordinate system
    origin = v.limits(1,:);

    %% lattice dimensions (excluding PML)
    xmax = abs( v.limits(2,1) - v.limits(1,1) );
    ymax = abs( v.limits(2,2) - v.limits(1,2) );
    zmax = abs( v.limits(2,3) - v.limits(1,3) );

    %% PML parameters
    m_s       = v.m_s;
    m_k       = v.m_k;
    m_a       = v.m_a;
    sig_max   = v.sig_max;
    kappa_max = v.kappa_max;
    a_max     = v.a_max;
    
    pmax   = v.pmax;   if ~pml_flag,    pmax   = 0; end;
    e_step = v.e_step; if ~energy_flag, e_step = 0; end;

    pmax2 = 2*pmax;
    
    %% lattice dimensions (including PML)
    xmax_pml = xmax + pmax2;
    ymax_pml = ymax + pmax2;
    zmax_pml = zmax + pmax2;

    %% near to far field transform
    freq_nf = v.freq_nf; %% frequencies of near fields
    freq_ff = v.freq_ff; %% frequencies of radiation pattern
    freq_nf_size = length(freq_nf);
    freq_ff_size = length(freq_ff);
    
    x_nfff = v.x_nfff; %% x-distance from PML of the NFFF int. surface [(D)]
    y_nfff = v.y_nfff; %% y-distance from PML of the NFFF int. surface [(D)]
    z_nfff = v.z_nfff; %% z-distance from PML of the NFFF int. surface [(D)]
    phase_center = v.phcenter.*v.Dx; %% phase center [m]
    resolution = v.resolution; %% pattern resolution in theta and phi
    
    %% phase center with respect to the FDTD domain
    phase_center = phase_center - origin.*Dx;

    %% integration surface -- real coordinates
    x_nfff(2) = xmax - x_nfff(2);
    y_nfff(2) = ymax - y_nfff(2);
    z_nfff(2) = zmax - z_nfff(2);

    %% refresh the excitation pulse
    if isequal(v.ExcitationPulse,'gauss_sine'),
        switch v.GaussSineShape
            case 'auto'
                [fc,fs,offset] = ComputeGaussSinePars( ...
                    v.fmin,v.fmax,v.td_thres,v.fd_thres);
                v.Vs = GaussSinePulse(v.mag,fc,fs,offset);
                if drop2_flag, source_delay = round(10/fc/Dt); end;
            case 'manual'
                v.Vs = GaussSinePulse(v.mag,v.fc,v.fs,v.offset);
                if drop2_flag, source_delay = round(10/v.fc/Dt); end;
                fc = v.fc;
                fs = v.fs;
                offset = v.offset;
        end;
    end;

    %% checking the number of time steps
    if Nt < length(v.Vs),
        h = errordlg({['Number of time steps is lower than ', ...
                       'the length of the excitation pulse.'], ...
                      'Aborting simulation.'},'Error','modal');
        uiwait;
        return;
    end;

    %% adjusting frequency resolution if it is too coarse
    if (v.Dt*v.Df*v.Nt)>1,
        h = errordlg('Frequency resolution too coarse - adjusting', ...
                     'Frequency resolution','modal');
        uiwait;
        Df_orig = v.Df;
        Df_factor = 1;
        Df_scaling = [2 2.5 2]; % stepping by 1, 2, 5, 10, 20, 50, ...
        i = 0;
        while (v.Dt*v.Df*v.Nt)>1,
            Df_factor = round(Df_factor*Df_scaling(mod(i,3)+1));
            v.Df = Df_orig/Df_factor;
            i = i + 1;
        end;
        set(Df_Edit,'String',num2str(v.Df,'%.4g'));
    end;

    %% wave impedance and wave number
    eta = sqrt(mu0/eps0);
    alpha = (2*pi).*freq_ff./c;
    
    %% delay for the termination condition should be the excitation pulse
    %% length
    delay = length(v.Vs);

    %% number of probes
    probe_count = 0;
    for i = 1:length(v.Geometry),
        if strcmp(v.Geometry{i}.type,'probe'),
            probe_count = probe_count + 1;
        end;
    end;
    probe_flag = probe_count > 0;
    
    if strcmp(arg,'fortran') && pref.NewInOutFiles,

        in = v.Vs;
        %% padding the excitation pulse with zeros
        if Nt > length(in), in(Nt) = 0; end;

        FortranDeploy;
        return;
        
    end;

    %% source type
    SourceType = [];
    
    %% lumped source variables
    SourceSubs = [];
    SourceDir = [];

    %% field source variables
    NoOfFieldSources = 0;
    FieldSources = {};
    
    %% relative permittivity tensor
    eps_rx  = ones(xmax_pml,ymax_pml,zmax_pml);
    eps_ry  = ones(xmax_pml,ymax_pml,zmax_pml);
    eps_rz  = ones(xmax_pml,ymax_pml,zmax_pml);

    %% conductivity tensor
    sig_x  = zeros(xmax_pml,ymax_pml,zmax_pml);
    sig_y  = zeros(xmax_pml,ymax_pml,zmax_pml);
    sig_z  = zeros(xmax_pml,ymax_pml,zmax_pml);

    %% relative permeability tensor
    mu_rx  = ones(xmax_pml,ymax_pml,zmax_pml);
    mu_ry  = ones(xmax_pml,ymax_pml,zmax_pml);
    mu_rz  = ones(xmax_pml,ymax_pml,zmax_pml);

    %% in the following, the z-arrays serve as voxels arrays
    %% (until the averaging)
    
    for GeomPointer = 1:length(v.Geometry),
        Type = v.Geometry{GeomPointer}.type;
        switch Type

            %% EXTERNAL OBJECT - volume arrays this time
            case 'extobj'

                VolArray = v.Geometry{GeomPointer}.VolArray;
                VolEps   = v.Geometry{GeomPointer}.VolEps;
                VolMu    = v.Geometry{GeomPointer}.VolMu;
                VolSig   = v.Geometry{GeomPointer}.VolSig;
                point1   = v.Geometry{GeomPointer}.point1;

                %% if volume does not exist, next iteration
                if isempty(VolArray), continue; end;

                %% cropping the array to the FDTD domain
                VolArraySize = size(VolArray);
                vol_limits = [ point1 ; point1+VolArraySize ];
                crop = zeros(2,3);

                crop(1,:) = v.limits(1,:) - vol_limits(1,:);
                crop(2,:) = vol_limits(2,:) - v.limits(2,:);
                crop(crop<0) = 0;

                VolArray = VolArray( 1+crop(1,1):end-crop(2,1), ...
                                     1+crop(1,2):end-crop(2,2), ...
                                     1+crop(1,3):end-crop(2,3) );
                point1 = point1 + crop(1,:) - origin;
                VolArraySize = size(VolArray);

                %% test for overlapping objects
                mask = false(size(eps_rz));
                mask(1+point1(1):point1(1)+VolArraySize(1), ...
                     1+point1(2):point1(2)+VolArraySize(2), ...
                     1+point1(3):point1(3)+VolArraySize(3)) = true;

                if any( (eps_rz(mask)~=1 | mu_rz(mask)~=1 | sig_z(mask)~=0 ) & VolArray(:) ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;

                %% filling the arrays
                mask(mask) = VolArray~=0;
                eps_rz(mask) = VolEps(VolArray(VolArray~=0));
                sig_z(mask) = VolSig(VolArray(VolArray~=0));
                mu_rz(mask) = VolMu(VolArray(VolArray~=0));
                
                
            %% BRICK
            case 'brick'
                
                point1 = round( v.Geometry{GeomPointer}.point1 );
                point2 = round( v.Geometry{GeomPointer}.point2 );

                if v.Geometry{GeomPointer}.PEC,
                    permittivity = Inf;
                    conductivity = 0;
                else
                    permittivity = v.Geometry{GeomPointer}.permittivity;
                    conductivity = v.Geometry{GeomPointer}.conductivity;
                end;

                %% brick extent, cropped
                if pref.CheckBrickLimits % strict checking
                    brick_limits = [ max( min(point1,point2), v.limits(1,:) ) ; ... 
                                     min( max(point1,point2), v.limits(2,:) ) ];
                else % bricks allowed in PML
                    brick_limits = [ max( min(point1,point2), v.limits(1,:) ) ; ... 
                                     min( max(point1,point2), v.limits(2,:)+2*pmax ) ];
                end;
                brick_limits = brick_limits - origin([1 1],:);
                             
                %% test for overlapping objects
                mask = false(size(eps_rz));
                mask( brick_limits(1,1)+1:brick_limits(2,1), ...
                      brick_limits(1,2)+1:brick_limits(2,2), ...
                      brick_limits(1,3)+1:brick_limits(2,3) ) = true;
                if any( eps_rz(mask)~=1 | mu_rz(mask)~=1 | sig_z(mask)~=0 ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;
                
                %% filling the arrays
                eps_rz(mask) = permittivity;
                sig_z(mask)  = conductivity;
                
                
        end; %% end element switch
    end; %% end for

    
    %% checking the FDTD dispersion limit in materials
    %% finding maximum of eps_r*mu_r : MAX does not work, as it counts Inf
    size_em = size(eps_rz);
    max_em  = 1;
    max_ind = 1;
    for i = 1:prod(size_em),
        em = eps_rz(i)*mu_rz(i);
        if isfinite(em),
            if em > max_em,
                max_em = em;
                max_ind = i;
            end;
        end;
    end;
    c_min = c/sqrt(max_em);
    
    if strcmp(v.ExcitationPulse,'gauss_sine') && ...
       strcmp(v.GaussSineShape,'auto') && ...
       v.fmax > (c_min/(10*v.Dx)),
        [max_x,max_y,max_z] = ind2sub(size_em,max_ind);
        Button = questdlg( ...
            {'Maximum frequency of the excitation pulse results in', ...
             'less than ten cells per wavelength in the solid materials.', ...
             ['(Approximate coordinates: ',mat2str(origin + [max_x,max_y,max_z]),')'], ...
             '', ...
             'Proceed?'}, ...
            'Warning','No');
        if ~strcmp(Button,'Yes'),
            return;
        end;
    end;


    %% averaging - z-arrays are voxel arrays
    eps_rx = ( eps_rz + ...
               eps_rz(:,[end 1:end-1],:) + ...
               eps_rz(:,:,[end 1:end-1]) + ...
               eps_rz(:,[end 1:end-1],[end 1:end-1]) )./4;
    eps_ry = ( eps_rz + ...
               eps_rz(:,:,[end 1:end-1]) + ...
               eps_rz([end 1:end-1],:,:) + ...
               eps_rz([end 1:end-1],:,[end 1:end-1]) )./4;
    eps_rz = ( eps_rz + ...
               eps_rz([end 1:end-1],:,:) + ...
               eps_rz(:,[end 1:end-1],:) + ...
               eps_rz([end 1:end-1],[end 1:end-1],:) )./4;

    sig_x = ( sig_z + ...
              sig_z(:,[end 1:end-1],:) + ...
              sig_z(:,:,[end 1:end-1]) + ...
              sig_z(:,[end 1:end-1],[end 1:end-1]) )./4;
    sig_y = ( sig_z + ...
              sig_z(:,:,[end 1:end-1]) + ...
              sig_z([end 1:end-1],:,:) + ...
              sig_z([end 1:end-1],:,[end 1:end-1]) )./4;
    sig_z = ( sig_z + ...
              sig_z([end 1:end-1],:,:) + ...
              sig_z(:,[end 1:end-1],:) + ...
              sig_z([end 1:end-1],[end 1:end-1],:) )./4;

    if v.MuAveraging,
        %% original mu average
        mu_rx = 2./( 1./mu_rz + 1./mu_rz([end 1:end-1],:,:) );
        mu_ry = 2./( 1./mu_rz + 1./mu_rz(:,[end 1:end-1],:) );
        mu_rz = 2./( 1./mu_rz + 1./mu_rz(:,:,[end 1:end-1]) );
    else
        %% no mu average
        mu_rx = max( mu_rz, mu_rz([end 1:end-1],:,:) );
        mu_ry = max( mu_rz, mu_rz(:,[end 1:end-1],:) );
        mu_rz = max( mu_rz, mu_rz(:,:,[end 1:end-1]) );
    end;
    
    
    
    for GeomPointer = 1:length(v.Geometry),
        Type = v.Geometry{GeomPointer}.type;
        switch Type

            %% RESISTOR
            case 'resistor'

                %% capacitor endpoints rounded to grid
                point1 = round( v.Geometry{GeomPointer}.point1 );
                point2 = round( v.Geometry{GeomPointer}.point2 );
                
                %% number of segments is a 1-norm of vector
                NumberOfSegments = sum( abs( point1 - point2 ) );
                
                %% equivalent conductivity of the grid segment
                sig_R = NumberOfSegments / v.Geometry{GeomPointer}.resistance / Dx;
                
                handle = v.Geometry{GeomPointer}.handle;

                %% retrieving vertices from the axes
                vertices = [ get(handle,'XData').' , ...
                             get(handle,'YData').' , ...
                             get(handle,'ZData').' ];

                %% shifting vertices into FDTD grid coordinates
                vertices = vertices - origin(ones(size(vertices,1),1),:);

                %% lengths of wire straight segments
                diff_vertices = diff(vertices);

                for i = 1:size(vertices,1)-1,
                    switch find(diff_vertices(i,:)) %% which direction
                        case 1, %% x-dir
                            sub = vertices(i,1):sign(diff_vertices(i,1)):vertices(i+1,1);
                            sub = sort(sub);
                            sub(1) = [];
                            sig_x(sub,vertices(i,2)+1,vertices(i,3)+1) = sig_R;
                        case 2, %% y-dir
                            sub = vertices(i,2):sign(diff_vertices(i,2)):vertices(i+1,2);
                            sub = sort(sub);
                            sub(1) = [];
                            sig_y(vertices(i,1)+1,sub,vertices(i,3)+1) = sig_R;
                        case 3, %% z-dir
                            sub = vertices(i,3):sign(diff_vertices(i,3)):vertices(i+1,3);
                            sub = sort(sub);
                            sub(1) = [];
                            sig_z(vertices(i,1)+1,vertices(i,2)+1,sub) = sig_R;
                    end;
                end;
                %% Recistor subscripts
                SourceSubs2 = v.Geometry{GeomPointer}.point1 - origin + [1 1 1] - ...
                    ( v.Geometry{GeomPointer}.point2 );

            %% CAPACITOR
            case 'capacitor'

                %% capacitor endpoints rounded to grid
                point1 = round( v.Geometry{GeomPointer}.point1 );
                point2 = round( v.Geometry{GeomPointer}.point2 );
                
                %% number of segments is a 1-norm of vector
                NumberOfSegments = sum( abs( point1 - point2 ) );
                
                %% equivalent relative permittivity of the grid segment
                eps_r_C = 1 + ...
                    v.Geometry{GeomPointer}.capacitance * 1e-9 * NumberOfSegments ...
                    / eps0 / Dx;
                
                handle = v.Geometry{GeomPointer}.handle;

                %% retrieving vertices from the axes
                vertices = [ get(handle,'XData').' , ...
                             get(handle,'YData').' , ...
                             get(handle,'ZData').' ];

                %% shifting vertices into FDTD grid coordinates
                vertices = vertices - origin(ones(size(vertices,1),1),:);

                %% lengths of wire straight segments
                diff_vertices = diff(vertices);

                for i = 1:size(vertices,1)-1,
                    switch find(diff_vertices(i,:)) %% which direction
                        case 1, %% x-dir
                            sub = vertices(i,1):sign(diff_vertices(i,1)):vertices(i+1,1);
                            sub = sort(sub);
                            sub(1) = [];
                            eps_rx(sub,vertices(i,2)+1,vertices(i,3)+1) = eps_r_C;
                        case 2, %% y-dir
                            sub = vertices(i,2):sign(diff_vertices(i,2)):vertices(i+1,2);
                            sub = sort(sub);
                            sub(1) = [];
                            eps_ry(vertices(i,1)+1,sub,vertices(i,3)+1) = eps_r_C;
                        case 3, %% z-dir
                            sub = vertices(i,3):sign(diff_vertices(i,3)):vertices(i+1,3);
                            sub = sort(sub);
                            sub(1) = [];
                            eps_rz(vertices(i,1)+1,vertices(i,2)+1,sub) = eps_r_C;
                    end;
                end;


        end; %% end element switch

    end; %% end for

    

    %% PEC plates
    pl_x = false(xmax_pml,ymax_pml,zmax_pml);
    pl_y = false(xmax_pml,ymax_pml,zmax_pml);
    pl_z = false(xmax_pml,ymax_pml,zmax_pml);

    for GeomPointer = 1:length(v.Geometry),
        Type = v.Geometry{GeomPointer}.type;
        switch Type

            %% WIRE
            case 'wire'

                handle = v.Geometry{GeomPointer}.handle;

                %% retrieving vertices from the axes
                vertices = [ get(handle,'XData').' , ...
                             get(handle,'YData').' , ...
                             get(handle,'ZData').' ];

                %% shifting vertices into FDTD grid coordinates
                vertices = vertices - origin(ones(size(vertices,1),1),:);

                %% lengths of wire straight segments
                diff_vertices = diff(vertices);

                for i = 1:size(vertices,1)-1,
                    switch find(diff_vertices(i,:)) %% which direction
                        case 1, %% x-dir
                            sub = vertices(i,1):sign(diff_vertices(i,1)):vertices(i+1,1);
                            sub = sort(sub);
                            sub(1) = [];
                            eps_rx(sub,vertices(i,2)+1,vertices(i,3)+1) = Inf;
                        case 2, %% y-dir
                            sub = vertices(i,2):sign(diff_vertices(i,2)):vertices(i+1,2);
                            sub = sort(sub);
                            sub(1) = [];
                            eps_ry(vertices(i,1)+1,sub,vertices(i,3)+1) = Inf;
                        case 3, %% z-dir
                            sub = vertices(i,3):sign(diff_vertices(i,3)):vertices(i+1,3);
                            sub = sort(sub);
                            sub(1) = [];
                            eps_rz(vertices(i,1)+1,vertices(i,2)+1,sub) = Inf;
                    end;
                end;
                
                
            %% TRIANGLE and RECTANGLE
            case {'triangle','rectangle'}
                
                handle = v.Geometry{GeomPointer}.handle;

                %% retrieving the geometrical data from the axes
                faces = get(handle,'Faces');
                vertices = get(handle,'Vertices');
                patch_count = get(handle,'UserData');
                
                %% shifting vertices into FDTD grid coordinates
                vertices(:,1) = vertices(:,1) - origin(1) + 1;
                vertices(:,2) = vertices(:,2) - origin(2) + 1;
                vertices(:,3) = vertices(:,3) - origin(3) + 1;
                
                vertices_ind = sub2ind([xmax_pml ymax_pml zmax_pml], ...
                    vertices(:,1),vertices(:,2),vertices(:,3));

                %% selecting the first vertex of each face
                vertices_ind = vertices_ind(faces(:,1),:);
                
                %% setting corresponding PEC plates
                patch_count = cumsum(patch_count);
                pl_x(vertices_ind(               1:patch_count(1),1)) = 1;
                pl_y(vertices_ind(patch_count(1)+1:patch_count(2),1)) = 1;
                pl_z(vertices_ind(patch_count(2)+1:patch_count(3),1)) = 1;


            %% SOURCE
            case 'source'

                %% source conductivity
                Rs = v.Geometry{GeomPointer}.resistance;
                hardsource_flag = Rs == 0;
                if hardsource_flag, sig_s = 0; else sig_s = 1/(Rs*Dx); end;

                %% source subscripts
                SourceSubs = v.Geometry{GeomPointer}.gndpoint - origin + [1 1 1] - ...
                    ( v.Geometry{GeomPointer}.direction == [2 4 6] );

                %% source direction
                SourceDir = round( v.Geometry{GeomPointer}.direction/2 );

                switch SourceDir
                    case 1,
                        sig_x(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = sig_s;
                    case 2,
                        sig_y(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = sig_s;
                    case 3,
                        sig_z(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = sig_s;
                end;
                
                %% sets source type to "lumped" only if field source is not
                %% already defined
                if isempty(SourceType), SourceType = 'lumped'; end;


            %% EXTERNAL OBJECT - plate arrays this time
            case 'extobj'

                PlxArray = v.Geometry{GeomPointer}.PlxArray;
                PlyArray = v.Geometry{GeomPointer}.PlyArray;
                PlzArray = v.Geometry{GeomPointer}.PlzArray;
                PlMetal  = v.Geometry{GeomPointer}.PlMetal;
                GNDPoint = v.Geometry{GeomPointer}.GNDPoint;
                Direction = v.Geometry{GeomPointer}.Direction;
                point1 = v.Geometry{GeomPointer}.point1;

                %% cropping
                VolArraySize = size(v.Geometry{GeomPointer}.VolArray);
                vol_limits = [ point1 ; point1+VolArraySize ];
                crop = zeros(2,3);

                crop(1,:) = v.limits(1,:) - vol_limits(1,:);
                crop(2,:) = vol_limits(2,:) - v.limits(2,:);
                crop(crop<0) = 0;

                point1 = point1 + crop(1,:) - origin;

                if ~isempty(PlxArray),
                    PlxArray = PlxArray( 1+crop(1,1):end-crop(2,1), ...
                                         1+crop(1,2):end-crop(2,2), ...
                                         1+crop(1,3):end-crop(2,3) );
                    mask = logical(zeros(size(pl_x)));
                    mask(1+point1(1):point1(1)+VolArraySize(1)+1, ...
                         1+point1(2):point1(2)+VolArraySize(2), ...
                         1+point1(3):point1(3)+VolArraySize(3)) = 1;
                    mask(mask) = PlxArray~=0;
                    pl_x(mask) = pl_x(mask) | PlMetal(PlxArray(PlxArray~=0));
                end;
                if ~isempty(PlyArray),
                    PlyArray = PlyArray( 1+crop(1,1):end-crop(2,1), ...
                                         1+crop(1,2):end-crop(2,2), ...
                                         1+crop(1,3):end-crop(2,3) );
                    mask = logical(zeros(size(pl_y)));
                    mask(1+point1(1):point1(1)+VolArraySize(1), ...
                         1+point1(2):point1(2)+VolArraySize(2)+1, ...
                         1+point1(3):point1(3)+VolArraySize(3)) = 1;
                    mask(mask) = PlyArray~=0;
                    pl_y(mask) = pl_y(mask) | PlMetal(PlyArray(PlyArray~=0));
                end;
                if ~isempty(PlzArray),
                    PlzArray = PlzArray( 1+crop(1,1):end-crop(2,1), ...
                                         1+crop(1,2):end-crop(2,2), ...
                                         1+crop(1,3):end-crop(2,3) );
                    mask = logical(zeros(size(pl_z)));
                    mask(1+point1(1):point1(1)+VolArraySize(1), ...
                         1+point1(2):point1(2)+VolArraySize(2), ...
                         1+point1(3):point1(3)+VolArraySize(3)+1) = 1;
                    mask(mask) = PlzArray~=0;
                    pl_z(mask) = pl_z(mask) | PlMetal(PlzArray(PlzArray~=0));
                end;

                %% source arrangement
                if ~isempty(GNDPoint),

                    %% source conductivity
                    Rs = DefaultSourceResistance;
                    hardsource_flag = Rs == 0;
                    if hardsource_flag, sig_s = 0; else sig_s = 1/(Rs*Dx); end;

                    %% source subscripts
                    SourceSubs = GNDPoint + point1 + [1 1 1] - ...
                        ( Direction == [2 4 6] );

                    %% source direction
                    SourceDir = round( Direction/2 );

                    switch SourceDir
                        case 1,
                            sig_x(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = sig_s;
                        case 2,
                            sig_y(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = sig_s;
                        case 3,
                            sig_z(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = sig_s;
                    end;
                    
                    %% sets source type to "lumped" only if field source is not
                    %% already defined
                    if isempty(SourceType), SourceType = 'lumped'; end;

                end;


            %% FIELD SOURCE
            case 'field source'

                %% array corners
                corners = round( [ v.Geometry{GeomPointer}.point1 ; ...
                                   v.Geometry{GeomPointer}.point2 ] );

                %% extremes of the coordinates = boundaries of the source array
                minima = min( corners ) - origin + 1;
                maxima = max( corners ) - origin + 1;
                
                %% check if out of bounds
                if any( minima < 1 ) || any( maxima > [xmax,ymax,zmax] ),
                    h = errordlg({'Field source exceeding the computational domain', ...
                                  'Aborting simulation.'},'Error','modal');
                    uiwait;
                    return;
                end;

                %% boundaries in the particular dimensions
                s.xlim = [ minima(1), maxima(1) ];
                s.ylim = [ minima(2), maxima(2) ];
                s.zlim = [ minima(3), maxima(3) ];
                
                phase = angle( v.Geometry{GeomPointer}.data );
                
                %% phase should be always negative
                phase(phase>0) = phase(phase>0) - 2*pi;

                %% other structure elements
                s.mag       = abs( v.Geometry{GeomPointer}.data );
                s.shift     = phase./(2*pi*fc*v.Dt);
                s.component = v.Geometry{GeomPointer}.component;
                s.update    = v.Geometry{GeomPointer}.update;
                
                NoOfFieldSources = NoOfFieldSources + 1;
                
                FieldSources{NoOfFieldSources} = s;
                
                SourceType = 'field';
                
                
        end; %% end element switch

    end; %% end for

    if isempty(SourceType),
        h = errordlg('Source not specified!','Source error','modal');
        waitfor(h);
        return;
    end;

    %% edges of the metallic surfaces are PEC
    PEC_x = pl_y | pl_y(:,:,[end 1:end-1]) | pl_z | pl_z(:,[end 1:end-1],:);
    PEC_y = pl_z | pl_z([end 1:end-1],:,:) | pl_x | pl_x(:,:,[end 1:end-1]);
    PEC_z = pl_x | pl_x(:,[end 1:end-1],:) | pl_y | pl_y([end 1:end-1],:,:);

    eps_rx(PEC_x) = Inf;
    eps_ry(PEC_y) = Inf;
    eps_rz(PEC_z) = Inf;


    %% Tweaks to material matrices
    for GeomPointer = 1:length(v.Geometry),
        if strcmp(v.Geometry{GeomPointer}.type,'extobj'),
            Tweaks = v.Geometry{GeomPointer}.Tweaks;
            point1 = v.Geometry{GeomPointer}.point1;
            if ~isempty(Tweaks),
                if ~isempty(Tweaks.eps_rx),
                    subs = Tweaks.eps_rx(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    eps_rx(ind) = eps_rx(ind).*Tweaks.eps_rx(:,5) + Tweaks.eps_rx(:,4);
                end;
                if ~isempty(Tweaks.eps_ry),
                    subs = Tweaks.eps_ry(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    eps_ry(ind) = eps_ry(ind).*Tweaks.eps_ry(:,5) + Tweaks.eps_ry(:,4);
                end;
                if ~isempty(Tweaks.eps_rz),
                    subs = Tweaks.eps_rz(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    eps_rz(ind) = eps_rz(ind).*Tweaks.eps_rz(:,5) + Tweaks.eps_rz(:,4);
                end;
                if ~isempty(Tweaks.mu_rx),
                    subs = Tweaks.mu_rx(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    mu_rx(ind) = mu_rx(ind).*Tweaks.mu_rx(:,5) + Tweaks.mu_rx(:,4);
                end;
                if ~isempty(Tweaks.mu_ry),
                    subs = Tweaks.mu_ry(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    mu_ry(ind) = mu_ry(ind).*Tweaks.mu_ry(:,5) + Tweaks.mu_ry(:,4);
                end;
                if ~isempty(Tweaks.mu_rz),
                    subs = Tweaks.mu_rz(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    mu_rz(ind) = mu_rz(ind).*Tweaks.mu_rz(:,5) + Tweaks.mu_rz(:,4);
                end;
                if ~isempty(Tweaks.sig_x),
                    subs = Tweaks.sig_x(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    sig_x(ind) = sig_x(ind).*Tweaks.sig_x(:,5) + Tweaks.sig_x(:,4);
                end;
                if ~isempty(Tweaks.sig_y),
                    subs = Tweaks.sig_y(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    sig_y(ind) = sig_y(ind).*Tweaks.sig_y(:,5) + Tweaks.sig_y(:,4);
                end;
                if ~isempty(Tweaks.sig_z),
                    subs = Tweaks.sig_z(:,1:3);
                    subs(:,1) = subs(:,1) + point1(1) - origin(1);
                    subs(:,2) = subs(:,2) + point1(2) - origin(2);
                    subs(:,3) = subs(:,3) + point1(3) - origin(3);
                    ind = sub2ind([xmax_pml ymax_pml zmax_pml],subs(:,1),subs(:,2),subs(:,3));
                    sig_z(ind) = sig_z(ind).*Tweaks.sig_z(:,5) + Tweaks.sig_z(:,4);
                end;
            end;
        end;
    end;

    %% free space coefficients
    Cb0 = Dt./Dx./eps0;
    Db0 = Dt./Dx./mu0;

    %% material dependent coefficients - (3.28), (3.29) in [Taflove]
    s_x = sig_x.*Dt./2./eps0./eps_rx;
    Cax = (1 - s_x)./(1 + s_x);
    Cbx = Cb0./eps_rx./(1 + s_x);
    if strcmp(SourceType,'lumped') && SourceDir==1, %% creating particular input
        in = ( sig_x(SourceSubs(1),SourceSubs(2),SourceSubs(3)) .* ...
            Cbx(SourceSubs(1),SourceSubs(2),SourceSubs(3)) ).*v.Vs;
    end;
    clear eps_rx sig_x s_x;

    s_y = sig_y.*Dt./2./eps0./eps_ry;
    Cay = (1 - s_y)./(1 + s_y);
    Cby = Cb0./eps_ry./(1 + s_y);
    if strcmp(SourceType,'lumped') && SourceDir==2,
        in = ( sig_y(SourceSubs(1),SourceSubs(2),SourceSubs(3)) .* ...
            Cby(SourceSubs(1),SourceSubs(2),SourceSubs(3)) ).*v.Vs;
    end;
    clear eps_ry sig_y s_y;

    s_z = sig_z.*Dt./2./eps0./eps_rz;
    Caz = (1 - s_z)./(1 + s_z);
    Cbz = Cb0./eps_rz./(1 + s_z);
    if strcmp(SourceType,'lumped') && SourceDir==3,
        in = ( sig_z(SourceSubs(1),SourceSubs(2),SourceSubs(3)) .* ...
            Cbz(SourceSubs(1),SourceSubs(2),SourceSubs(3)) ).*v.Vs;
    end;
    clear eps_rz sig_z s_z;
    
    switch SourceType
        
        case 'lumped'
            
            %% hard source case
            if hardsource_flag, in = -v.Vs./Dx; end;
            
        case 'field'
            
            if drop2_flag,
                h = errordlg(['Source energy termination condition ', ...
                    'cannot be used with field source'],'Source error','modal');
                waitfor(h);
                return;
            end;
            
            %% calculate DFT of Vs at fc and divide (normalize) the
            %% magnitudes
            Nt_fft = round( 1/(v.Dt*v.Df) );
            Vs_fft = 2.*fft(v.Vs,Nt_fft)./Nt_fft;
            f_sample = round(fc./v.Df) + 1;
            for i = 1:length(FieldSources),
                FieldSources{i}.mag = FieldSources{i}.mag ./ abs(Vs_fft(f_sample));
            end;

            %% find the highest shift and add to delay for the termination
            %% condition
            max_shift = 0;
            for i = 1:length(FieldSources),
                max_shift = max( max_shift, max( abs(FieldSources{i}.shift(:))) );
            end;
            max_shift = ceil(max_shift);
            delay = delay + max_shift;

            %% offset of the center of the pulse in time steps
            offset_Dt = round( offset/v.Dt );

            %% prepare short indices
            [fx1,fx2,fy1,fy2,fz1,fz2] ...
                = deal( zeros(NoOfFieldSources,1) );
            
            for l = 1:NoOfFieldSources,
                fx1(l) = FieldSources{l}.xlim(1);
                fx2(l) = FieldSources{l}.xlim(2);
                fy1(l) = FieldSources{l}.ylim(1);
                fy2(l) = FieldSources{l}.ylim(2);
                fz1(l) = FieldSources{l}.zlim(1);
                fz2(l) = FieldSources{l}.zlim(2);
            end;
            
    end;
    
    %% definition of the Gauss-Sine pulse
    function out = inG(n_in)
        
        %% span of the Gauss-Sine pulse in time steps
        n = n_in - 1 - offset_Dt;

        %% Gauss-weighted Sine
        if ( n < -offset_Dt ) | ( n > offset_Dt ),
            out = 0;
        else
            out = v.mag .* sin((2*pi*fc*v.Dt).*n) .* exp( -((fs*v.Dt).*n).^2 );
        end;

    end
    
    lumpedsource_flag = strcmp(SourceType,'lumped');
    fieldsource_flag  = strcmp(SourceType,'field');

    Dbx = Db0./mu_rx; clear mu_rx;
    Dby = Db0./mu_ry; clear mu_ry;
    Dbz = Db0./mu_rz; clear mu_rz;

    %% padding the excitation pulse with zeros
    if Nt > length(in), in(Nt) = 0; end;

    %% auxiliary indices
    indx_pml = xmax+1:xmax_pml;
    indy_pml = ymax+1:ymax_pml;
    indz_pml = zmax+1:zmax_pml;

    ix = ones(xmax_pml,1);
    iy = ones(ymax_pml,1);
    iz = ones(zmax_pml,1);


    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% PML initialization %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    if pml_flag,
    
        %% PML conductivity profiles
        sigx1 = sig_max.*([(0  :pmax-1  ) (pmax    :-1:1  )].'./pmax).^m_s;
        sigx2 = sig_max.*([(0.5:pmax-0.5) (pmax-0.5:-1:0.5)].'./pmax).^m_s;

        %% PML kappa profiles
        kappax1 = 1 + ( (kappa_max - 1).*([(0  :pmax-1  ) (pmax    :-1:1  )].'./pmax).^m_k );
        kappax2 = 1 + ( (kappa_max - 1).*([(0.5:pmax-0.5) (pmax-0.5:-1:0.5)].'./pmax).^m_k );

        %% PML a profiles
        ax1 = a_max.*((pmax - [(0  :pmax-1  ) (pmax    :-1:1  )].')./pmax).^m_a;
        ax2 = a_max.*((pmax - [(0.5:pmax-0.5) (pmax-0.5:-1:0.5)].')./pmax).^m_a;

        %% PML coefficients
        bx1 = exp( -(sigx1./kappax1 + ax1).*Dt./eps0 );
        bx2 = exp( -(sigx2./kappax2 + ax2).*Dt./eps0 );

        cx1 = ( bx1 - 1 )./( 1 + kappax1.*ax1.*sigx1.^(-1) );
        cx2 = ( bx2 - 1 )./( 1 + kappax2.*ax2.*sigx2.^(-1) );

        cx1(isnan(cx1)) = 0;
        
        kappay1 = kappax1.';
        kappay2 = kappax2.';

        by1 = bx1.';
        by2 = bx2.';

        cy1 = cx1.';
        cy2 = cx2.';

        kappaz1 = shiftdim(kappay1,-1);
        kappaz2 = shiftdim(kappay2,-1);

        bz1 = shiftdim(by1,-1);
        bz2 = shiftdim(by2,-1);

        cz1 = shiftdim(cy1,-1);
        cz2 = shiftdim(cy2,-1);

        %% changing the Cb terms due to PML
        Cbx(indx_pml,:,:) = Cbx(indx_pml,:,:).*    kappax2(:,iy,iz);
        Cbx(:,indy_pml,:) = Cbx(:,indy_pml,:).*(1./kappay1(ix,:,iz));
        Cbx(:,:,indz_pml) = Cbx(:,:,indz_pml).*(1./kappaz1(ix,iy,:));

        Cby(indx_pml,:,:) = Cby(indx_pml,:,:).*(1./kappax1(:,iy,iz));
        Cby(:,indy_pml,:) = Cby(:,indy_pml,:).*    kappay2(ix,:,iz);
        Cby(:,:,indz_pml) = Cby(:,:,indz_pml).*(1./kappaz1(ix,iy,:));

        Cbz(indx_pml,:,:) = Cbz(indx_pml,:,:).*(1./kappax1(:,iy,iz));
        Cbz(:,indy_pml,:) = Cbz(:,indy_pml,:).*(1./kappay1(ix,:,iz));
        Cbz(:,:,indz_pml) = Cbz(:,:,indz_pml).*    kappaz2(ix,iy,:);

        %% changing the Db terms due to PML
        Dbx(indx_pml,:,:) = Dbx(indx_pml,:,:).*    kappax1(:,iy,iz);
        Dbx(:,indy_pml,:) = Dbx(:,indy_pml,:).*(1./kappay2(ix,:,iz));
        Dbx(:,:,indz_pml) = Dbx(:,:,indz_pml).*(1./kappaz2(ix,iy,:));

        Dby(indx_pml,:,:) = Dby(indx_pml,:,:).*(1./kappax2(:,iy,iz));
        Dby(:,indy_pml,:) = Dby(:,indy_pml,:).*    kappay1(ix,:,iz);
        Dby(:,:,indz_pml) = Dby(:,:,indz_pml).*(1./kappaz2(ix,iy,:));

        Dbz(indx_pml,:,:) = Dbz(indx_pml,:,:).*(1./kappax2(:,iy,iz));
        Dbz(:,indy_pml,:) = Dbz(:,indy_pml,:).*(1./kappay2(ix,:,iz));
        Dbz(:,:,indz_pml) = Dbz(:,:,indz_pml).*    kappaz1(ix,iy,:);

    end;
    
    

    %%%%%%%%%%%%%%%%%%%%
    %% Media indexing %% (optimized for speed)
    %%%%%%%%%%%%%%%%%%%%
    
    Cax = Cax(:);
    Cay = Cay(:);
    Caz = Caz(:);
    Cbx = Cbx(:);
    Cby = Cby(:);
    Cbz = Cbz(:);
    Dbx = Dbx(:);
    Dby = Dby(:);
    Dbz = Dbz(:);


    %% Reducing arrays by excluding repeated materials %%

    rows = size(Cax,1); %% length of array

    d = false(rows-1,1); %% mask where changes in material

    d = d | (Cax(1:rows-1)~=Cax(2:rows));
    d = d | (Cay(1:rows-1)~=Cay(2:rows));
    d = d | (Caz(1:rows-1)~=Caz(2:rows));
    d = d | (Cbx(1:rows-1)~=Cbx(2:rows));
    d = d | (Cby(1:rows-1)~=Cby(2:rows));
    d = d | (Cbz(1:rows-1)~=Cbz(2:rows));
    d = d | (Dbx(1:rows-1)~=Dbx(2:rows));
    d = d | (Dby(1:rows-1)~=Dby(2:rows));
    d = d | (Dbz(1:rows-1)~=Dbz(2:rows));

    d = [true ; d]; %% first line is unique, of course

    Cax = Cax(d,:);
    Cay = Cay(d,:);
    Caz = Caz(d,:);
    Cbx = Cbx(d,:);
    Cby = Cby(d,:);
    Cbz = Cbz(d,:);
    Dbx = Dbx(d,:);
    Dby = Dby(d,:);
    Dbz = Dbz(d,:);

    media = cumsum(d); %% temporary media array


    %% Sorting %%

    rows = size(Cax,1); %% new length of array

    ndx = (1:rows)';

    %% sorting the coefficient from back to front
    [ignore,ind] = sort(Dbz(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Dby(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Dbx(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Cbz(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Cby(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Cbx(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Caz(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Cay(ndx)); ndx = ndx(ind);
    [ignore,ind] = sort(Cax(ndx)); ndx = ndx(ind);

    clear ignore ind;

    %% rearranging in sorted order
    Cax = Cax(ndx,:);
    Cay = Cay(ndx,:);
    Caz = Caz(ndx,:);
    Cbx = Cbx(ndx,:);
    Cby = Cby(ndx,:);
    Cbz = Cbz(ndx,:);
    Dbx = Dbx(ndx,:);
    Dby = Dby(ndx,:);
    Dbz = Dbz(ndx,:);


    %% Creating the lookup table %%

    d = false(rows-1,1); %% mask where changes

    d = d | (Cax(1:rows-1,1)~=Cax(2:rows,1));
    d = d | (Cay(1:rows-1,1)~=Cay(2:rows,1));
    d = d | (Caz(1:rows-1,1)~=Caz(2:rows,1));
    d = d | (Cbx(1:rows-1,1)~=Cbx(2:rows,1));
    d = d | (Cby(1:rows-1,1)~=Cby(2:rows,1));
    d = d | (Cbz(1:rows-1,1)~=Cbz(2:rows,1));
    d = d | (Dbx(1:rows-1,1)~=Dbx(2:rows,1));
    d = d | (Dby(1:rows-1,1)~=Dby(2:rows,1));
    d = d | (Dbz(1:rows-1,1)~=Dbz(2:rows,1));

    d = [true ; d]; %% first line unique

    %% create lookup table by indexing into sorted list.
    Cax = Cax(d,:);
    Cay = Cay(d,:);
    Caz = Caz(d,:);
    Cbx = Cbx(d,:);
    Cby = Cby(d,:);
    Cbz = Cbz(d,:);
    Dbx = Dbx(d,:);
    Dby = Dby(d,:);
    Dbz = Dbz(d,:);

    pos = cumsum(d); %% position mapping vector

    
    %% Creating the media array %%

    pos(ndx) = pos; %% re-reference POS to indexing of SORT.
    media = pos(media); %% real indexing in media array
    clear ndx pos d;

    switch arg
        case 'matlab'
            if size(Cax,1)<2^8,
                media = uint8(media);
            elseif size(Cax,1)<2^16,
                media = uint16(media);
            elseif size(Cax,1)<2^32,
                media = uint32(media);
            elseif size(Cax,1)<2^64,
                media = uint64(media);
            else
                disp('Cannot create media array'); keyboard;
            end;
        case 'fortran' %% fortran needs to pass media in INTEGER*4
            if size(Cax,1)<(2^32/2),
                media = int32(media);
            else
                disp('Cannot create media array'); keyboard;
            end;
    end;

    media = reshape(media,[xmax_pml ymax_pml zmax_pml]);

    
    %% theta and phi steps
    Dth =   pi/resolution(1);
    Dph = 2*pi/resolution(2);

    
    switch arg
        case 'matlab'
            MatlabFDTD;
            ComputeResults; %% compute FFTs and impedance
        case 'fortran'
            FortranDeploy;
    end;



    %% FDTD Matlab kernel
    function MatlabFDTD

        v.time_start = clock;
        tic;
        
        %% array allocation
        Ex = zeros(xmax_pml,ymax_pml,zmax_pml);
        Ey = zeros(xmax_pml,ymax_pml,zmax_pml);
        Ez = zeros(xmax_pml,ymax_pml,zmax_pml);

        Hx = zeros(xmax_pml,ymax_pml,zmax_pml);
        Hy = zeros(xmax_pml,ymax_pml,zmax_pml);
        Hz = zeros(xmax_pml,ymax_pml,zmax_pml);

        %% PML array allocation
        if pml_flag,
            
            Psi_Ey_x = zeros(pmax2,ymax_pml,zmax_pml);
            Psi_Ez_x = zeros(pmax2,ymax_pml,zmax_pml);

            Psi_Ez_y = zeros(xmax_pml,pmax2,zmax_pml);
            Psi_Ex_y = zeros(xmax_pml,pmax2,zmax_pml);

            Psi_Ex_z = zeros(xmax_pml,ymax_pml,pmax2);
            Psi_Ey_z = zeros(xmax_pml,ymax_pml,pmax2);

            Psi_Hy_x = zeros(pmax2,ymax_pml,zmax_pml);
            Psi_Hz_x = zeros(pmax2,ymax_pml,zmax_pml);

            Psi_Hz_y = zeros(xmax_pml,pmax2,zmax_pml);
            Psi_Hx_y = zeros(xmax_pml,pmax2,zmax_pml);

            Psi_Hx_z = zeros(xmax_pml,ymax_pml,pmax2);
            Psi_Hy_z = zeros(xmax_pml,ymax_pml,pmax2);

            %% auxiliary indices
            indx1 = [xmax_pml (1:xmax_pml-1)];
            indy1 = [ymax_pml (1:ymax_pml-1)];
            indz1 = [zmax_pml (1:zmax_pml-1)];

            indx2 = [(2:xmax_pml) 1];
            indy2 = [(2:ymax_pml) 1];
            indz2 = [(2:zmax_pml) 1];

            indx_pml1 = [xmax:xmax_pml-1];
            indy_pml1 = [ymax:ymax_pml-1];
            indz_pml1 = [zmax:zmax_pml-1];

            indx_pml2 = [xmax+2:xmax_pml 1];
            indy_pml2 = [ymax+2:ymax_pml 1];
            indz_pml2 = [zmax+2:zmax_pml 1];
            
            x_wall = xmax + 1 + pmax;
            y_wall = ymax + 1 + pmax;
            z_wall = zmax + 1 + pmax;

        else
            
            indx1 = [1 (1:xmax-1)];
            indy1 = [1 (1:ymax-1)];
            indz1 = [1 (1:zmax-1)];

            indx2 = [(2:xmax) xmax];
            indy2 = [(2:ymax) ymax];
            indz2 = [(2:zmax) zmax];
            
        end;
        
        %% auxiliary indices
        indx = 1:xmax;
        indy = 1:ymax;
        indz = 1:zmax;


        %% NFFF initialization
        if nfff_flag,
        
            nfff_size = [x_nfff(2)-x_nfff(1); ...
                         y_nfff(2)-y_nfff(1); ...
                         z_nfff(2)-z_nfff(1)];
            Wff = exp(-1j.*2.*pi.*((0:Nt-1).'*freq_ff).*Dt);

            %% x-y plane
            Ex_xy_dft  = zeros(nfff_size(1)  ,nfff_size(2)+1,2,freq_ff_size);
            Ey_xy_dft  = zeros(nfff_size(1)+1,nfff_size(2)  ,2,freq_ff_size);
            Hx_xy_dft1 = zeros(nfff_size(1)+1,nfff_size(2)  ,2,freq_ff_size);
            Hx_xy_dft2 = zeros(nfff_size(1)+1,nfff_size(2)  ,2,freq_ff_size);
            Hy_xy_dft1 = zeros(nfff_size(1)  ,nfff_size(2)+1,2,freq_ff_size);
            Hy_xy_dft2 = zeros(nfff_size(1)  ,nfff_size(2)+1,2,freq_ff_size);

            %% y-z plane
            Ey_yz_dft  = zeros(2,nfff_size(2)  ,nfff_size(3)+1,freq_ff_size);
            Ez_yz_dft  = zeros(2,nfff_size(2)+1,nfff_size(3)  ,freq_ff_size);
            Hy_yz_dft1 = zeros(2,nfff_size(2)+1,nfff_size(3)  ,freq_ff_size);
            Hy_yz_dft2 = zeros(2,nfff_size(2)+1,nfff_size(3)  ,freq_ff_size);
            Hz_yz_dft1 = zeros(2,nfff_size(2)  ,nfff_size(3)+1,freq_ff_size);
            Hz_yz_dft2 = zeros(2,nfff_size(2)  ,nfff_size(3)+1,freq_ff_size);

            %% z-x plane
            Ez_zx_dft  = zeros(nfff_size(1)+1,2,nfff_size(3)  ,freq_ff_size);
            Ex_zx_dft  = zeros(nfff_size(1)  ,2,nfff_size(3)+1,freq_ff_size);
            Hz_zx_dft1 = zeros(nfff_size(1)  ,2,nfff_size(3)+1,freq_ff_size);
            Hz_zx_dft2 = zeros(nfff_size(1)  ,2,nfff_size(3)+1,freq_ff_size);
            Hx_zx_dft1 = zeros(nfff_size(1)+1,2,nfff_size(3)  ,freq_ff_size);
            Hx_zx_dft2 = zeros(nfff_size(1)+1,2,nfff_size(3)  ,freq_ff_size);

        end;


        %% Near field initialization
        if field_flag,
            
            v.Ex_field = zeros(xmax,ymax,zmax,freq_nf_size);
            v.Ey_field = zeros(xmax,ymax,zmax,freq_nf_size);
            v.Ez_field = zeros(xmax,ymax,zmax,freq_nf_size);
            
            Wnf = exp(-1j.*2.*pi.*((0:Nt-1).'*freq_nf).*Dt);
            
        end;

        
        %% energy initialization
        if energy_flag,
            Cbx0 = Cbx.*(2./(1 + Cax)); Cbx0(Cbx0==0) = 1;
            Cby0 = Cby.*(2./(1 + Cay)); Cby0(Cby0==0) = 1;
            Cbz0 = Cbz.*(2./(1 + Caz)); Cbz0(Cbz0==0) = 1;
            energy_mult = Dt*(Dx^2)/2; %% multiplier
            e_comp = e_step; %% variable for comparing
            v.energy = zeros(Nt,1);
            energy_now = 0;
            energy_percent = 0;
            energy_max = 0;
            Drop_pct = 100*10^(v.Drop_dB/10);
        else
            energy_now = '#';
            energy_percent = '#';
        end;
        
        
        %% source termination condition initialization
        if drop2_flag,
            source_now = 0;
            source_max = 0;
            source_count = 0;
            source_percent = 0;
            Drop2_pct = 100*10^(v.Drop2_dB/10);
        end;
        
        
        drop_cond = true;
        drop2_cond = true;

        
        %% output (response)
        out = zeros(Nt,1);
 
        
        x_probe = zeros(probe_count,1);
        y_probe = zeros(probe_count,1);
        z_probe = zeros(probe_count,1);
        probe_type = zeros(probe_count,1);
        
        j = 0;
        for i = 1:length(v.Geometry),
            if strcmp(v.Geometry{i}.type,'probe'),
                j = j + 1;
                point1 = v.Geometry{i}.point1 - origin + 1;
                x_probe(j) = point1(1);
                y_probe(j) = point1(2);
                z_probe(j) = point1(3);
                probe_type(j) = v.Geometry{i}.component;
            end;
        end;
        
        if probe_flag,
            v.Probes = zeros(Nt,probe_count);
        end;
        
        
        %% timing
        v.partial_times(1) = toc;
        tic;

        
        %%%%%%%%%%%%%%%%
        %% MAIN CYCLE %%
        %%%%%%%%%%%%%%%%

        for t = 1:Nt,


            %%%%%%%%%%%%%%%%%%%%%%%%
            %% E component update %%
            %%%%%%%%%%%%%%%%%%%%%%%%

            %% regular region - (3.31) in [Taflove]
            for k = 1:zmax_pml,
                iz1 = indz1(k);
                Ex(:,:,k) = Cax(media(:,:,k)).*Ex(:,:,k) ...
                    + Cbx(media(:,:,k)).* ...
                    ( Hz(:,:,k) - Hz(:,indy1,k) ...
                    - Hy(:,:,k) + Hy(:,:,iz1) );
                Ey(:,:,k) = Cay(media(:,:,k)).*Ey(:,:,k) ...
                    + Cby(media(:,:,k)).* ...
                    ( Hx(:,:,k) - Hx(:,:,iz1) ...
                    - Hz(:,:,k) + Hz(indx1,:,k) );
                Ez(:,:,k) = Caz(media(:,:,k)).*Ez(:,:,k) ...
                    + Cbz(media(:,:,k)).* ...
                    ( Hy(:,:,k) - Hy(indx1,:,k) ...
                    - Hx(:,:,k) + Hx(:,indy1,k) );
            end;

            %% PML region
            if pml_flag,

                %% Ex component
                Psi_Ex_y = by1(ix,:,iz).*Psi_Ex_y ...
                    + cy1(ix,:,iz).*( Hz(:,indy_pml,:) - Hz(:,indy_pml1,:) );
                Psi_Ex_z = bz1(ix,iy,:).*Psi_Ex_z ...
                    + cz1(ix,iy,:).*( Hy(:,:,indz_pml) - Hy(:,:,indz_pml1) );

                Ex(:,indy_pml,:) = Ex(:,indy_pml,:) ...
                    + Cbx(media(:,indy_pml,:)).*Psi_Ex_y;
                Ex(:,:,indz_pml) = Ex(:,:,indz_pml) ...
                    - Cbx(media(:,:,indz_pml)).*Psi_Ex_z;

                %% Ey component
                Psi_Ey_z = bz1(ix,iy,:).*Psi_Ey_z ...
                    + cz1(ix,iy,:).*( Hx(:,:,indz_pml) - Hx(:,:,indz_pml1) );
                Psi_Ey_x = bx1(:,iy,iz).*Psi_Ey_x ...
                    + cx1(:,iy,iz).*( Hz(indx_pml,:,:) - Hz(indx_pml1,:,:) );

                Ey(indx_pml,:,:) = Ey(indx_pml,:,:) ...
                    - Cby(media(indx_pml,:,:)).*Psi_Ey_x;
                Ey(:,:,indz_pml) = Ey(:,:,indz_pml) ...
                    + Cby(media(:,:,indz_pml)).*Psi_Ey_z;

                %% Ez component
                Psi_Ez_x = bx1(:,iy,iz).*Psi_Ez_x ...
                    + cx1(:,iy,iz).*( Hy(indx_pml,:,:) - Hy(indx_pml1,:,:) );
                Psi_Ez_y = by1(ix,:,iz).*Psi_Ez_y ...
                    + cy1(ix,:,iz).*( Hx(:,indy_pml,:) - Hx(:,indy_pml1,:) );

                Ez(indx_pml,:,:) = Ez(indx_pml,:,:) ...
                    + Cbz(media(indx_pml,:,:)).*Psi_Ez_x;
                Ez(:,indy_pml,:) = Ez(:,indy_pml,:) ...
                    - Cbz(media(:,indy_pml,:)).*Psi_Ez_y;

                %% PEC walls -- maybe not necessary (or useful)
                Ey(x_wall,:,:) = 0;
                Ez(x_wall,:,:) = 0;
                Ez(:,y_wall,:) = 0;
                Ex(:,y_wall,:) = 0;
                Ex(:,:,z_wall) = 0;
                Ey(:,:,z_wall) = 0;
                
            end;
            
            
            
            %% field source
            if fieldsource_flag,
                
                for l = 1:NoOfFieldSources,
                    
                    switch FieldSources{l}.component
                        
                    %% Ex
                    case 1,

                        %% if replacing, the old value is zeroed
                        if FieldSources{l}.update == 1,
                            Ex(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0;
                        end;

                        Ex(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = ...
                        Ex(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + ...
                        FieldSources{l}.mag .* inG( t + FieldSources{l}.shift );                          

                    %% Ey
                    case 2,

                        %% if replacing, the old value is zeroed
                        if FieldSources{l}.update == 1,
                            Ey(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0;
                        end;

                        Ey(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = ...
                        Ey(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + ...
                        FieldSources{l}.mag .* inG( t + FieldSources{l}.shift );                          

                    %% Ez
                    case 3,

                        %% if replacing, the old value is zeroed
                        if FieldSources{l}.update == 1,
                            Ez(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0;
                        end;

                        Ez(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = ...
                        Ez(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + ...
                        FieldSources{l}.mag .* inG( t + FieldSources{l}.shift );                          
                            
                    end;
                    
                end;
                
            end;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Excitation and readout %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist = 8;
t0 = 14;
            if lumpedsource_flag,
                if hardsource_flag,

                    %% Hard source
                    switch SourceDir
                        case 1, %% x dir
                            out(t) = Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3));
                            Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = in(t);
                            
                             
                            
                            
                        case 2, %% y dir
                            out(t) = Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3));

                            Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = in(t);
                            
%                             Ey(49- origin(1),25- origin(2),57- origin(3)) = in(t+15);
                        case 3, %% z dir
                            out(t) = Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3));
                            Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = in(t);
                            
                    end;

                else

                    %% Resistive source
                    switch SourceDir
                        case 1, %% x dir
                            Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = ...
                                Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)) - in(t);
                            out(t) = Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3));
                        case 2, %% y dir
                            Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = ...
                                Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)) - in(t);
                            out(t) = Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3));
case 3, %% z dir
Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = ...
Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)) - in(t);
 %% ---  added field update --- %% 
Ez(SourceSubs(1)-dist,SourceSubs(2),SourceSubs(3)) = in(t+t0);
 %% --- --- %%
out(t) = Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3));
                    end;

                end;
            end;
            


            %% DFT of NFFF surfaces
            if nfff_flag,
                for p = 1:freq_ff_size,

                    %% x-y plane
                    Ex_xy_dft(:,:,:,p) = Ex_xy_dft(:,:,:,p) ...
                        + Wff(t,p).*Ex( x_nfff(1)+1:x_nfff(2), ...
                        y_nfff(1)+1:y_nfff(2)+1, ...
                        z_nfff + [1 1] );
                    Ey_xy_dft(:,:,:,p) = Ey_xy_dft(:,:,:,p) ...
                        + Wff(t,p).*Ey( x_nfff(1)+1:x_nfff(2)+1, ...
                        y_nfff(1)+1:y_nfff(2), ...
                        z_nfff + [1 1] );

                    %% y-z plane
                    Ey_yz_dft(:,:,:,p) = Ey_yz_dft(:,:,:,p) ...
                        + Wff(t,p).*Ey( x_nfff + [1 1], ...
                        y_nfff(1)+1:y_nfff(2), ...
                        z_nfff(1)+1:z_nfff(2)+1 );
                    Ez_yz_dft(:,:,:,p) = Ez_yz_dft(:,:,:,p) ...
                        + Wff(t,p).*Ez( x_nfff + [1 1], ...
                        y_nfff(1)+1:y_nfff(2)+1, ...
                        z_nfff(1)+1:z_nfff(2) );

                    %% z-x plane
                    Ez_zx_dft(:,:,:,p) = Ez_zx_dft(:,:,:,p) ...
                        + Wff(t,p).*Ez( x_nfff(1)+1:x_nfff(2)+1, ...
                        y_nfff + [1 1], ...
                        z_nfff(1)+1:z_nfff(2) );
                    Ex_zx_dft(:,:,:,p) = Ex_zx_dft(:,:,:,p) ...
                        + Wff(t,p).*Ex( x_nfff(1)+1:x_nfff(2), ...
                        y_nfff + [1 1], ...
                        z_nfff(1)+1:z_nfff(2)+1 );
                    
                end;
            end;

            %% near field DFT
            if field_flag,
                for p = 1:freq_nf_size,

                    %% E-field store
                    v.Ex_field(:,:,:,p) = v.Ex_field(:,:,:,p) + Wnf(t,p).*Ex(indx,indy,indz);
                    v.Ey_field(:,:,:,p) = v.Ey_field(:,:,:,p) + Wnf(t,p).*Ey(indx,indy,indz);
                    v.Ez_field(:,:,:,p) = v.Ez_field(:,:,:,p) + Wnf(t,p).*Ez(indx,indy,indz);

                end;
            end;

            

            %%%%%%%%%%%%%%%%%%%%%%%%
            %% H component update %%
            %%%%%%%%%%%%%%%%%%%%%%%%

            %% regular region - (3.30) in [Taflove] (without magnetic losses)
            for k = 1:zmax_pml,
                iz2 = indz2(k);
                Hx(:,:,k) = Hx(:,:,k) + Dbx(media(:,:,k)).* ...
                    ( Ey(:,:,iz2) - Ey(:,:,k) ...
                    - Ez(:,indy2,k) + Ez(:,:,k) );
                Hy(:,:,k) = Hy(:,:,k) + Dby(media(:,:,k)).* ...
                    ( Ez(indx2,:,k) - Ez(:,:,k) ...
                    - Ex(:,:,iz2) + Ex(:,:,k) );
                Hz(:,:,k) = Hz(:,:,k) + Dbz(media(:,:,k)).* ...
                    ( Ex(:,indy2,k) - Ex(:,:,k) ...
                    - Ey(indx2,:,k) + Ey(:,:,k) );
            end;

            %% PML region
            if pml_flag,

                %% Hx component
                Psi_Hx_y = by2(ix,:,iz).*Psi_Hx_y ...
                    + cy2(ix,:,iz).*( Ez(:,indy_pml2,:) - Ez(:,indy_pml,:) );
                Psi_Hx_z = bz2(ix,iy,:).*Psi_Hx_z ...
                    + cz2(ix,iy,:).*( Ey(:,:,indz_pml2) - Ey(:,:,indz_pml) );

                Hx(:,indy_pml,:) = Hx(:,indy_pml,:) ...
                    - Dbx(media(:,indy_pml,:)).*Psi_Hx_y;
                Hx(:,:,indz_pml) = Hx(:,:,indz_pml) ...
                    + Dbx(media(:,:,indz_pml)).*Psi_Hx_z;

                %% Hy component
                Psi_Hy_z = bz2(ix,iy,:).*Psi_Hy_z ...
                    + cz2(ix,iy,:).*( Ex(:,:,indz_pml2) - Ex(:,:,indz_pml) );
                Psi_Hy_x = bx2(:,iy,iz).*Psi_Hy_x ...
                    + cx2(:,iy,iz).*( Ez(indx_pml2,:,:) - Ez(indx_pml,:,:) );

                Hy(indx_pml,:,:) = Hy(indx_pml,:,:) ...
                    + Dby(media(indx_pml,:,:)).*Psi_Hy_x;
                Hy(:,:,indz_pml) = Hy(:,:,indz_pml) ...
                    - Dby(media(:,:,indz_pml)).*Psi_Hy_z;

                %% Hz component
                Psi_Hz_x = bx2(:,iy,iz).*Psi_Hz_x ...
                    + cx2(:,iy,iz).*( Ey(indx_pml2,:,:) - Ey(indx_pml,:,:) );
                Psi_Hz_y = by2(ix,:,iz).*Psi_Hz_y ...
                    + cy2(ix,:,iz).*( Ex(:,indy_pml2,:) - Ex(:,indy_pml,:) );

                Hz(indx_pml,:,:) = Hz(indx_pml,:,:) ...
                    - Dbz(media(indx_pml,:,:)).*Psi_Hz_x;
                Hz(:,indy_pml,:) = Hz(:,indy_pml,:) ...
                    + Dbz(media(:,indy_pml,:)).*Psi_Hz_y;

            end;

            
            
            %% field source
            if fieldsource_flag,
                
                for l = 1:NoOfFieldSources,
                    
                    switch FieldSources{l}.component
                        
                    %% Hx
                    case 4,

                        %% if replacing, the old value is zeroed
                        if FieldSources{l}.update == 1,
                            Hx(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0;
                        end;

                        Hx(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = ...
                        Hx(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + ...
                        FieldSources{l}.mag .* inG( t + FieldSources{l}.shift );                          

                    %% Hy
                    case 5,
                        
                        %% if replacing, the old value is zeroed
                        if FieldSources{l}.update == 1,
                            Hy(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0;
                        end;

                        Hy(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = ...
                        Hy(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + ...
                        FieldSources{l}.mag .* inG( t + FieldSources{l}.shift );                          

                    %% Hz
                    case 6,

                        %% if replacing, the old value is zeroed
                        if FieldSources{l}.update == 1,
                            Hz(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0;
                        end;

                        Hz(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = ...
                        Hz(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + ...
                        FieldSources{l}.mag .* inG( t + FieldSources{l}.shift );                          
                            
                    end;
                    
                end;
                
            end;

            

            %% DFT of NFFF surfaces
            if nfff_flag,
                for p = 1:freq_ff_size,

                    %% x-y plane
                    Hy_xy_dft1(:,:,:,p) = Hy_xy_dft1(:,:,:,p) ...
                        + Wff(t,p).*Hy( x_nfff(1)+1:x_nfff(2), ...
                        y_nfff(1)+1:y_nfff(2)+1, ...
                        z_nfff );
                    Hx_xy_dft1(:,:,:,p) = Hx_xy_dft1(:,:,:,p) ...
                        + Wff(t,p).*Hx( x_nfff(1)+1:x_nfff(2)+1, ...
                        y_nfff(1)+1:y_nfff(2), ...
                        z_nfff );
                    Hy_xy_dft2(:,:,:,p) = Hy_xy_dft2(:,:,:,p) ...
                        + Wff(t,p).*Hy( x_nfff(1)+1:x_nfff(2), ...
                        y_nfff(1)+1:y_nfff(2)+1, ...
                        z_nfff + [1 1] );
                    Hx_xy_dft2(:,:,:,p) = Hx_xy_dft2(:,:,:,p) ...
                        + Wff(t,p).*Hx( x_nfff(1)+1:x_nfff(2)+1, ...
                        y_nfff(1)+1:y_nfff(2), ...
                        z_nfff + [1 1] );

                    %% y-z plane
                    Hz_yz_dft1(:,:,:,p) = Hz_yz_dft1(:,:,:,p) ...
                        + Wff(t,p).*Hz( x_nfff, ...
                        y_nfff(1)+1:y_nfff(2), ...
                        z_nfff(1)+1:z_nfff(2)+1 );
                    Hy_yz_dft1(:,:,:,p) = Hy_yz_dft1(:,:,:,p) ...
                        + Wff(t,p).*Hy( x_nfff, ...
                        y_nfff(1)+1:y_nfff(2)+1, ...
                        z_nfff(1)+1:z_nfff(2) );
                    Hz_yz_dft2(:,:,:,p) = Hz_yz_dft2(:,:,:,p) ...
                        + Wff(t,p).*Hz( x_nfff + [1 1], ...
                        y_nfff(1)+1:y_nfff(2), ...
                        z_nfff(1)+1:z_nfff(2)+1 );
                    Hy_yz_dft2(:,:,:,p) = Hy_yz_dft2(:,:,:,p) ...
                        + Wff(t,p).*Hy( x_nfff + [1 1], ...
                        y_nfff(1)+1:y_nfff(2)+1, ...
                        z_nfff(1)+1:z_nfff(2) );

                    %% z-x plane
                    Hx_zx_dft1(:,:,:,p) = Hx_zx_dft1(:,:,:,p) ...
                        + Wff(t,p).*Hx( x_nfff(1)+1:x_nfff(2)+1, ...
                        y_nfff, ...
                        z_nfff(1)+1:z_nfff(2) );
                    Hz_zx_dft1(:,:,:,p) = Hz_zx_dft1(:,:,:,p) ...
                        + Wff(t,p).*Hz( x_nfff(1)+1:x_nfff(2), ...
                        y_nfff, ...
                        z_nfff(1)+1:z_nfff(2)+1 );
                    Hx_zx_dft2(:,:,:,p) = Hx_zx_dft2(:,:,:,p) ...
                        + Wff(t,p).*Hx( x_nfff(1)+1:x_nfff(2)+1, ...
                        y_nfff + [1 1], ...
                        z_nfff(1)+1:z_nfff(2) );
                    Hz_zx_dft2(:,:,:,p) = Hz_zx_dft2(:,:,:,p) ...
                        + Wff(t,p).*Hz( x_nfff(1)+1:x_nfff(2), ...
                        y_nfff + [1 1], ...
                        z_nfff(1)+1:z_nfff(2)+1 );

                end;
            end;
            
            
            %% Probes
            if probe_flag,
                for i = 1:probe_count,
                    switch probe_type(i)
                        case 1, v.Probes(t,i) = Ex(x_probe(i),y_probe(i),z_probe(i));
                        case 2, v.Probes(t,i) = Ey(x_probe(i),y_probe(i),z_probe(i));
                        case 3, v.Probes(t,i) = Ez(x_probe(i),y_probe(i),z_probe(i));
                        case 4, v.Probes(t,i) = Hx(x_probe(i),y_probe(i),z_probe(i));
                        case 5, v.Probes(t,i) = Hy(x_probe(i),y_probe(i),z_probe(i));
                        case 6, v.Probes(t,i) = Hz(x_probe(i),y_probe(i),z_probe(i));
                    end;
                end;
            end;


            %% Energy in the domain - experimental
            if energy_flag && (t==e_comp),
                energy_now = ( ...
                    sum(sum(sum( ...
                    (Ex(indx,indy,indz).^2)./Cbx0(media(indx,indy,indz)) + ...
                    (Ey(indx,indy,indz).^2)./Cby0(media(indx,indy,indz)) + ...
                    (Ez(indx,indy,indz).^2)./Cbz0(media(indx,indy,indz)) ))) ...
                    + sum(sum(sum( ...
                    (Hx(indx,indy,indz).^2)./Dbx(media(indx,indy,indz)) + ...
                    (Hy(indx,indy,indz).^2)./Dby(media(indx,indy,indz)) + ...
                    (Hz(indx,indy,indz).^2)./Dbz(media(indx,indy,indz)) ))) );
                if lumpedsource_flag,
                    switch SourceDir %% subtraction of source
                        case 1, %% x dir
                            energy_now = energy_now ...
                                - (Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)).^2) ...
                                ./Cbx0(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)));
                        case 2, %% y dir
                            energy_now = energy_now ...
                                - (Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)).^2) ...
                                ./Cby0(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)));
                        case 3, %% z dir
                            energy_now = energy_now ...
                                - (Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)).^2) ...
                                ./Cbz0(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)));
                    end;
                end;
                energy_now = energy_now*energy_mult;
                v.energy(t) = energy_now;
                e_comp = e_comp + e_step; %% next energy computation
                
                if energy_now > energy_max, energy_max = energy_now; end;
                energy_percent = energy_now/energy_max*100;
                
            end;
            
            
            %% Energy in the domain termination condition
            if drop_flag,
                drop_cond = t > delay && energy_percent < Drop_pct;
            end;
            
            
            %% Energy in the source termination condition
            if drop2_flag,
                
                if hardsource_flag, source_now = in(t)^2;
                else                source_now = out(t)^2;
                end;
                
                if source_now > source_max, source_max = source_now; end;
                
                source_percent = source_now/source_max*100;

                if source_percent < Drop2_pct,
                    source_count = source_count + 1;
                else
                    source_count = 0;
                end;
                
                drop2_cond = t > delay && source_count > source_delay;
                
            end;

            
            %% info display
            overall = toc/t*Nt;
            disp(['n=',int2str(t),'/',int2str(Nt), ...
                '; time=',int2str(overall-toc),'/',int2str(overall), ...
                'sec; energy=',num2str(energy_now),'/', ...
                num2str(energy_percent,'%0.2f'),'%']);
            drawnow;


            %% break on termination condition
            if ( drop_flag || drop2_flag ) && ( drop_cond && drop2_cond ),
                break;
            end;
            

        end;
        %% end of cycle

        v.overall = overall; %% CPU time elapsed
        v.partial_times(2) = toc;
        tic;
        
        
        %% NFFF Radiation pattern computation
        if nfff_flag,

            D = Dx;


            %% spatial geometric averaging of H-fields

            %% x-y plane
            Hx_xy_dft = zeros(nfff_size(1)+1,nfff_size(2)  ,2,freq_ff_size);
            Hx_xy_dft1(Hx_xy_dft1==0) = 1e-20; %% avoiding div by zero
            Hx_xy_dft = sqrt( abs(Hx_xy_dft1) .* abs(Hx_xy_dft2) ) ...
                .*exp( 1j.*( angle(Hx_xy_dft1) ...
                + angle(Hx_xy_dft2./Hx_xy_dft1)./2 ) );
            Hy_xy_dft = zeros(nfff_size(1)  ,nfff_size(2)+1,2,freq_ff_size);
            Hy_xy_dft1(Hy_xy_dft1==0) = 1e-20; %% avoiding div by zero
            Hy_xy_dft = sqrt( abs(Hy_xy_dft1) .* abs(Hy_xy_dft2) ) ...
                .*exp( 1j.*( angle(Hy_xy_dft1) ...
                + angle(Hy_xy_dft2./Hy_xy_dft1)./2 ) );

            %% y-z plane
            Hy_yz_dft = zeros(2,nfff_size(2)+1,nfff_size(3)  ,freq_ff_size);
            Hy_yz_dft1(Hy_yz_dft1==0) = 1e-20; %% avoiding div by zero
            Hy_yz_dft = sqrt( abs(Hy_yz_dft1) .* abs(Hy_yz_dft2) ) ...
                .*exp( 1j.*( angle(Hy_yz_dft1) ...
                + angle(Hy_yz_dft2./Hy_yz_dft1)./2 ) );
            Hz_yz_dft = zeros(2,nfff_size(2)  ,nfff_size(3)+1,freq_ff_size);
            Hz_yz_dft1(Hz_yz_dft1==0) = 1e-20; %% avoiding div by zero
            Hz_yz_dft = sqrt( abs(Hz_yz_dft1) .* abs(Hz_yz_dft2) ) ...
                .*exp( 1j.*( angle(Hz_yz_dft1) ...
                + angle(Hz_yz_dft2./Hz_yz_dft1)./2 ) );

            %% z-x plane
            Hz_zx_dft = zeros(nfff_size(1)  ,2,nfff_size(3)+1,freq_ff_size);
            Hz_zx_dft1(Hz_zx_dft1==0) = 1e-20; %% avoiding div by zero
            Hz_zx_dft = sqrt( abs(Hz_zx_dft1) .* abs(Hz_zx_dft2) ) ...
                .*exp( 1j.*( angle(Hz_zx_dft1) ...
                + angle(Hz_zx_dft2./Hz_zx_dft1)./2 ) );
            Hx_zx_dft = zeros(nfff_size(1)+1,2,nfff_size(3)  ,freq_ff_size);
            Hx_zx_dft1(Hx_zx_dft1==0) = 1e-20; %% avoiding div by zero
            Hx_zx_dft = sqrt( abs(Hx_zx_dft1) .* abs(Hx_zx_dft2) ) ...
                .*exp( 1j.*( angle(Hx_zx_dft1) ...
                + angle(Hx_zx_dft2./Hx_zx_dft1)./2 ) );

            clear Hx_xy_dft1 Hx_xy_dft2 Hy_xy_dft1 Hy_xy_dft2 ...
                  Hy_yz_dft1 Hy_yz_dft2 Hz_yz_dft1 Hz_yz_dft2 ...
                  Hz_zx_dft1 Hz_zx_dft2 Hx_zx_dft1 Hx_zx_dft2 ;

              
            %% zeroing switched off planes
            
            %% x-y plane
            Ex_xy_dft(:,:,~v.nfff_faces(5:6),:) = 0;
            Ey_xy_dft(:,:,~v.nfff_faces(5:6),:) = 0;
            Hx_xy_dft(:,:,~v.nfff_faces(5:6),:) = 0;
            Hy_xy_dft(:,:,~v.nfff_faces(5:6),:) = 0;
            
            %% y-z plane
            Ey_yz_dft(~v.nfff_faces(1:2),:,:,:) = 0;
            Ez_yz_dft(~v.nfff_faces(1:2),:,:,:) = 0;
            Hy_yz_dft(~v.nfff_faces(1:2),:,:,:) = 0;
            Hz_yz_dft(~v.nfff_faces(1:2),:,:,:) = 0;

            %% z-x plane
            Ez_zx_dft(:,~v.nfff_faces(3:4),:,:) = 0;
            Ex_zx_dft(:,~v.nfff_faces(3:4),:,:) = 0;
            Hz_zx_dft(:,~v.nfff_faces(3:4),:,:) = 0;
            Hx_zx_dft(:,~v.nfff_faces(3:4),:,:) = 0;
            
              
            %% half-timestep shift and eta multiplication
            for p = 1:freq_ff_size,

                phase_factor = eta*exp(1j*Dt*pi*freq_ff(p));

                Hx_xy_dft(:,:,:,p) = Hx_xy_dft(:,:,:,p).*phase_factor;
                Hy_xy_dft(:,:,:,p) = Hy_xy_dft(:,:,:,p).*phase_factor;
                Hy_yz_dft(:,:,:,p) = Hy_yz_dft(:,:,:,p).*phase_factor;
                Hz_yz_dft(:,:,:,p) = Hz_yz_dft(:,:,:,p).*phase_factor;
                Hz_zx_dft(:,:,:,p) = Hz_zx_dft(:,:,:,p).*phase_factor;
                Hx_zx_dft(:,:,:,p) = Hx_zx_dft(:,:,:,p).*phase_factor;

            end;

            %% surfaces of some boundary elements are halved
            %% x-y plane
            Ex_xy_dft(:,[1 end],:,:) = Ex_xy_dft(:,[1 end],:,:) ./ 2;
            Ey_xy_dft([1 end],:,:,:) = Ey_xy_dft([1 end],:,:,:) ./ 2;
            Hx_xy_dft([1 end],:,:,:) = Hx_xy_dft([1 end],:,:,:) ./ 2;
            Hy_xy_dft(:,[1 end],:,:) = Hy_xy_dft(:,[1 end],:,:) ./ 2;
            %% y-z plane
            Ey_yz_dft(:,:,[1 end],:) = Ey_yz_dft(:,:,[1 end],:) ./ 2;
            Ez_yz_dft(:,[1 end],:,:) = Ez_yz_dft(:,[1 end],:,:) ./ 2;
            Hy_yz_dft(:,[1 end],:,:) = Hy_yz_dft(:,[1 end],:,:) ./ 2;
            Hz_yz_dft(:,:,[1 end],:) = Hz_yz_dft(:,:,[1 end],:) ./ 2;
            %% z-x plane
            Ez_zx_dft([1 end],:,:,:) = Ez_zx_dft([1 end],:,:,:) ./ 2;
            Ex_zx_dft(:,:,[1 end],:) = Ex_zx_dft(:,:,[1 end],:) ./ 2;
            Hz_zx_dft(:,:,[1 end],:) = Hz_zx_dft(:,:,[1 end],:) ./ 2;
            Hx_zx_dft([1 end],:,:,:) = Hx_zx_dft([1 end],:,:,:) ./ 2;

            %% radiation patterns in both polarizations
            v.Eth = complex(zeros(resolution(1)+1,resolution(2),freq_ff_size));
            v.Eph = complex(zeros(resolution(1)+1,resolution(2),freq_ff_size));

            %% geometrical center
            geom_center = [ (x_nfff(1) + x_nfff(2))/2 , ...
                            (y_nfff(1) + y_nfff(2))/2 , ...
                            (z_nfff(1) + z_nfff(2))/2 ].*D;

            %% difference vector between geometrical and phase centers
            r_GP = geom_center - phase_center;

            %% position vectors on surfaces - originals
            x1_orig = (x_nfff(1)    :x_nfff(2)    ).*D - geom_center(1);
            x2_orig = (x_nfff(1)+0.5:x_nfff(2)-0.5).*D - geom_center(1);
            y1_orig = (y_nfff(1)    :y_nfff(2)    ).*D - geom_center(2);
            y2_orig = (y_nfff(1)+0.5:y_nfff(2)-0.5).*D - geom_center(2);
            z1_orig = (z_nfff(1)    :z_nfff(2)    ).*D - geom_center(3);
            z2_orig = (z_nfff(1)+0.5:z_nfff(2)-0.5).*D - geom_center(3);

            %% auxiliary indices
            ones_x1 = ones(size(x1_orig));
            ones_x2 = ones(size(x2_orig));
            ones_y1 = ones(size(y1_orig));
            ones_y2 = ones(size(y2_orig));
            ones_z1 = ones(size(z1_orig));
            ones_z2 = ones(size(z2_orig));

            %% shifting dimensions
            x1_orig = shiftdim(x1_orig,1);
            x2_orig = shiftdim(x2_orig,1);
            z1_orig = shiftdim(z1_orig,-1);
            z2_orig = shiftdim(z2_orig,-1);

            %% overall coefficient
            coef = -1j.*alpha.*(D^2)./(4*pi);


            for n = 1:resolution(2),

                disp(['azimuth step = ', num2str(n)]); drawnow;

                %% phi angle
                ph = (n-1)*Dph;
                sin_ph = sin(ph);
                cos_ph = cos(ph);

                %% phi-vector components
                ph_x = -sin_ph;
                ph_y = cos_ph;
                ph_z = 0;

                for m = 1:resolution(1)+1,

                    %% theta angle
                    th = (m-1)*Dth;
                    sin_th = sin(th);
                    cos_th = cos(th);

                    %% r-vector components
                    r_x  = sin_th*cos_ph;
                    r_y  = sin_th*sin_ph;
                    r_z  = cos_th;

                    %% theta-vector components
                    th_x = cos_th*cos_ph;
                    th_y = cos_th*sin_ph;
                    th_z = -sin_th;

                    for p = 1:freq_ff_size,

                        %% position vectors - preparation
                        x1 = (alpha(p)*r_x).*x1_orig;
                        x2 = (alpha(p)*r_x).*x2_orig;
                        y1 = (alpha(p)*r_y).*y1_orig;
                        y2 = (alpha(p)*r_y).*y2_orig;
                        z1 = (alpha(p)*r_z).*z1_orig;
                        z2 = (alpha(p)*r_z).*z2_orig;


                        %% x-y plane

                        %% Hx, Ey
                        x = x1(:,ones_y2,1);
                        y = y2(ones_x1,:,1);

                        ex1 = exp(1j.*(x + y + z1(1)));
                        ex2 = exp(1j.*(x + y + z1(end)));

                        v.Eth(m,n,p) = v.Eth(m,n,p) + th_y.*( ...
                            sum(sum(Hx_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Hx_xy_dft(:,:,1,p).*ex1)) ) ...
                            + ph_x.*( ...
                            sum(sum(Ey_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Ey_xy_dft(:,:,1,p).*ex1)) );

                        v.Eph(m,n,p) = v.Eph(m,n,p) + ph_y.*( ...
                            sum(sum(Hx_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Hx_xy_dft(:,:,1,p).*ex1)) ) ...
                            - th_x.*( ...
                            sum(sum(Ey_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Ey_xy_dft(:,:,1,p).*ex1)) );

                        %% Hy, Ex
                        x = x2(:,ones_y1,1);
                        y = y1(ones_x2,:,1);

                        ex1 = exp(1j.*(x + y + z1(1)));
                        ex2 = exp(1j.*(x + y + z1(end)));

                        v.Eth(m,n,p) = v.Eth(m,n,p) - th_x.*( ...
                            sum(sum(Hy_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Hy_xy_dft(:,:,1,p).*ex1)) ) ...
                            - ph_y.*( ...
                            sum(sum(Ex_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Ex_xy_dft(:,:,1,p).*ex1)) );

                        v.Eph(m,n,p) = v.Eph(m,n,p) - ph_x.*( ...
                            sum(sum(Hy_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Hy_xy_dft(:,:,1,p).*ex1)) ) ...
                            + th_y.*( ...
                            sum(sum(Ex_xy_dft(:,:,2,p).*ex2)) - ...
                            sum(sum(Ex_xy_dft(:,:,1,p).*ex1)) );


                        %% y-z plane

                        %% Hy, Ez
                        y = y1(1,:,ones_z2);
                        z = z2(1,ones_y1,:);

                        ex1 = exp(1j.*(x1(1)   + y + z));
                        ex2 = exp(1j.*(x1(end) + y + z));

                        v.Eth(m,n,p) = v.Eth(m,n,p) + th_z.*( ...
                            sum(sum(Hy_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Hy_yz_dft(1,:,:,p).*ex1)) ) ...
                            + ph_y.*( ...
                            sum(sum(Ez_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Ez_yz_dft(1,:,:,p).*ex1)) );

                        v.Eph(m,n,p) = v.Eph(m,n,p) + ph_z.*( ...
                            sum(sum(Hy_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Hy_yz_dft(1,:,:,p).*ex1)) ) ...
                            - th_y.*( ...
                            sum(sum(Ez_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Ez_yz_dft(1,:,:,p).*ex1)) );

                        %% Hz, Ey
                        y = y2(1,:,ones_z1);
                        z = z1(1,ones_y2,:);

                        ex1 = exp(1j.*(x1(1)   + y + z));
                        ex2 = exp(1j.*(x1(end) + y + z));

                        v.Eth(m,n,p) = v.Eth(m,n,p) - th_y.*( ...
                            sum(sum(Hz_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Hz_yz_dft(1,:,:,p).*ex1)) ) ...
                            - ph_z.*( ...
                            sum(sum(Ey_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Ey_yz_dft(1,:,:,p).*ex1)) );

                        v.Eph(m,n,p) = v.Eph(m,n,p) - ph_y.*( ...
                            sum(sum(Hz_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Hz_yz_dft(1,:,:,p).*ex1)) ) ...
                            + th_z.*( ...
                            sum(sum(Ey_yz_dft(2,:,:,p).*ex2)) - ...
                            sum(sum(Ey_yz_dft(1,:,:,p).*ex1)) );


                        %% z-x plane

                        %% Hz, Ex
                        z = z1(ones_x2,1,:);
                        x = x2(:,1,ones_z1);

                        ex1 = exp(1j.*(x + y1(1)   + z));
                        ex2 = exp(1j.*(x + y1(end) + z));

                        v.Eth(m,n,p) = v.Eth(m,n,p) + th_x.*( ...
                            sum(sum(Hz_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Hz_zx_dft(:,1,:,p).*ex1)) ) ...
                            + ph_z.*( ...
                            sum(sum(Ex_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Ex_zx_dft(:,1,:,p).*ex1)) );

                        v.Eph(m,n,p) = v.Eph(m,n,p) + ph_x.*( ...
                            sum(sum(Hz_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Hz_zx_dft(:,1,:,p).*ex1)) ) ...
                            - th_z.*( ...
                            sum(sum(Ex_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Ex_zx_dft(:,1,:,p).*ex1)) );

                        %% Hx, Ez
                        z = z2(ones_x1,1,:);
                        x = x1(:,1,ones_z2);

                        ex1 = exp(1j.*(x + y1(1)   + z));
                        ex2 = exp(1j.*(x + y1(end) + z));

                        v.Eth(m,n,p) = v.Eth(m,n,p) - th_z.*( ...
                            sum(sum(Hx_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Hx_zx_dft(:,1,:,p).*ex1)) ) ...
                            - ph_x.*( ...
                            sum(sum(Ez_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Ez_zx_dft(:,1,:,p).*ex1)) );

                        v.Eph(m,n,p) = v.Eph(m,n,p) - ph_z.*( ...
                            sum(sum(Hx_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Hx_zx_dft(:,1,:,p).*ex1)) ) ...
                            + th_x.*( ...
                            sum(sum(Ez_zx_dft(:,2,:,p).*ex2)) - ...
                            sum(sum(Ez_zx_dft(:,1,:,p).*ex1)) );


                        %% phase of the G-P difference
                        ex1 = exp( 1j*alpha(p)*( r_GP(1)*r_x + r_GP(2)*r_y + r_GP(3)*r_z ) );

                        %% final multiplication
                        v.Eth(m,n,p) = v.Eth(m,n,p)*coef(p)*ex1;
                        v.Eph(m,n,p) = v.Eph(m,n,p)*coef(p)*ex1;

                    end; %% end p

                end; %% end m (th)

            end; %% end n (ph)

        end; %% end nfff
        
        v.Rs = Rs;
        
        v.partial_times(3) = toc;
        v.partial_times(4) = 0;
        v.time_end = clock;

    end %% end MatlabFDTD

    
    
    %% Fortran FDTD deployment to file
    function FortranDeploy

        if pref.CacheOptim, % cache optimizations for opteron
            
            precision =  8;      %% double precision: 8 bytes
            overhead  = 80;      %% allocation overhead in bytes: 80
            pagesize  =  4*1024; %% memory page size in bytes: 4k
            cachesize = 64*1024; %% cache size in bytes: 64k
            assoc     =  2;      %% cache associativity: 2-way
            cacheline = 64;      %% length of the cache line in bytes: 64

            if ~pref.NewInOutFiles,
            
                cachepages = cachesize/assoc/pagesize; %% unique pages in cache
                if cachepages~=round(cachepages),
                    error('Cache size not multiple of page size.');
                end;

                %% padding arrays
                pad  = [0 0 0];
                padf = [0 0 0];

                %% Ex, Ey, Ez, Hx, Hy, Hz padding
                count = 0;
                %% index from the largest dimension to the shortest
                [dummy,ind] = sort([xmax_pml;ymax_pml;zmax_pml],'descend');

                while 1,

                    %% span across four arrays
                    span = 4.*ceil((( ...
                        (xmax_pml + pad(1)).*(ymax_pml + pad(2)).*(zmax_pml + pad(3)) ...
                        ).*precision + overhead)./pagesize);

                    %% problem occurs when span is divisible
                    problem = mod(span,cachepages)==0;

                    if ~problem, break; end;

                    pad = [0 0 0];

                    %% trying to increase dimension by 1 from the largest to the shortest
                    %% (overhead increases), if not success, do the same with increase by
                    %% 2, etc...
                    pad(ind(rem(count,3)+1)) = floor(count/3)+1;
                    count = count + 1;

                end;

                %% nearfield padding
                if field_flag,

                    count = 0;
                    %% index from the largest dimension to the shortest
                    [dummy,ind] = sort([xmax;ymax;zmax],'descend');

                    while 1,

                        %% span across two complex arrays
                        span = 2.*ceil((( ...
                            (xmax + padf(1)).*(ymax + padf(2)).*(zmax + padf(3)).*freq_nf_size ...
                            ).*(2*precision) + overhead)./pagesize);

                        %% problem occurs when span is divisible
                        problem = mod(span,cachepages)==0;

                        if ~problem, break; end;

                        padf = [0 0 0];

                        %% trying to increase dimension by 1 from the largest to the shortest
                        %% (overhead increases), if not success, do the same with increase by
                        %% 2, etc...
                        padf(ind(rem(count,3)+1)) = floor(count/3)+1;
                        count = count + 1;

                    end;

                end;

            end;

        end;

        %% name of the file serving as input for Fortran FDTD code
        FileName = strcat(ActualFileName,'.in');

        %% making the file
        [fid,message] = fopen(fullfile(DeployDir,FileName),'w');
        if isequal(fid,-1),
            errordlg({'Error while opening .IN file',message},'Open error','modal');
            uiwait;
            return;
        end;
        
        if pref.NewInOutFiles, % new deploy
            
            %% number of field sources
            fieldsource_count = 0; max_shift = 0;
            for i = 1:length(v.Geometry),
                if strcmp(v.Geometry{i}.type,'field source'),
                    fieldsource_count = fieldsource_count + 1;
                    phase = angle( v.Geometry{i}.data );
                    %% phase should be always negative
                    phase(phase>0) = phase(phase>0) - 2*pi;
                    max_shift = max( max_shift, max(abs( phase(:)./(2*pi*fc*v.Dt) )) );
                end;
            end;
            %% find the highest shift and add to delay for the termination
            %% condition
            max_shift = ceil(max_shift);
            delay = delay + max_shift;

            offset_Dt = 0;
            
            if fieldsource_count > 0,

                %% calculate DFT of Vs at fc and divide (normalize) the
                %% magnitudes
                Nt_fft = round( 1/(v.Dt*v.Df) );
                Vs_fft = 2.*fft(v.Vs,Nt_fft)./Nt_fft;
                f_sample = round(fc./v.Df) + 1;
                mag_norm = abs(Vs_fft(f_sample));
                
                %% offset of the center of the pulse in time steps
                offset_Dt = round( offset/v.Dt );
                
            end;

            cache_flag = pref.CacheOptim;
            muaver_flag = v.MuAveraging;
            
            %% common part
            fwrite(fid,int32([v.limits(1,:),v.limits(2,:),pmax, ...
                nfff_flag,v.nfff_faces,field_flag, ...
                e_step,drop_flag,drop2_flag,cache_flag,muaver_flag, ...
                probe_count,fieldsource_count, ...
                Nt,x_nfff,y_nfff,z_nfff,freq_ff_size,freq_nf_size]),'int32');
            fwrite(fid,[Dx,Dt,pi,DefaultSourceResistance,eps0,mu0, ...
                fc,fs,v.mag,offset_Dt],'double');
            fwrite(fid,in,'double');

            %% optional parts
            if pml_flag,
                fwrite(fid,[sig_max,kappa_max,a_max,m_s,m_k,m_a],'double');
            end;
            if nfff_flag,
                fwrite(fid,resolution,'int32');
                fwrite(fid,[eta,c,phase_center,freq_ff],'double');
            end;
            if field_flag,
                fwrite(fid,freq_nf,'double');
            end;
            if drop_flag || drop2_flag,
                fwrite(fid,delay,'int32');
            end;
            if drop_flag,
                fwrite(fid,v.Drop_dB,'double');
            end;
            if drop2_flag,
                fwrite(fid,v.Drop2_dB,'double');
                fwrite(fid,source_delay,'int32');
            end;
            if cache_flag,
                fwrite(fid,[precision,overhead,pagesize,cachesize,assoc,cacheline],'int32');
            end;

            for i = 1:length(v.Geometry),
                switch v.Geometry{i}.type

                    %% WIRE
                    case 'wire'

                        fwrite(fid,[1,48],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.point2],'double');

                    %% TRIANGLE
                    case 'triangle'

                        fwrite(fid,[2,72],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.point2, ...
                                    v.Geometry{i}.point3],'double');

                    %% RECTANGLE
                    case 'rectangle'

                        fwrite(fid,[3,72],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.point2, ...
                                    v.Geometry{i}.point3],'double');

                    %% BRICK
                    case 'brick'

                        eps_r = v.Geometry{i}.permittivity;
                        sig   = v.Geometry{i}.conductivity;
                        if v.Geometry{i}.PEC || isinf(eps_r), eps_r = 0; end; % PEC flag

                        fwrite(fid,[4,64],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.point2, ...
                                    eps_r,sig],'double');

                    %% SOURCE
                    case 'source'

                        fwrite(fid,[5,40],'int32');
                        fwrite(fid,[v.Geometry{i}.gndpoint, ...
                                    v.Geometry{i}.direction, ...
                                    v.Geometry{i}.resistance],'double');
                                
                    %% RESISTOR
                    case 'resistor'

                        fwrite(fid,[6,56],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.point2, ...
                                    v.Geometry{i}.resistance],'double');

                    %% CAPACITOR
                    case 'capacitor'

                        fwrite(fid,[7,56],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.point2, ...
                                    v.Geometry{i}.capacitance],'double');

                    %% EXTERNAL OBJECT
                    case 'extobj'

                        sm = zeros(4,3); % size matrix

                        [sm(1,1),sm(1,2),sm(1,3)] = size(v.Geometry{i}.VolArray);
                        [sm(2,1),sm(2,2),sm(2,3)] = size(v.Geometry{i}.PlxArray);
                        [sm(3,1),sm(3,2),sm(3,3)] = size(v.Geometry{i}.PlyArray);
                        [sm(4,1),sm(4,2),sm(4,3)] = size(v.Geometry{i}.PlzArray);

                        flags = ~~prod(sm,2); % flags non-empty arrays

                        if ~flags, % if all empty
                            h = errordlg({'External object found empty!', ...
                                'Aborting deployment.'},'Error','modal');
                            uiwait;
                            fclose(fid);
                            return;
                        end;

                        % normalize the sizes of the plate arrays, if not empty
                        sm(2:4,:) = sm(2:4,:) - diag(flags(2:4));

                        sm = sm(flags,:); % only non-empty sizes

                        smd = diff(sm,1,1); % differences between sizes

                        if ~isempty(smd) && any(smd(:)~=0),
                            h = errordlg({'External object arrays do not have corresponding sizes!', ...
                                'Aborting deployment.'},'Error','modal');
                            uiwait;
                            fclose(fid);
                            return;
                        end;

                        ArraySize = sm(1,:);
                        ArrayFlags = flags.';

                        VolMaxValue = 0;
                        PlMaxValue = 0;

                        if flags(1),
                            VolMaxValue = size(v.Geometry{i}.VolEps,1);
                        end;
                        if any(flags(2:4)),
                            PlMaxValue = size(v.Geometry{i}.PlMetal,1);
                        end;

                        %% number of materials cannot exceed capacity of INT8
                        if VolMaxValue > 127 || PlMaxValue > 127,
                            h = errordlg({'Too many materials in external object.', ...
                                'Aborting deployment.'},'Error','modal');
                            uiwait;
                            fclose(fid);
                            return;
                        end;

                        %% when there's no source
                        if isempty(v.Geometry{i}.GNDPoint) || ...
                           isempty(v.Geometry{i}.Direction),
                            SourceDir = 0;
                        else
                            SourceDir = v.Geometry{i}.Direction;
                        end;

                        sm = zeros(9,2);

                        if ~isempty(v.Geometry{i}.Tweaks),
                            [sm(1,1),sm(1,2)] = size(v.Geometry{i}.Tweaks.eps_rx);
                            [sm(2,1),sm(2,2)] = size(v.Geometry{i}.Tweaks.eps_ry);
                            [sm(3,1),sm(3,2)] = size(v.Geometry{i}.Tweaks.eps_rz);
                            [sm(4,1),sm(4,2)] = size(v.Geometry{i}.Tweaks.mu_rx);
                            [sm(5,1),sm(5,2)] = size(v.Geometry{i}.Tweaks.mu_ry);
                            [sm(6,1),sm(6,2)] = size(v.Geometry{i}.Tweaks.mu_rz);
                            [sm(7,1),sm(7,2)] = size(v.Geometry{i}.Tweaks.sig_x);
                            [sm(8,1),sm(8,2)] = size(v.Geometry{i}.Tweaks.sig_y);
                            [sm(9,1),sm(9,2)] = size(v.Geometry{i}.Tweaks.sig_z);
                        end;

                        flags = ~~prod(sm,2); % flags non-empty arrays

                        TweaksSizes = sm(:,1).*flags;

                        sm = sm(flags,2); % only widths of non-empty

                        if ~isempty(sm) && any(sm~=5),
                            h = errordlg({'External object tweaks do not have correct sizes!', ...
                                'Aborting deployment.'},'Error','modal');
                            uiwait;
                            fclose(fid);
                            return;
                        end;

                        %% length of the jump over this record
                        JumpLength = ...
                            3*4 + ... % position point
                            4*4 + ... % 4 flags
                            3*4 + ... % array size
                            prod(ArraySize)*ArrayFlags(1) + ... % volume array
                            prod(ArraySize+[1,0,0])*ArrayFlags(2) + ...
                            prod(ArraySize+[0,1,0])*ArrayFlags(3) + ...
                            prod(ArraySize+[0,0,1])*ArrayFlags(4) + ...
                            ArrayFlags(1)*(VolMaxValue*3*8) + ... % volume LUTs
                            any(ArrayFlags(2:4))*(PlMaxValue*4) + ... % metal LUT
                            4 + ... % source direction
                            (~~SourceDir)*3*4 + ... % source GND point
                            9*4 + ... % tweaks sizes
                            sum(TweaksSizes)*5*8; % tweaks 

                        fwrite(fid,[8,JumpLength],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                            (ArrayFlags.*[VolMaxValue,PlMaxValue,PlMaxValue,PlMaxValue]), ...
                            ArraySize],'int32');

                        %% writing arrays
                        if ArrayFlags(1), fwrite(fid,v.Geometry{i}.VolArray,'int8'); end;
                        if ArrayFlags(2), fwrite(fid,v.Geometry{i}.PlxArray,'int8'); end;
                        if ArrayFlags(3), fwrite(fid,v.Geometry{i}.PlyArray,'int8'); end;
                        if ArrayFlags(4), fwrite(fid,v.Geometry{i}.PlzArray,'int8'); end;

                        %% volume array lookup tables
                        if ArrayFlags(1),
                            VolEps = v.Geometry{i}.VolEps;
                            VolEps(isinf(VolEps)) = 0; % cannot pass Inf
                            fwrite(fid,VolEps,'double');
                            fwrite(fid,v.Geometry{i}.VolMu ,'double');
                            fwrite(fid,v.Geometry{i}.VolSig,'double');
                        end;

                        %% metal lookup table
                        if any(ArrayFlags(2:4)),
                            fwrite(fid,~~v.Geometry{i}.PlMetal,'int32');
                        end;

                        %% source direction -- if 0 then there's no source
                        fwrite(fid,SourceDir,'int32');

                        %% source GND point coordinates
                        if SourceDir,
                            fwrite(fid,v.Geometry{i}.GNDPoint,'int32');
                        end;

                        %% sizes of tweaks -- can also be all zero
                        fwrite(fid,TweaksSizes,'int32');

                        %% tweaks if exist
                        if TweaksSizes(1), fwrite(fid,v.Geometry{i}.Tweaks.eps_rx.','double'); end;
                        if TweaksSizes(2), fwrite(fid,v.Geometry{i}.Tweaks.eps_ry.','double'); end;
                        if TweaksSizes(3), fwrite(fid,v.Geometry{i}.Tweaks.eps_rz.','double'); end;
                        if TweaksSizes(4), fwrite(fid,v.Geometry{i}.Tweaks.mu_rx.' ,'double'); end;
                        if TweaksSizes(5), fwrite(fid,v.Geometry{i}.Tweaks.mu_ry.' ,'double'); end;
                        if TweaksSizes(6), fwrite(fid,v.Geometry{i}.Tweaks.mu_rz.' ,'double'); end;
                        if TweaksSizes(7), fwrite(fid,v.Geometry{i}.Tweaks.sig_x.' ,'double'); end;
                        if TweaksSizes(8), fwrite(fid,v.Geometry{i}.Tweaks.sig_y.' ,'double'); end;
                        if TweaksSizes(9), fwrite(fid,v.Geometry{i}.Tweaks.sig_z.' ,'double'); end;

                    %% PROBE
                    case 'probe'

                        fwrite(fid,[9,16],'int32');
                        fwrite(fid,[v.Geometry{i}.point1, ...
                                    v.Geometry{i}.component],'int32');

                    %% FIELD SOURCE
                    case 'field source'
                        
                        %% phase should be always negative
                        phase = angle( v.Geometry{i}.data );
                        phase(phase>0) = phase(phase>0) - 2*pi;

                        ArraySize = prod(size(v.Geometry{i}.data));

                        JumpLength = ...
                            6*4 +           ... % corner points
                            2*8*ArraySize + ... % magnitude and shift
                            2*4;                % component and update type
                        
                        fwrite(fid,[10,JumpLength],'int32');
                        
                        %% corner points
                        fwrite(fid,[ v.Geometry{i}.point1 , ...
                                     v.Geometry{i}.point2 ],'int32');
                                 
                        %% magnitude
                        fwrite(fid,abs(v.Geometry{i}.data)./mag_norm,'double');
                        
                        %% shift
                        fwrite(fid,phase./(2*pi*fc*v.Dt),'double');
                        
                        %% component and update type
                        fwrite(fid,[ v.Geometry{i}.component , ...
                                     v.Geometry{i}.update ],'int32');

                end; % end element switch
            end; % end for

            fwrite(fid,[0,0],'int32'); % end marks
            fclose(fid);
            
        else % old deploy
        
            %% length of the lookup tables
            mat_len = length(Cax);
            media_class = 'int32';

            fwrite(fid,int32([xmax_pml,ymax_pml,zmax_pml,pmax, ...
                nfff_flag,field_flag,e_step,drop_flag,drop2_flag,mat_len,Nt]),'int32');
            fwrite(fid,int32([SourceSubs,SourceDir, ...
                x_nfff,y_nfff,z_nfff,freq_ff_size]),'int32');
            fwrite(fid,[Dx Dt],'double');
            fwrite(fid,pi,'double');
            fwrite(fid,[pad padf],'int32');
            fwrite(fid,[Cax Cay Caz Cbx Cby Cbz Dbx Dby Dbz],'double');
            fwrite(fid,media,media_class);
            fwrite(fid,in,'double');
            if pml_flag, fwrite(fid,[bx1 bx2 cx1 cx2],'double'); end;
            if nfff_flag,
                fwrite(fid,[eta c],'double');
                fwrite(fid,int32(resolution),'int32');
                fwrite(fid,phase_center,'double');
                fwrite(fid,freq_ff,'double');
            end;
            if field_flag,
                fwrite(fid,freq_nf,'double');
            end;
            if drop_flag || drop2_flag,
                fwrite(fid,int32(delay),'int32');
            end;
            if drop_flag,
                Drop_pct = 100*10^(v.Drop_dB/10); %% in percent
                fwrite(fid,Drop_pct,'double');
            end;
            if drop2_flag,
                Drop2_pct = 100*10^(v.Drop2_dB/10); %% in percent
                fwrite(fid,Drop2_pct,'double');
                fwrite(fid,int32(source_delay),'int32');
            end;
            fwrite(fid,int32(hardsource_flag),'int32');
            fclose(fid);

        end;
        
        %% writing the file name into the queue
        %% (opening file in text mode and append)
        [fid,message] = fopen(fullfile(DeployDir,TasklistName),'at');
        if isequal(fid,-1),
            errordlg({'Error while opening tasklist file',message},'Open error','modal');
            uiwait;
            return;
        end;
        
        fprintf(fid,'%s\n',FileName);
        fclose(fid);
        
        SuccessFlag = true;
        
    end % end FortranDeploy

end %% end Simulate


    
function SimulationDeploy(src,eventdata)

    option = src;
    
    %% if deploy multi
    switch option, case MenuDeployMulti
    
        %% which files to take
        [FileNames,PathName] = uigetfile( ...
            {'*.aaf','AAU3 FDTD Files'; ...
             '*','All Files (*.*)'}, ...
            'Deploy Multiple Files for Fortran Processing',[WorkDir,filesep], ...
            'MultiSelect','on');

        if isequal(FileNames,0) || isequal(PathName,0), return; end;
        WorkDir = PathName;
        
        if ~iscell(FileNames), FileNames = cellstr(FileNames); end;

    end;
    
    %% the window
    DeployDialog = dialog('Visible','off', ...
        'Name','Job Folder', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenSize(3:4)./2-[180 40] 360 80], ...
        'KeyPressFcn',{@RetEscPressed,@DeployProceed,@CloseFigure});

    %% the job folder edit
    uitext(DeployDialog,'Deploy file(s) to folder:',[10 45 120 17]);
    JobDirEdit = uiedit(DeployDialog,Engine.DefJobFolder, ...
        [130 45 195 21],'HorizontalAlignment','left');
    JobDirPush = uipushbutton(DeployDialog,'...',[329 45 21 21], ...
        {@GetDir,JobDirEdit,'Deploy Path...'});

    %% if the directory is remote, disable the pushbutton
    switch Engine.Location, case 'Remote'
        set(JobDirPush,'Enable','off');
    end;
    
    %% default behavior of Enter and Escape on controls
    set(JobDirEdit,'KeyPressFcn',{@RetEscPressed,@DeployProceed,@CloseFigure});

    %% buttons
    OKButton     = Button(DeployDialog,'OK'    ,2,@DeployProceed);
    CancelButton = Button(DeployDialog,'Cancel',1,@CloseFigure);

    uicontrol(JobDirEdit); %% give focus to the edit field
    set(DeployDialog,'Visible','on'); %% display figure

    function DeployProceed(src,eventdata)

        JobDir = get(JobDirEdit,'String');

        %% check if the directory is empty
        if isempty(JobDir),
            errordlg('Input cannot be empty','Bad input','modal');
            pause(1);
            uicontrol(JobDirEdit);
            return;
        end;

        %% if Local, check if the directory is valid
        if strcmp(Engine.Location,'Local') && ~isdir(JobDir),
            errordlg('Path must point to a valid local folder', ...
                     'Bad input','modal');
            pause(1);
            uicontrol(JobDirEdit);
            return;
        end;

        %% close the window
        close(DeployDialog);

        switch Engine.Location

            case 'Local',
                DeployDir = JobDir;
                
            case 'Remote'
                DeployDir = Engine.LocalDeployFolder;
                
                %% draw the tasklist file from the remote server
                %% (if it does not exist, it does not matter)
                PuttyDownload(fullfile(JobDir,'tasklist.txt'),DeployDir);
                
        end;
        
        %% make sure that the deploy folder exists (and suppress warning)
        [null,null,null] = mkdir(DeployDir);
        
        %% if deploy multi
        switch option, case MenuDeployMulti
            
            SuccessCount = 0;

            for n = 1:length(FileNames),

                Open(MenuDeployMulti,FileNames{n});

                if v.locked_flag,
                    h = msgbox({['File ',FileNames{n},' already contains results.'], ...
                                'Not deploying.'},'Error','error','modal');
                    uiwait(h);
                    continue;
                end;

                %% add the job folder to v and save it to the AAF file
                v.JobDir = JobDir;
                Save;

                SuccessFlag = false;

                %% call the simulation routine
                Simulate([],[],'fortran',DeployDir);

                %% if remote...
                switch Engine.Location, case 'Remote',
                    
                    %% first make sure that the job folder exists
                    PuttyRun(['mkdir -p "',JobDir,'"']);
                    
                    %% upload the file and delete it
                    PuttyUpload(fullfile(DeployDir,[ActualFileName,'.in']),JobDir);
                    delete(     fullfile(DeployDir,[ActualFileName,'.in']));
                    
                end;

                if SuccessFlag, SuccessCount = SuccessCount + 1; end;

            end;
            
            %% if remote, upload the tasklist and delete it
            switch Engine.Location, case 'Remote',
                PuttyUpload(fullfile(DeployDir,'tasklist.txt'),JobDir);
                delete(     fullfile(DeployDir,'tasklist.txt'));
            end;

            msgbox([num2str(SuccessCount),' file(s) successfully deployed.'],'Info','modal');

        otherwise

            %% add the job folder to v and save it to the AAF file
            v.JobDir = JobDir;
            Save;

            %% call the simulation routine
            Simulate([],[],'fortran',DeployDir);

            %% if remote...
            switch Engine.Location, case 'Remote',

                %% first make sure that the job folder exists
                PuttyRun(['mkdir -p "',JobDir,'"']);
                
                %% upload the deployed file and tasklist into the job folder
                PuttyUpload(fullfile(DeployDir,[ActualFileName,'.in']),JobDir);
                PuttyUpload(fullfile(DeployDir,'tasklist.txt'),JobDir);

                %% delete the files
                delete(fullfile(DeployDir,[ActualFileName,'.in']));
                delete(fullfile(DeployDir,'tasklist.txt'));

            end;

        end;
        
    end;

end



function SimulationStart(src,eventdata)

    %% if Matlab engine, simulate directly
    if pref.ActiveEngine == 1,
        Simulate([],[],'matlab','');
        return;
    end;

    %% original path of the script
    OrigScriptPath = Engine.RunScript;

    switch Engine.Script
        
        case 'Ready to use'
            
            %% path to the script unchanged
            ScriptPath = OrigScriptPath;
            
            %% go directly to proceed
            SimulationProceed;
            
        case 'Edit from template'
            
            %% prepare the filename from parts
            [null,FileName,Extension] = fileparts(OrigScriptPath);
            ScriptName = [FileName,Extension];
            
            %% script name in job folder
            ScriptPath_JobDir = fullfile(v.JobDir,ScriptName);

            switch Engine.Location
                
                case 'Local'
                    
                    %% the script template will be in the original position
                    ScriptPath = OrigScriptPath;
                    
                    %% if option for checking if the script is in the job
                    %% folder and it is really there, then take it from there
                    if Engine.CheckScript && exist(ScriptPath_JobDir,'file'),
                        ScriptPath = ScriptPath_JobDir;
                    end;

                case 'Remote'
                    
                    %% define the deploy folder
                    DeployDir = Engine.LocalDeployFolder;
                    
                    %% the script template will be in the deploy folder
                    ScriptPath = fullfile(DeployDir,ScriptName);                   

                    %% if option for checking if the script is in the job
                    %% folder, then try to download it locally
                    if Engine.CheckScript,
                        PuttyDownload(ScriptPath_JobDir,DeployDir);
                    end;
                    
                    %% if the template is not yet in the deploy folder, 
                    %% download the default one instead
                    if ~exist(ScriptPath,'file'),
                        PuttyDownload(OrigScriptPath,DeployDir);
                    end;

            end;
            
            %% open the file with the script template
            ScriptString = fileread(ScriptPath);

            %% the window
            ScriptDialog = dialog('Visible','off', ...
                'Name','Edit Script', ...
                'Units','pixels', ...
                'Position',[ScreenSize(3:4)./2-[300 200] 600 400], ...
                'KeyPressFcn',{@RetEscPressed,@SaveScript,@CloseFigure});

            ScriptEdit = uimultiedit(ScriptDialog,ScriptString,[10 45 580 345], ...
                'FontName','FixedWidth');

            %% buttons
            OKButton     = Button(ScriptDialog,'OK'    ,2,@SaveScript);
            CancelButton = Button(ScriptDialog,'Cancel',1,@CloseFigure);

            uicontrol(ScriptEdit); %% give focus to the edit field

            set(ScriptDialog,'Visible','on'); %% display figure
            
    end


    function SaveScript(src,eventdata)
        
        %% retrieve the string from the edit box
        ScriptString = get(ScriptEdit,'String');
        
        %% change the location where the script will be written
        switch Engine.Location
            case 'Local',  ScriptPath = fullfile( v.JobDir,ScriptName);
            case 'Remote', ScriptPath = fullfile(DeployDir,ScriptName);
        end;

        %% open the edited script in the job folder
        [fid,message] = fopen(ScriptPath,'w');
        if isequal(fid,-1),
            errordlg({'Error while opening the script file',message},'Open error','modal');
            uiwait;
            return;
        end;

        %% go line by line and write the string to the file; then close
        for i = 1:size(ScriptString,1),
            fprintf(fid,'%s\n',ScriptString(i,:));
        end;
        fclose(fid);
        
        %% on unix, set the attributes to executable
        if isunix, fileattrib(ScriptPath,'+x'); end;

        %% if remote, we need to upload it and clean
        switch Engine.Location, case 'Remote'
            
            PuttyUpload(ScriptPath,v.JobDir);
            
            %% set the attributes to executable
            switch Engine.Processing, case 'Run immediately'
                PuttyRun(['chmod 100 ',fullfile(v.JobDir,ScriptName)]);
            end;
                
            delete(ScriptPath);
            ScriptPath = ScriptPath_JobDir;
            
        end;

        %% close the window
        close(ScriptDialog);
        
        %% go to proceed
        SimulationProceed;
        
    end


    function SimulationProceed
        
        %% Processing type
        switch Engine.Processing
            case 'Run immediately' %% run in background
                CommandSequence = ['"',ScriptPath,'" &'];
            case 'Submit to queue' %% submit command first
                CommandSequence = ['cd ',v.JobDir,'; ', ...
                    Engine.SubmitCommand,' "',ScriptPath,'"'];
        end;
        
        %% Location
        switch Engine.Location
            
            case 'Local'
                
                %% save the original folder and switch to the job folder
                OrigDir = cd;
                cd(v.JobDir);

                %% start the command/script
                if     ispc,   system(                       CommandSequence);
                elseif isunix, system(['xterm -rv -hold -e ',CommandSequence]);
                else   warning('Where are we?');
                end;

                %% restore the original folder
                cd(OrigDir);

            case 'Remote'

                %% run remotely
                PuttyRun(['xterm -rv -hold -e ''',CommandSequence,'''']);

        end;

    end;
        
end



function SimulationCollect(src,eventdata)
    
    %% flags
    nfff_flag   = v.nfff_flag;
    field_flag  = v.field_flag;
    energy_flag = v.energy_flag;
    
    %% number of probes
    probe_count = 0;
    for i = 1:length(v.Geometry),
        if strcmp(v.Geometry{i}.type,'probe'),
            probe_count = probe_count + 1;
        end;
    end;
    probe_flag = probe_count > 0;

    Nt = v.Nt;
    resolution = v.resolution; %% pattern resolution in theta and phi

    %% lattice dimensions (excluding PML)
    xmax  = abs( v.limits(2,1) - v.limits(1,1) );
    ymax  = abs( v.limits(2,2) - v.limits(1,2) );
    zmax  = abs( v.limits(2,3) - v.limits(1,3) );
    nfmax = length(v.freq_nf);
    ffmax = length(v.freq_ff);
    
    x_nfff = v.x_nfff; %% x-distance from PML of the NFFF int. surface [(D)]
    y_nfff = v.y_nfff; %% y-distance from PML of the NFFF int. surface [(D)]
    z_nfff = v.z_nfff; %% z-distance from PML of the NFFF int. surface [(D)]
    
    patt_length = prod(resolution + [1 0])*ffmax;
    field_length = xmax*ymax*zmax*nfmax;

    %% name of the file serving as output from the Fortran FDTD code
    FileName = [ActualFileName,'.in.out'];
    
    switch Engine.Location
        
        case 'Local' %% open from the job folder
    
            [fid,message] = fopen(fullfile(v.JobDir,FileName),'r');

        case 'Remote' %% download first to the deploy folder
            
            DeployDir = Engine.LocalDeployFolder;
            PuttyDownload(fullfile(v.JobDir,FileName),DeployDir);
            [fid,message] = fopen(fullfile(DeployDir,FileName),'r');
    
    end;

    if isequal(fid,-1),
        errordlg({'Error while opening .OUT file',message},'Open error','modal');
        uiwait;
        return;
    end;
    
    out = fread(fid,Nt,'double');
    overall = fread(fid,1,'int32');
    v.overall = double(overall);

    if nfff_flag,
        Eth_real = fread(fid,patt_length,'double');
        Eth_imag = fread(fid,patt_length,'double');
        Eph_real = fread(fid,patt_length,'double');
        Eph_imag = fread(fid,patt_length,'double');
    end;

    if field_flag,
        Ex_field_real = fread(fid,field_length,'double');
        Ex_field_imag = fread(fid,field_length,'double');
        Ey_field_real = fread(fid,field_length,'double');
        Ey_field_imag = fread(fid,field_length,'double');
        Ez_field_real = fread(fid,field_length,'double');
        Ez_field_imag = fread(fid,field_length,'double');
    end;

    if energy_flag,
        v.energy = fread(fid,Nt,'double');
    end;
    
    if pref.NewInOutFiles,
        if probe_flag,
            v.Probes = fread(fid,Nt*probe_count,'double');
            v.Probes = reshape(v.Probes,[Nt,probe_count]);
        end;
        v.Rs = fread(fid,1,'double');
    else
        v.Rs = DefaultSourceResistance;
    end;
    
    v.time_start = fread(fid,8,'int32');
    v.time_end   = fread(fid,8,'int32');
    v.partial_times = fread(fid,4,'int32');

    fclose(fid);
    
    %% if remote operation, delete the local out file
    switch Engine.Location, case 'Remote'
        delete(fullfile(DeployDir,FileName));
    end;

    if nfff_flag,
        v.Eth = []; v.Eph = [];
        v.Eth = reshape(Eth_real + 1j.*Eth_imag,[resolution + [1 0] ffmax]);
        v.Eph = reshape(Eph_real + 1j.*Eph_imag,[resolution + [1 0] ffmax]);
        clear Eth_real Eth_imag Eph_real Eph_imag;
    end;

    if field_flag,
        v.Ex_field = []; v.Ey_field = []; v.Ez_field = [];
        v.Ex_field = reshape(Ex_field_real + 1j.*Ex_field_imag,[xmax ymax zmax nfmax]);
        v.Ey_field = reshape(Ey_field_real + 1j.*Ey_field_imag,[xmax ymax zmax nfmax]);
        v.Ez_field = reshape(Ez_field_real + 1j.*Ez_field_imag,[xmax ymax zmax nfmax]);
        clear Ex_field_real Ex_field_imag ...
              Ey_field_real Ey_field_imag ...
              Ez_field_real Ez_field_imag;
    end;
    
    v.time_start = double(v.time_start([1:3 5:7]));
    v.time_end   = double(v.time_end([1:3 5:7]));
    v.partial_times = double(v.partial_times);
    
    %% compute FFTs and impedance
    ComputeResults;
    
end


function ComputeResults

    %% output 'out' is E-field, it needs multiplication with cell span
    %% to get voltage
    v.V = out.*v.Dx;

    %% Number of elements fitting the frequency resolution
    Nt_fft = round( 1/(v.Dt*v.Df) );

    %% extended version of input pulse (Vs) to match the output (V)
    Vs_ext = v.Vs;
    if v.Nt > length(Vs_ext), Vs_ext(v.Nt) = 0; end;


    if v.Rs == 0,
        
        %% source current
        Cb_factor = v.Dt./v.Dx./v.eps0;
        v.I = ( v.V + [ 0 ; Vs_ext(1:end-1) ] )./Cb_factor;

        %% current is shifted half of a time step before start
        v.I = ( v.I + v.I([2:end 1]) )./2;

        %% FFT with given frequency resolution Df is done by zero padding
        %% to length Nt_fft;
        %% 2 and Nt_fft follow the Parseval theorem
        v.Vs_fft = 2.*fft(v.Vs,Nt_fft)./Nt_fft;
        v.I_fft  = 2.*fft(v.I ,Nt_fft)./Nt_fft;
        
        %% taking only the unique part of the spectra
        v.Vs_fft = v.Vs_fft(1:floor(Nt_fft/2)+1);
        v.I_fft  = v.I_fft( 1:floor(Nt_fft/2)+1);
        
        %% with hard source, voltage on load is equivalent to
        %% source voltage
        v.V     = Vs_ext;
        v.V_fft = v.Vs_fft;
        
        %% frequency domain impedance
        v.Z_fft = v.Vs_fft./v.I_fft;
        
    else
        
        %% output voltage is actually in half steps, interpolation needed
        %% (not quite sure about this)
        v.V = ( v.V + [ 0 ; v.V(1:end-1,1) ] )./2;

        %% source current
        v.I = ( v.V + Vs_ext )./v.Rs;
        
        %% FFT with given frequency resolution Df is done by zero padding
        %% to length Nt_fft;
        %% 2 and Nt_fft follow the Parseval theorem
        v.Vs_fft = 2.*fft(v.Vs,Nt_fft)./Nt_fft;
        v.V_fft  = 2.*fft(v.V ,Nt_fft)./Nt_fft;

        %% taking only the unique part of the spectra
        v.Vs_fft = v.Vs_fft(1:floor(Nt_fft/2)+1);
        v.V_fft  = v.V_fft( 1:floor(Nt_fft/2)+1);
        
        %% frequency domain current
        v.I_fft = ( v.V_fft + v.Vs_fft )./v.Rs;
        
        %% frequency domain impedance
        v.Z_fft = -v.Rs./( 1 + v.Vs_fft./v.V_fft );
        
    end;

    %% power entering the domain
    v.P_in = real( abs(v.V_fft).^2 ./ conj(v.Z_fft) )./2;
    
    %% far fields harmonic value
    if v.nfff_flag,

        v.Eth = 2.*v.Eth./Nt_fft;
        v.Eph = 2.*v.Eph./Nt_fft;

        %% theta and phi steps
        Dth =   pi/v.resolution(1);
        Dph = 2*pi/v.resolution(2);

        %% surface matrix
        S = 2.*sin(Dth/2).*Dph.*sin( (0:v.resolution(1)).'.*Dth );
        S([1 end]) = 2.*pi.*(1 - cos(Dth/2))/v.resolution(2); %% cap sector areas at poles
        S = S(:,ones(1,v.resolution(2)));

        %% Poynting vectors in two polarizations
        Pi_theta = ( abs(v.Eth).^2 )./(120*pi)./2;
        Pi_phi   = ( abs(v.Eph).^2 )./(120*pi)./2;
        Pi       = Pi_theta + Pi_phi; %% total

        %% radiated power in far field
        v.P_rad_ff       = zeros(size(v.freq_ff));
        v.P_rad_ff_theta = zeros(size(v.freq_ff));
        v.P_rad_ff_phi   = zeros(size(v.freq_ff));
        for p = 1:length(v.freq_ff),
            v.P_rad_ff(p)       = sum(sum(Pi(:,:,p)      .*S));
            v.P_rad_ff_theta(p) = sum(sum(Pi_theta(:,:,p).*S));
            v.P_rad_ff_phi(p)   = sum(sum(Pi_phi(:,:,p)  .*S));
        end;

    end;
    
    %% near fields harmonic value
    if v.field_flag,
        
        v.Ex_field = 2.*v.Ex_field./Nt_fft;
        v.Ey_field = 2.*v.Ey_field./Nt_fft;
        v.Ez_field = 2.*v.Ez_field./Nt_fft;
        
    end;
    
    %% enabling appropriate menu items
    set([MenuSourceVoltage,MenuInputVoltage,MenuInputCurrent, ...
         MenuProbeResults,MenuImpedance,MenuQandBW,MenuCPUTime],'Enable','on');
    if v.nfff_flag, set([MenuPatt,MenuEfficiency],'Enable','on'); end;
    if ~isempty(v.energy), set(MenuEnergy,'Enable','on'); end;
    if v.field_flag, set([MenuNearFields,MenuCurrents,MenuPostprocessing],'Enable','on'); end;
    
    %% locking the controls
    v.locked_flag = true;
    LockMainFig;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%% RESULTS MENU %%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Source voltage, input voltage and input current display
function ResponseDisplay(src,eventdata,quantity,type)

    %% if not adding, create new figure
    if ~isequal(eventdata,'add'),
        figure('NumberTitle','off');
        plot_axes = axes;
        if src == MenuCompareResults,
            hold all;
        end;
    end;
    
    switch quantity
        
        case 'Source Voltage'
            VariableString = 'Vs';
            VertLabel = 'Voltage [V]';
            
        case 'Input Voltage'
            VariableString = 'V';
            VertLabel = 'Voltage [V]';
            
        case 'Input Current'
            VariableString = 'I';
            VertLabel = 'Current [A]';
            
    end;
    
    switch type
        
        case 'Time'
            
            ydata = v.(VariableString);
            xdata = 1:length(ydata);
            plot(plot_axes,xdata,ydata);
            
            if isequal(eventdata,'add'), return; end;
            
            xlabel('time steps');
            ylabel(VertLabel);
            TitleString = [quantity,' Time-Domain Waveform'];
            title(TitleString);
            set(gcf,'Name',TitleString);
            
        otherwise % in case of 'Magnitude' and 'Phase'
            
            %% Number of elements fitting the frequency resolution
            Nt_fft = round( 1/(v.Dt*v.Df) );

            %% "real" Df
            Df_real = 1/(v.Dt*Nt_fft);
            
            %% frequency axis up to the Nyquist frequency
            f_axis = (0:floor(Nt_fft/2)).*Df_real;

            switch type
                case 'Magnitude'
                    plot(plot_axes,f_axis,abs(v.([VariableString,'_fft'])));
                    if isequal(eventdata,'add'), return; end;
                    ylabel(VertLabel);
                case 'Phase'
                    plot(plot_axes,f_axis,angle(v.([VariableString,'_fft'])));
                    if isequal(eventdata,'add'), return; end;
                    ylabel('phase [rad]');
                    set(gca,'YLim',[-pi pi]);
            end;

            xlabel('f [Hz]');

            TitleString = [quantity,' Frequency-Domain ',type];
            title(TitleString);
            set(gcf,'Name',TitleString);

            %% zoom into the range of interest
            if strcmp(v.ExcitationPulse,'gauss_sine') && ...
               strcmp(v.GaussSineShape,'auto'),
                set(gca,'XLim',[v.fmin v.fmax]);
            else
                %% mask only where the excitation signal is higher than 10 %
                mask = abs(v.Vs_fft) > (0.1*max(abs(v.Vs_fft)));
                fmin = ( find(mask,1,'first') - 1 )*Df_real;
                fmax = ( find(mask,1,'last')  - 1 )*Df_real;
                set(gca,'XLim',[fmin fmax]);
            end;

            h = uimenu(gcf,'Label','Display');
            uimenu(h,'Label','Zoom Out', ...
                     'Callback','axis auto; zoom reset;');

    end;

end



%% Probe results display
function ProbeResults(src,eventdata,type)
    
    j = 0;
    GeomPointer = [];
    %% Probe list
    for i = 1:length(v.Geometry),
        if strcmp(v.Geometry{i}.type,'probe'),
            
            %% adding a space for every newline character in comment
            CommentString = strcat(v.Geometry{i}.comment,{' '});
            CommentString = cat(2,CommentString{:}); %% in one line

            j = j + 1;
            ProbeList{j} = [Components{v.Geometry{i}.component}, ...
                            '(',num2str(v.Geometry{i}.point1(1)),',', ...
                                num2str(v.Geometry{i}.point1(2)),',', ...
                                num2str(v.Geometry{i}.point1(3)),')', ...
                            SeparationString,CommentString];
            GeomPointer(j) = i;
                        
        end;
    end;

    if j == 0,
        errordlg('There are no probes defined', ...
                 'Error','modal');
        uiwait;
        if src == MenuCompareResults,
            error('Error: no probes defined');
        end;
        return;
    end;

    %% if not adding...
    if ~isequal(eventdata,'add'),
        
        %% Probe selection
        [pointer,ok] = listdlg('ListString',ProbeList, ...
                               'SelectionMode','single', ...
                               'ListSize',[160 160], ...
                               'Name','Probes', ...
                               'PromptString','Select probes:');
        if not(ok) || isempty(pointer), return; end;
        
        %% create new figure
        figure('NumberTitle','off');
        plot_axes = axes;
        if src == MenuCompareResults,
            hold all;
            probe_pointer = pointer;
        end;
        
    else
        
        pointer = probe_pointer;
        
        if isempty(GeomPointer) || ( pointer > length(GeomPointer) ),
            error('Error: probes mismatch');
        end;
        
    end;
    
    if v.Geometry{GeomPointer(pointer)}.component < 4, % E-fields
        VertLabel = 'E-field [V/m]';
    else % H-fields
        VertLabel = 'H-field [A/m]';
    end;
    
    switch type
        
        case 'Time'
            
            ydata = v.Probes(:,pointer);
            xdata = 1:length(ydata);
            plot(plot_axes,xdata,ydata);
            
            if isequal(eventdata,'add'), return; end;
            
            xlabel('time steps');
            ylabel(VertLabel);
            TitleString = [ProbeList{pointer},' - Time-Domain Waveform'];
            title(TitleString);
            set(gcf,'Name',TitleString);
            
        otherwise % in case of 'Magnitude' and 'Phase'
            
            %% Number of elements fitting the frequency resolution
            Nt_fft = round( 1/(v.Dt*v.Df) );

            %% "real" Df
            Df_real = 1/(v.Dt*Nt_fft);
            
            %% frequency axis up to the Nyquist frequency
            f_axis = (0:floor(Nt_fft/2)).*Df_real;

            Probe_fft = 2.*fft(v.Probes(:,pointer),Nt_fft)./Nt_fft;
            Probe_fft = Probe_fft(1:floor(Nt_fft/2)+1);
            
            switch type
                case 'Magnitude'
                    plot(plot_axes,f_axis,abs(Probe_fft));
                    if isequal(eventdata,'add'), return; end;
                    ylabel(VertLabel);
                case 'Phase'
                    plot(plot_axes,f_axis,angle(Probe_fft));
                    if isequal(eventdata,'add'), return; end;
                    ylabel('phase [rad]');
                    set(gca,'YLim',[-pi pi]);
            end;

            xlabel('f [Hz]');

            TitleString = [ProbeList{pointer},' - Frequency-Domain ',type];
            title(TitleString);
            set(gcf,'Name',TitleString);

            %% zoom into the range of interest
            if strcmp(v.ExcitationPulse,'gauss_sine') && ...
               strcmp(v.GaussSineShape,'auto'),
                set(gca,'XLim',[v.fmin v.fmax]);
            else
                %% mask only where the excitation signal is higher than 10 %
                mask = abs(v.Vs_fft) > (0.1*max(abs(v.Vs_fft)));
                fmin = ( find(mask,1,'first') - 1 )*Df_real;
                fmax = ( find(mask,1,'last')  - 1 )*Df_real;
                set(gca,'XLim',[fmin fmax]);
            end;

            h = uimenu(gcf,'Label','Display');
            uimenu(h,'Label','Zoom Out', ...
                     'Callback','axis auto; zoom reset;');

    end;
    
end

    

%% Input impedance (and S11) display
function Impedance(src,eventdata,type)
        
    %% if not adding, create new figure
    if ~isequal(eventdata,'add'),
        figure('NumberTitle','off');
        plot_axes = axes;
        if src == MenuCompareResults,
            hold all;
        end;
    end;

    %% Number of elements fitting the frequency resolution
    Nt_fft = round( 1/(v.Dt*v.Df) );
    
    %% "real" Df
    Df_real = 1/(v.Dt*Nt_fft);
    
    %% frequency axis up to the Nyquist frequency
    f_axis = (0:floor(Nt_fft/2)).'.*Df_real;
    
    Z_disp = v.Z_fft; %% local variable -- possible to edit
    Z_disp(1) = 0; %% impedance at zero frequency is usually nonsense

    switch type
        
        case 'real'
            plot(plot_axes,f_axis,real(Z_disp));
            if isequal(eventdata,'add'), return; end;
            ylabel('Re{Z} [Ohm]');
            title('Input Impedance - Real Part');
            set(gcf,'Name','Input Impedance - Real Part');
            
        case 'imag'
            plot(plot_axes,f_axis,imag(Z_disp));
            if isequal(eventdata,'add'), return; end;
            ylabel('Im{Z} [Ohm]');
            title('Input Impedance - Imaginary Part');
            set(gcf,'Name','Input Impedance - Imaginary Part');
            
        case 'mag'
            plot(plot_axes,f_axis,abs(Z_disp));
            if isequal(eventdata,'add'), return; end;
            ylabel('|Z| [Ohm]');
            title('Input Impedance - Magnitude');
            set(gcf,'Name','Input Impedance - Magnitude');
            
        case 'phase'
            plot(plot_axes,f_axis,angle(Z_disp));
            if isequal(eventdata,'add'), return; end;
            ylabel('phase [rad]');
            title('Input Impedance - Phase');
            set(gcf,'Name','Input Impedance - Phase');
            set(gca,'YLim',[-pi pi]);
            
        case 'Q'
            imZ = imag(Z_disp);
            Q = abs( ( imZ([2:end end]) - imZ([1 1:end-1]) ) .* f_axis ...
                ./ ((4*v.Df) .* real(Z_disp)) );
            plot(plot_axes,f_axis,Q);
            if isequal(eventdata,'add'), return; end;
            ylabel('Q');
            title('Quality Factor (Q)');
            set(gcf,'Name','Quality Factor');
            
        case 's11'
            s11 = ( Z_disp - v.Rs )./( Z_disp + v.Rs );
            s11db = 20.*log10(abs(s11));
            plot(plot_axes,f_axis,s11db);
            if isequal(eventdata,'add'), return; end;
            ylabel('|s11| [dB]');
            title('S11 Parameter');
            set(gcf,'Name','S11 Parameter');
            
        case 'vswr'
            s11 = ( Z_disp - v.Rs )./( Z_disp + v.Rs );
            vswr = ( 1 + abs(s11) )./( 1 - abs(s11) );
            plot(plot_axes,f_axis,vswr);
            if isequal(eventdata,'add'), return; end;
            ylabel('VSWR');
            title('Voltage Standing Wave Ratio');
            set(gcf,'Name','Voltage Standing Wave Ratio');
            
    end;
    
    xlabel('f [Hz]');

    %% zoom into the range of interest
    if strcmp(v.ExcitationPulse,'gauss_sine') && ...
       strcmp(v.GaussSineShape,'auto'),
        set(gca,'XLim',[v.fmin v.fmax]);
    else
        %% mask only where the excitation signal is higher than 10 %
        mask = abs(v.Vs_fft) > (0.1*max(abs(v.Vs_fft)));
        fmin = ( find(mask,1,'first') - 1 )*Df_real;
        fmax = ( find(mask,1,'last')  - 1 )*Df_real;
        set(gca,'XLim',[fmin fmax]);
    end;

    h = uimenu(gcf,'Label','Display');
    uimenu(h,'Label','Zoom Out', ...
             'Callback','axis auto; zoom reset;');
    
end



%% Smith chart
function SmithChart(src,eventdata)

    %% if not adding, create new figure
    if ~isequal(eventdata,'add'),
        figure('NumberTitle','off', ...
               'Name','Smith Chart', ...
               'Color','w');
        plot_axes = axes;
        hold on;
    end;

    %% Number of elements fitting the frequency resolution
    Nt_fft = round( 1/(v.Dt*v.Df) );
    
    %% "real" Df
    Df_real = 1/(v.Dt*Nt_fft);

    %% frequency axis up to the Nyquist frequency
    f_axis = (0:floor(Nt_fft/2)).'.*Df_real;
    
    Z_disp = v.Z_fft; %% local variable -- possible to edit
    Z_disp(1) = 0; %% impedance at zero frequency is usually nonsense
    
    %% s11 parameter dimensionless
    s11 = (Z_disp - v.Rs)./(Z_disp + v.Rs);
    
    if strcmp(v.ExcitationPulse,'gauss_sine') && ...
       strcmp(v.GaussSineShape,'auto'),
        %% mask only the frequency range of interest
        mask = ( f_axis > v.fmin ) & ( f_axis < v.fmax );
    else
        %% mask only where the excitation signal is higher than 10 %
        mask = abs(v.Vs_fft) > (0.1*max(abs(v.Vs_fft)));
    end;

    f_axis = f_axis(mask);
    s11 = s11(mask);

    if ~isequal(eventdata,'add'),

        %% drawing of the Smith chart
        axis equal;
        axis([-1 1 -1 1]);
        axis off;

        %% outer circle
        phi = linspace(-pi,pi,201);
        plot(cos(phi)+1j*sin(phi),'k');

        %% real axis
        plot([-1 1],[0 0],'k');

        %% exponential vector from (almost) zero to (almost) infinity
        vec = exp( (-10:.2:10).' );

        %% ( Re(Z) = const ) curves
        for r = [10,25,50,100],
            plot(gam(r+1j*[vec -vec]),'k');
        end;

        %% ( Im(Z) = const ) curves
        for i = [-100,-50,-25,25,50,100],
            plot(gam(vec+1j*i),'k');
        end;

        %% axis labels
        pos = gam( [0,10,25,50,100,100j,50j,25j,-25j,-50j,-100j] );
        h = text(real(pos),imag(pos), ...
            {' 0',' 10',' 25',' 50',' 100',' j100','j50','j25','-j25','-j50','-j100'} );

        set(h(1:8) ,'VerticalAlignment','bottom');
        set(h(9:11),'VerticalAlignment','top');

        set(h([1:6 11]),'HorizontalAlignment','left');
        set(h([7,10])  ,'HorizontalAlignment','center');
        set(h([8, 9])  ,'HorizontalAlignment','right');

        if src == MenuCompareResults,
            hold all;
        end;
        
    end;
    
    %% s11 in the chart -- bold line
    plot(plot_axes,s11,'LineWidth',2);

    function out = gam(in)
        out = (in - 50)./(in + 50);
    end

end



%% Radiation pattern
function Pattern(src,eventdata,type,type2)

    if ~isequal(eventdata,'add'),

        %% select frequency if more than one
        if length(v.freq_ff) == 1,
            p = 1;
        else
            [p,ok] = listdlg('ListString',num2str(v.freq_ff.','%.4g'), ...
                             'SelectionMode','single', ...
                             'ListSize',[160 160], ...
                             'Name','Frequency', ...
                             'PromptString','Select frequency:');
            if not(ok) || isempty(p), return; end;
        end;
        freq_pointer = p;
        
    else
        
        p = freq_pointer;
        if p > length(v.freq_ff),
            error('Error: frequency mismatch');
        end;
    
    end;
    
    TitleFreq = [', f = ',num2str(v.freq_ff(p),'%.4g'),' Hz'];
    
    %% theta and phi steps
    Dth =   pi/v.resolution(1);
    Dph = 2*pi/v.resolution(2);
    
    %% surface matrix
    S = 2.*sin(Dth/2).*Dph.*sin( (0:v.resolution(1)).'.*Dth );
    S([1 end]) = 2.*pi.*(1 - cos(Dth/2))/v.resolution(2); %% cap sector areas at poles
    S = S(:,ones(1,v.resolution(2)));

    %% Poynting vectors in two polarizations
    Pi_theta = ( abs(v.Eth(:,:,p)).^2 )./(120*pi)./2;
    Pi_phi   = ( abs(v.Eph(:,:,p)).^2 )./(120*pi)./2;
    Pi       = Pi_theta + Pi_phi; %% total
    
    if src == MenuCompareResults,
        label = pol_label;
    else
        label = get(get(src,'Parent'),'Label');
    end;
    
    %% Select polarization
    switch label,
        case 'Total',
            Pi_temp = Pi;
            TitlePol = ', total';
        case 'Theta polarized',
            Pi_temp = Pi_theta;
            TitlePol = ', theta polarized';
        case 'Phi polarized',
            Pi_temp = Pi_phi;
            TitlePol = ', phi polarized';
    end;
    
    G  = 4.*pi.*Pi_temp./sum(sum(Pi.*S)); %% total gain

    Gdb = 10.*log10(G); %% gain in dB
    
    switch type
        case 'flat',      TitleString = 'Flat';
        case 'flat_db',   TitleString = 'Flat [dB]';
        case 'sphere',    TitleString = '3D';
        case 'sphere_db', TitleString = '3D [dB]';
        case 'xy',        TitleString = 'X-Y Plane';
        case 'xz',        TitleString = 'X-Z Plane';
        case 'yz',        TitleString = 'Y-Z Plane';
        case 'xy_db',     TitleString = 'X-Y Plane [dB]';
        case 'xz_db',     TitleString = 'X-Z Plane [dB]';
        case 'yz_db',     TitleString = 'Y-Z Plane [dB]';
    end;
    TitleString = ['Radiation Pattern ',TitleString,TitlePol,TitleFreq];
    
    %% if not adding, create new figure
    if ~isequal(eventdata,'add'),
        figure('NumberTitle','off','Name',TitleString);
        plot_axes = axes;
        if src == MenuCompareResults,
            hold all;
        end;
    end;

    switch type
        
        case 'flat'

            imagesc([0,360],[0,180],G(:,[1:end 1]));
            axis image; colorbar;
            xlabel('phi');
            ylabel('theta');
            set(gca,'XTick',0:45:360,'YTick',0:45:180);
            title('Gain');
            
        case 'flat_db'

            %% clipping everything lower than -30dB
            max_Gdb = max(Gdb(:));
            Gdb(Gdb < (max_Gdb-30)) = (max_Gdb-30);

            imagesc([0,360],[0,180],Gdb(:,[1:end 1]),[max_Gdb-30 max_Gdb]);
            axis image; colorbar;
            xlabel('phi');
            ylabel('theta');
            set(gca,'XTick',0:45:360,'YTick',0:45:180);
            title('Gain [dB]');
            
        case 'sphere'
            
            G(:,end+1) = G(:,1); %% closing the sphere
            r = G./max(G(:)); %% radius
            [th,ph] = ndgrid((0:v.resolution(1))*Dth,(0:v.resolution(2))*Dph);
            
            [xx,yy,zz] = sph2cart(ph,pi/2-th,r);
            surf(xx,yy,zz,G);
            
            axis equal;
            axis([-1 1 -1 1 -1 1]);
            axis vis3d;
            colorbar;
            
            xlabel('x');
            ylabel('y');
            zlabel('z');
            title('Gain');

        case 'sphere_db'
            
            %% clipping everything lower than -30dB
            max_Gdb = max(Gdb(:));
            Gdb(Gdb < (max_Gdb-30)) = (max_Gdb-30);
            
            Gdb(:,end+1) = Gdb(:,1); %% closing the sphere
            r = (Gdb - max_Gdb)/30 + 1; %% radius
            [th,ph] = ndgrid((0:v.resolution(1))*Dth,(0:v.resolution(2))*Dph);
            
            [xx,yy,zz] = sph2cart(ph,pi/2-th,r);
            surf(xx,yy,zz,Gdb);
            
            axis equal;
            axis([-1 1 -1 1 -1 1]);
            axis vis3d;
            colorbar;
            
            xlabel('x');
            ylabel('y');
            zlabel('z');
            title('Gain [dB]');
            
        case 'xy'

            PlotPlanePattern(G(round(v.resolution(1)/2)+1,[1:end,1]).');
            
        case 'xz'

            stripe = [ G(end:-1:1,1) ; ...
                       G(2:end-1,round(v.resolution(2)/2)+1) ];
            stripe = circshift(stripe,-round(v.resolution(1)/2));
            
            PlotPlanePattern(stripe([1:end,1]));

        case 'yz'

            stripe = [ G(end:-1:1,round(v.resolution(2)/4)+1) ; ...
                       G(2:end-1, round(v.resolution(2)/4*3)+1) ];
            stripe = circshift(stripe,-round(v.resolution(1)/2));
            
            PlotPlanePattern(stripe([1:end,1]));
            
        case 'xy_db'

            PlotPlanePattern_dB(Gdb(round(v.resolution(1)/2)+1,[1:end,1]).');

        case 'xz_db'

            stripe = [ Gdb(end:-1:1,1) ; ...
                       Gdb(2:end-1,round(v.resolution(2)/2)+1) ];
            stripe = circshift(stripe,-round(v.resolution(1)/2));
            
            PlotPlanePattern_dB(stripe([1:end,1]));

        case 'yz_db'

            stripe = [ Gdb(end:-1:1,round(v.resolution(2)/4)+1) ; ...
                       Gdb(2:end-1, round(v.resolution(2)/4*3)+1) ];
            stripe = circshift(stripe,-round(v.resolution(1)/2));
            
            PlotPlanePattern_dB(stripe([1:end,1]));

    end;
    
    function PlotPlanePattern(y)
        
        len_y = length(y)-1;
        switch type2
            case 'polar'
                x = (0:len_y).'*(2*pi/len_y);
                polar(plot_axes,x,y);
                if isequal(eventdata,'add'), return; end;
            case 'cart'
                x = (0:len_y).'*(360/len_y);
                plot(plot_axes,x,y);
                if isequal(eventdata,'add'), return; end;
                xlim([0,360]);
                set(gca,'XTick',0:45:360);
                xlabel('Angle [deg]');
                ylabel('Gain');
        end;
        title(['Gain',TitlePol,TitleFreq]);
        
    end
    
    function PlotPlanePattern_dB(y)
        
        %% radial limits of the dB values
        max_limit = ceil(max(y(:))/5)*5; %% rounded to nearest multiple of 5
        min_limit = max_limit - 40; %% -40 dB below

        len_y = length(y)-1;
        switch type2
            
            case 'polar'
                
                theta  = (0:len_y).'*(2*pi/len_y);
                rho_db = y;

                delete(plot_axes);
                axes;
                hold on;
                axis equal;
                axis([-1 1 -1.15 1.15]);
                axis off;

                rho_db(rho_db < min_limit) = min_limit;

                %% outer circle and white background
                phi = linspace(-pi,pi,181).';
                plot(cos(phi),sin(phi),'k');
                patch('XData',cos(phi),'YData',sin(phi),'FaceColor','w');

                %% circles
                r = [0.25,0.5,0.75];
                plot(cos(phi)*r,sin(phi)*r,'k:');

                %% tick labels
                r = [0.25;0.5;0.75;1];
                phi = 80/180*pi;
                text(r.*cos(phi),r.*sin(phi),num2str(min_limit+(10:10:40).'), ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','bottom');

                %% radials and angle labels
                for phi_deg = 0:30:330,
                    phi = phi_deg/180*pi;
                    plot([0;cos(phi)],[0;sin(phi)],'k:');
                    text(1.1.*cos(phi),1.1.*sin(phi),num2str(phi_deg), ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle');
                end;

                %% scaling to 'real' rho and plotting
                rho = (rho_db - min_limit)./40;
                plot(rho.*cos(theta),rho.*sin(theta));
        
            case 'cart'
                
                x = (0:len_y).'*(360/len_y);
                plot(plot_axes,x,y);
                if isequal(eventdata,'add'), return; end;
                axis([0,360,min_limit,max_limit]);
                set(gca,'XTick',0:45:360);
                xlabel('Angle [deg]');
                ylabel('Gain [dB]');
                
        end;

        title(['Gain [dB]',TitlePol,TitleFreq]);

    end

end



%% computing SAR distribution in external objects and bricks
function SARDist(src,eventdata)

    %% origin of the FDTD coordinate system
    origin = v.limits(1,:);

    %% lattice dimensions (excluding PML)
    xmax = abs( v.limits(2,1) - v.limits(1,1) );
    ymax = abs( v.limits(2,2) - v.limits(1,2) );
    zmax = abs( v.limits(2,3) - v.limits(1,3) );

    %% total arrays of conductivity and mass density
    sig_array = zeros(xmax,ymax,zmax);
    rho_array = ones(xmax,ymax,zmax); % will be in denominator, hence 1

    for GeomPointer = 1:length(v.Geometry),
        
        switch v.Geometry{GeomPointer}.type
            
            case 'extobj'

                VolArray = v.Geometry{GeomPointer}.VolArray;
                VolSig   = v.Geometry{GeomPointer}.VolSig;
                if isfield(v.Geometry{GeomPointer},'VolRho'),
                    VolRho = v.Geometry{GeomPointer}.VolRho;
                else %% in case that Rho field does not exist, default value
                    VolRho = repmat(DefaultMassDensity,size(VolSig));
                end;
                point1 = v.Geometry{GeomPointer}.point1;

                %% if volume does not exist, next iteration
                if isempty(VolArray), continue; end;

                VolArraySize = size(VolArray);
                vol_limits = [ point1 ; point1+VolArraySize ];
                crop = zeros(2,3);

                %% cropping the object part extending beyond boundaries
                crop(1,:) = v.limits(1,:) - vol_limits(1,:);
                crop(2,:) = vol_limits(2,:) - v.limits(2,:);
                crop(crop<0) = 0;
                VolArray = VolArray( 1+crop(1,1):end-crop(2,1), ...
                                     1+crop(1,2):end-crop(2,2), ...
                                     1+crop(1,3):end-crop(2,3) );
                point1 = point1 + crop(1,:) - origin;
                VolArraySize = size(VolArray);

                %% test for overlapping objects
                mask = false(xmax,ymax,zmax);
                mask(1+point1(1):point1(1)+VolArraySize(1), ...
                     1+point1(2):point1(2)+VolArraySize(2), ...
                     1+point1(3):point1(3)+VolArraySize(3)) = true;
                if any( ( sig_array(mask)~=0 ) & VolArray(:) ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;

                %% filling the arrays
                mask(mask) = VolArray~=0;
                sig_array(mask) = VolSig(VolArray(VolArray~=0));
                rho_array(mask) = VolRho(VolArray(VolArray~=0));
                
            case 'brick'

                %% there is no SAR in PEC or if conductivity is 0 or
                %% if mass density is 0
                if ( v.Geometry{GeomPointer}.PEC ) || ...
                   ( v.Geometry{GeomPointer}.conductivity == 0 ) || ...
                   ( v.Geometry{GeomPointer}.massdensity  == 0 ),
                    continue;
                end;
                
                %% crop the brick into the limits
                point1 = round( v.Geometry{GeomPointer}.point1 );
                point2 = round( v.Geometry{GeomPointer}.point2 );

                brick_limits = [ max( min(point1,point2), v.limits(1,:) ) ; ... 
                                 min( max(point1,point2), v.limits(2,:) ) ] ...;
                               - origin([1 1],:);
                
                %% test for overlapping objects
                mask = false(xmax,ymax,zmax);
                mask( brick_limits(1,1)+1:brick_limits(2,1), ...
                      brick_limits(1,2)+1:brick_limits(2,2), ...
                      brick_limits(1,3)+1:brick_limits(2,3) ) = true;
                  
                if any( sig_array(mask)~=0 ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;
                
                %% filling the arrays
                sig_array(mask) = v.Geometry{GeomPointer}.conductivity;
                rho_array(mask) = v.Geometry{GeomPointer}.massdensity;

        end; %% end switch

    end; %% end for

    %% spatial in-cell averaging of E-fields as they are on the edges
    %% of the voxels (distributed, memory saving)
    v.SAR =         abs( v.Ex_field ...
                       + v.Ex_field(:,[2:end end],:,:) ...
                       + v.Ex_field(:,:,[2:end end],:) ...
                       + v.Ex_field(:,[2:end end],[2:end end],:) ).^2;

    v.SAR = v.SAR + abs( v.Ey_field ...
                       + v.Ey_field([2:end end],:,:,:) ...
                       + v.Ey_field(:,:,[2:end end],:) ...
                       + v.Ey_field([2:end end],:,[2:end end],:) ).^2;

    v.SAR = v.SAR + abs( v.Ez_field ...
                       + v.Ez_field([2:end end],:,:,:) ...
                       + v.Ez_field(:,[2:end end],:,:) ...
                       + v.Ez_field([2:end end],[2:end end],:,:) ).^2;

    %% frequencies of the near fields in Df resolution
    f_samples = round(v.freq_nf./v.Df) + 1;

    %% SAR
    for p = 1:length(v.freq_nf),
        
        v.SAR(:,:,:,p) = v.SAR(:,:,:,p).*sig_array./32./rho_array;
        %% explanation: 32 = 2 (from SAR formula) * 4^2 (from E averaging)
        
        %% normalizing to 1 W
        if v.FieldNormalizing,
            v.SAR(:,:,:,p) = v.SAR(:,:,:,p) ./ v.P_in(f_samples(p));
        end;
        
    end;

    %% enabling menu items for averaging
    set([MenuSARAver1g MenuSARAver10g],'Enable','on');
    set(MenuSARDisp,'Enable','on');

end



%% SAR averaging in external objects and bricks
function SAR_Average(src,eventdata)
    
    %% Geometry object selection
    [GeomPointer,ok] = listdlg('ListString',v.GeometryString, ...
                               'ListSize',[160 160], ...
                               'Name','SAR Averaging', ...
                               'PromptString','Select objects:');
    if not(ok), return; end;
    
    %% Tissue selection in external objects
    sar_button = questdlg('Select tissues in external objects?', ...
                          'SAR Averaging','No');
    switch sar_button
        case 'Yes', sel_tissues = true;
        case 'No',  sel_tissues = false;
        otherwise,  return;
    end;
    
    %% Ask about the limits
    [SARLimits,ok] = limitdlg('SAR Averaging Domain Limits',v.limits);
    if not(ok), return; end;
    
    %% origin of the FDTD coordinate system
    origin = v.limits(1,:);

    %% lattice dimensions (excluding PML)
    xmax = abs( v.limits(2,1) - v.limits(1,1) );
    ymax = abs( v.limits(2,2) - v.limits(1,2) );
    zmax = abs( v.limits(2,3) - v.limits(1,3) );
    
    %% SAR averaging limits in array coordinates
    SARLimits = SARLimits - origin([1 1],:) + [1 1 1;0 0 0];

    %% mass matrix
    mass = zeros(xmax,ymax,zmax);

    for i = GeomPointer,
        
        switch v.Geometry{i}.type
            
            case 'extobj'

                VolArray = v.Geometry{i}.VolArray;
                VolSig   = v.Geometry{i}.VolSig;
                if isfield(v.Geometry{i},'VolRho'),
                    VolRho = v.Geometry{i}.VolRho;
                else
                    VolRho = repmat(DefaultMassDensity,size(VolSig));
                end;
                point1 = v.Geometry{i}.point1;

                %% if volume does not exist, next iteration
                if isempty(VolArray), continue; end;

                VolArraySize = size(VolArray);
                vol_limits = [ point1 ; point1+VolArraySize ];
                crop = zeros(2,3);

                crop(1,:) = v.limits(1,:) - vol_limits(1,:);
                crop(2,:) = vol_limits(2,:) - v.limits(2,:);
                crop(crop<0) = 0;

                VolArray = VolArray( 1+crop(1,1):end-crop(2,1), ...
                                     1+crop(1,2):end-crop(2,2), ...
                                     1+crop(1,3):end-crop(2,3) );
                point1 = point1 + crop(1,:) - origin;
                VolArraySize = size(VolArray);

                %% test for overlapping objects
                mask = false(xmax,ymax,zmax);
                mask(1+point1(1):point1(1)+VolArraySize(1), ...
                     1+point1(2):point1(2)+VolArraySize(2), ...
                     1+point1(3):point1(3)+VolArraySize(3)) = true;

                if any( mass(mask)~=0 & VolArray(:) ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;

                % selecting a particular tissue
                if sel_tissues,
                    tiss_nrs = unique(VolArray(:));
                    tiss_nrs(tiss_nrs==0) = []; % background excluded
                    [tiss_sel,ok] = listdlg( ...
                        'ListString',num2str(tiss_nrs), ...
                        'ListSize',[160 160], ...
                        'Name','SAR Averaging', ...
                        'PromptString',['Select tissues for ', ...
                                        v.GeometryString(i,:),':']);
                    if ~ok, return; end;
                    % zero those which were not selected
                    tiss_nrs(tiss_sel) = [];
                    if ~isempty(tiss_nrs),
                        for ts = tiss_nrs.', VolArray(VolArray==ts) = 0; end;
                    end;
                end;

                %% filling the mass array
                mask(mask) = VolArray~=0;
                mass(mask) = VolRho(VolArray(VolArray~=0)); %% not yet mass...
                
            case 'brick'
            
                %% there is no SAR in PEC or if conductivity is 0 or
                %% if mass density is 0
                if ( v.Geometry{i}.PEC ) || ...
                   ( v.Geometry{i}.conductivity == 0 ) || ...
                   ( v.Geometry{i}.massdensity  == 0 ),
                    continue;
                end;
                
                %% crop the brick into the limits
                point1 = round( v.Geometry{i}.point1 );
                point2 = round( v.Geometry{i}.point2 );

                brick_limits = [ max( min(point1,point2), v.limits(1,:) ) ; ... 
                                 min( max(point1,point2), v.limits(2,:) ) ] ...;
                               - origin([1 1],:);
                
                %% test for overlapping objects
                mask = false(xmax,ymax,zmax);
                mask( brick_limits(1,1)+1:brick_limits(2,1), ...
                      brick_limits(1,2)+1:brick_limits(2,2), ...
                      brick_limits(1,3)+1:brick_limits(2,3) ) = true;
                  
                if any( mass(mask)~=0 ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;
                
                %% filling the mass array
                mass(mask) = v.Geometry{i}.massdensity; %% not yet mass...

            otherwise
            
                h = errordlg('Selected object is not an external object or a brick', ...
                             'Error','modal');
                waitfor(h);
                return;

        end; %% end switch
        
    end; %% end for

    mass = mass.*(v.Dx^3); %% ...now it's mass
    data = mass~=0; %% TRUE where material

    switch src, case MenuSARAver1g,  mass_aver = 1e-3;  %%  1 gram
                case MenuSARAver10g, mass_aver = 10e-3; %% 10 grams
    end;

    %% strings for displaying
    ZmaxString1 = 'slice z = ';
    ZmaxString2 = [' / ',int2str(SARLimits(2,3))];

    disp(' ');
    disp(['SAR ',int2str(mass_aver*1e3),'g averaging started...']);

    %% here it starts
    for p = 1:length(v.freq_nf),
        
        disp(['f = ',num2str(v.freq_nf(p),'%.4g'),' Hz']);
    
        %% mask storage
        F_mask_store = {};
        E_mask_store = {};
        C_mask_store = {};
        largest = 0; %% largest mask created so far

        border = false;

        SAR_weight = v.SAR(:,:,:,p).*mass; %% SAR values pre-weighted by corresponding mass
        SAR_avg = zeros(size(SAR_weight)); %% resulting array

        %% for all voxels in the given limits
        for z = SARLimits(1,3):SARLimits(2,3),
        for y = SARLimits(1,2):SARLimits(2,2),
        for x = SARLimits(1,1):SARLimits(2,1),

            if data(x,y,z), %% if material

                %% volume boundaries
                x1 = x; x2 = x; y1 = y; y2 = y; z1 = z; z2 = z;

                r = 0; %% "radius" of the cube

                fill = mass(x,y,z); %% total mass filling the cube
                fill_prev = 0; %% previous fill

                while fill < mass_aver, %% if the cube is heavier, end

                    r = r + 1;

                    %% expand the cube, but not at the border of the array
                    if x1>1, x1 = x1 - 1; else border = true; end;
                    if y1>1, y1 = y1 - 1; else border = true; end;
                    if z1>1, z1 = z1 - 1; else border = true; end;
                    if x2<xmax, x2 = x2 + 1; else border = true; end;
                    if y2<ymax, y2 = y2 + 1; else border = true; end;
                    if z2<zmax, z2 = z2 + 1; else border = true; end;

                    fill_prev = fill; %% store the previous value of fill (core)

                    mass_cut = mass(x1:x2,y1:y2,z1:z2);
                    fill = sum(mass_cut(:));

                    if fill==fill_prev, break; end; %% island!

                end;

                if fill==fill_prev, continue; end; %% island

                if border, %% if on border, cube is not completely filled

                    sc = 2*r + 1; %% size of the cube

                    r1 = r + 1;

                    %% indices in original array
                    indx = x1:x2;
                    indy = y1:y2;
                    indz = z1:z2;

                    %% indices in the cube
                    indxc = indx + (r1 - x);
                    indyc = indy + (r1 - y);
                    indzc = indz + (r1 - z);

                    %% cube of SAR is partially filled with weighted SAR
                    SAR_cube = zeros(sc,sc,sc);
                    SAR_cube(indxc,indyc,indzc) = SAR_weight(indx,indy,indz);

                    %% cube of mass is partially filled with mass cut
                    mass_cube = zeros(sc,sc,sc);
                    mass_cube(indxc,indyc,indzc) = mass_cut;

                    %% reset for the next cycle
                    border = false;

                else %% ...otherwise filled completely

                    SAR_cube = SAR_weight(x1:x2,y1:y2,z1:z2);
                    mass_cube = mass_cut;

                end;

                if r > largest, %% making the masks, if not done before

                    for i = largest+1:r, %% adding to previously created

                        sc = 2*i + 1; %% size of the cube

                        bd = [1 sc]; %% boundary indices
                        in = 2:sc-1; %% inner part indices

                        %% mask for face voxels
                        F_mask = false(sc,sc,sc);
                        F_mask(bd,in,in) = true;
                        F_mask(in,bd,in) = true;
                        F_mask(in,in,bd) = true;

                        %% mask for edge voxels
                        E_mask = false(sc,sc,sc);
                        E_mask(in,bd,bd) = true;
                        E_mask(bd,in,bd) = true;
                        E_mask(bd,bd,in) = true;

                        %% mask for corner voxels
                        C_mask = false(sc,sc,sc);
                        C_mask(bd,bd,bd) = true;

                        F_mask_store{i} = F_mask;
                        E_mask_store{i} = E_mask;
                        C_mask_store{i} = C_mask;

                    end;

                    largest = r;

                else %% masks are ready

                    F_mask = F_mask_store{r};
                    E_mask = E_mask_store{r};
                    C_mask = C_mask_store{r};

                end;

                %% summing contributions of appropriate masses
                F = sum( mass_cube(F_mask) );
                E = sum( mass_cube(E_mask) );
                C = sum( mass_cube(C_mask) );

                %% obtaining the fraction of the last layer
                if C > 0, %% cubic equation

                    a = E/C;
                    b = F/C;
                    c = (fill_prev - mass_aver)/C;

                    R = (9*a*b - 27*c - 2*a^3)/54;
                    Q = (3*b - a^2)/9;

                    D = R^2 + Q^3; %% discriminant

                    if D > 0,
                        S = ( R + sqrt(D) )^(1/3);
                        if R < sqrt(D),
                            T = ( sqrt(D) - R )^(1/3);
                            f = S - T - a/3;
                        else
                            T = ( R - sqrt(D) )^(1/3);
                            f = S + T - a/3;
                        end;
                    elseif D < 0,
                        S = ( R + 1j*sqrt( -D ) )^(1/3);
                        Re = real( S );
                        f = 2*Re - a/3;
                    else %% D == 0,
                        if R > 0,
                            S = ( R )^(1/3);
                            f = 2*S - a/3;
                        elseif R < 0,
                            S = ( -R )^(1/3);
                            f = S - a/3;
                        else
                            error('Error: R == 0');
                        end;
                    end;

                elseif E > 0, %% quadratic equation

                    a = E;
                    b = F;
                    c = fill_prev - mass_aver;

                    f = ( sqrt( b^2 - 4*a*c ) - b )/(2*a);

                else %% linear equation

                    f = ( mass_aver - fill_prev )/F;

                end;

                %% sometimes the f is out of <0,1> range due to
                %% limited precision
                if f < 0, f = 0; end;
                if f > 1, f = 1; end;

                %% weighting the appropriate parts with fractions
                SAR_cube(F_mask) = SAR_cube(F_mask).*f;
                SAR_cube(E_mask) = SAR_cube(E_mask).*(f^2);
                SAR_cube(C_mask) = SAR_cube(C_mask).*(f^3);

                %% average summing, division by m is done after the loop
                SAR_avg(x,y,z) = sum( SAR_cube(:) );

            end;

        end;
        end;
        disp([ZmaxString1,int2str(z),ZmaxString2]); drawnow;
        end;

        %% division by averaged mass
        SAR_avg = SAR_avg./mass_aver;

        switch src,
            case MenuSARAver1g,
                v.SAR_Aver_1g(:,:,:,p)  = SAR_avg; %% 1 gram averaged SAR
            case MenuSARAver10g
                v.SAR_Aver_10g(:,:,:,p) = SAR_avg; %% 10 gram averaged SAR
        end;
    
    end;

    disp(['SAR ',int2str(mass_aver*1e3),'g averaging finished.']);

    switch src,
        case MenuSARAver1g,  set(MenuSAR1gDisp ,'Enable','on');
        case MenuSARAver10g, set(MenuSAR10gDisp,'Enable','on');
    end;

end



%% computing E-field magnitude distribution
function Efield(src,eventdata)

    %% spatial in-cell averaging (distributed, memory saving)
    v.E_field =             abs( v.Ex_field ...
                               + v.Ex_field(:,[2:end end],:,:) ...
                               + v.Ex_field(:,:,[2:end end],:) ...
                               + v.Ex_field(:,[2:end end],[2:end end],:) ).^2;

    v.E_field = v.E_field + abs( v.Ey_field ...
                               + v.Ey_field([2:end end],:,:,:) ...
                               + v.Ey_field(:,:,[2:end end],:) ...
                               + v.Ey_field([2:end end],:,[2:end end],:) ).^2;
           
    v.E_field = v.E_field + abs( v.Ez_field ...
                               + v.Ez_field([2:end end],:,:,:) ...
                               + v.Ez_field(:,[2:end end],:,:) ...
                               + v.Ez_field([2:end end],[2:end end],:,:) ).^2;
    
    %% E-field magnitude
    v.E_field = sqrt( v.E_field )./4; %% 4 ... from averaging

    %% frequencies of the near fields in Df resolution
    f_samples = round(v.freq_nf./v.Df) + 1;
    
    %% normalizing to 1 W
    if v.FieldNormalizing,
        for p = 1:length(v.freq_nf),
            v.E_field(:,:,:,p) = v.E_field(:,:,:,p)./sqrt(v.P_in(f_samples(p)));
        end;
    end;
    
    set(MenuEmagField,'Enable','on');

end



%% computing H-field magnitude distribution
function Hfield(src,eventdata)
    
    size_E = size(v.Ex_field);

    %% H component fields only for one frequency
    Hx_field = zeros(size_E(1:3));
    Hy_field = zeros(size_E(1:3));
    Hz_field = zeros(size_E(1:3));
    
    %% final H field
    v.H_field = zeros(size_E);

    %% for all NF frequencies
    for p = 1:length(v.freq_nf),
        
        %% constant with phase shift - p.184
        C = v.Dt/(2j*v.mu0*v.Dx*sin(pi*v.freq_nf(p)*v.Dt));

        Hx_field = C.*( v.Ey_field(:,:,[2:end end],p) - v.Ey_field(:,:,:,p) ...
                      - v.Ez_field(:,[2:end end],:,p) + v.Ez_field(:,:,:,p) );
        Hy_field = C.*( v.Ez_field([2:end end],:,:,p) - v.Ez_field(:,:,:,p) ...
                      - v.Ex_field(:,:,[2:end end],p) + v.Ex_field(:,:,:,p) );
        Hz_field = C.*( v.Ex_field(:,[2:end end],:,p) - v.Ex_field(:,:,:,p) ...
                      - v.Ey_field([2:end end],:,:,p) + v.Ey_field(:,:,:,p) );

        %% spatial in-cell averaging
        Hx_field = ( Hx_field + Hx_field([2:end end],:,:) )./2;
        Hy_field = ( Hy_field + Hy_field(:,[2:end end],:) )./2;
        Hz_field = ( Hz_field + Hz_field(:,:,[2:end end]) )./2;

        v.H_field(:,:,:,p) = sqrt( abs(Hx_field).^2 + ...
                                   abs(Hy_field).^2 + ...
                                   abs(Hz_field).^2 );
    
    end;
    
    %% frequencies of the near fields in Df resolution
    f_samples = round(v.freq_nf./v.Df) + 1;
    
    %% normalizing to 1 W
    if v.FieldNormalizing,
        for p = 1:length(v.freq_nf),
            v.H_field(:,:,:,p) = v.H_field(:,:,:,p)./sqrt(v.P_in(f_samples(p)));
        end;
    end;
    
    clear Hx_field Hy_field Hz_field;

    set(MenuHmagField,'Enable','on');
    
end



%% computing current density (J) magnitude distribution
function Jfield(src,eventdata)
        
    %% origin of the FDTD coordinate system
    origin = v.limits(1,:);

    %% lattice dimensions (excluding PML)
    xmax = abs( v.limits(2,1) - v.limits(1,1) );
    ymax = abs( v.limits(2,2) - v.limits(1,2) );
    zmax = abs( v.limits(2,3) - v.limits(1,3) );

    sig_array = zeros(xmax,ymax,zmax);

    for GeomPointer = 1:length(v.Geometry),

        switch v.Geometry{GeomPointer}.type
            
            case 'extobj'

                VolArray = v.Geometry{GeomPointer}.VolArray;
                VolSig   = v.Geometry{GeomPointer}.VolSig;
                point1 = v.Geometry{GeomPointer}.point1;

                %% if volume does not exist, next iteration
                if isempty(VolArray), continue; end;

                VolArraySize = size(VolArray);
                vol_limits = [ point1 ; point1+VolArraySize ];
                crop = zeros(2,3);

                crop(1,:) = v.limits(1,:) - vol_limits(1,:);
                crop(2,:) = vol_limits(2,:) - v.limits(2,:);
                crop(crop<0) = 0;

                VolArray = VolArray( 1+crop(1,1):end-crop(2,1), ...
                                     1+crop(1,2):end-crop(2,2), ...
                                     1+crop(1,3):end-crop(2,3) );
                point1 = point1 + crop(1,:) - origin;
                VolArraySize = size(VolArray);

                %% test for overlapping objects
                mask = false(xmax,ymax,zmax);
                mask(1+point1(1):point1(1)+VolArraySize(1), ...
                     1+point1(2):point1(2)+VolArraySize(2), ...
                     1+point1(3):point1(3)+VolArraySize(3)) = true;

                if any( ( sig_array(mask)~=0 ) & VolArray(:) ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;

                %% filling the arrays
                mask(mask) = VolArray~=0;
                sig_array(mask) = VolSig(VolArray(VolArray~=0));

            case 'brick'

                %% there is no J inside PEC or if conductivity is 0
                if ( v.Geometry{GeomPointer}.PEC ) || ...
                   ( v.Geometry{GeomPointer}.conductivity == 0 ),
                    continue;
                end;
                
                %% crop the brick into the limits
                point1 = round( v.Geometry{GeomPointer}.point1 );
                point2 = round( v.Geometry{GeomPointer}.point2 );

                brick_limits = [ max( min(point1,point2), v.limits(1,:) ) ; ... 
                                 min( max(point1,point2), v.limits(2,:) ) ] ...;
                               - origin([1 1],:);
                
                %% test for overlapping objects
                mask = false(xmax,ymax,zmax);
                mask( brick_limits(1,1)+1:brick_limits(2,1), ...
                      brick_limits(1,2)+1:brick_limits(2,2), ...
                      brick_limits(1,3)+1:brick_limits(2,3) ) = true;
                  
                if any( sig_array(mask)~=0 ),
                    h = errordlg('Objects overlapping','Error','modal');
                    waitfor(h);
                    return;
                end;
                
                %% filling the arrays
                sig_array(mask) = v.Geometry{GeomPointer}.conductivity;
                
        end; %% end switch
        
    end; %% end for

    %% spatial in-cell averaging (distributed, memory saving)
    v.J_field =             abs( v.Ex_field ...
                               + v.Ex_field(:,[2:end end],:,:) ...
                               + v.Ex_field(:,:,[2:end end],:) ...
                               + v.Ex_field(:,[2:end end],[2:end end],:) ).^2;

    v.J_field = v.J_field + abs( v.Ey_field ...
                               + v.Ey_field([2:end end],:,:,:) ...
                               + v.Ey_field(:,:,[2:end end],:) ...
                               + v.Ey_field([2:end end],:,[2:end end],:) ).^2;
           
    v.J_field = v.J_field + abs( v.Ez_field ...
                               + v.Ez_field([2:end end],:,:,:) ...
                               + v.Ez_field(:,[2:end end],:,:) ...
                               + v.Ez_field([2:end end],[2:end end],:,:) ).^2;
    
    %% frequencies of the near fields in Df resolution
    f_samples = round(v.freq_nf./v.Df) + 1;

    %% for all frequencies
    for p = 1:length(v.freq_nf),

        %% current density magnitude
        v.J_field(:,:,:,p) = sqrt( v.J_field(:,:,:,p) ).*sig_array./4;
        %% 4 ... from averaging

        %% normalizing to 1 W
        if v.FieldNormalizing,
            v.J_field(:,:,:,p) = v.J_field(:,:,:,p)./sqrt(v.P_in(f_samples(p)));
        end;

    end;
    
    set(MenuJmagField,'Enable','on');
    
end



%% displaying energy
function EnergyDisplay(src,eventdata)

    %% if not adding, create new figure
    if ~isequal(eventdata,'add'),
        figure('NumberTitle','off', ...
               'Name','Total Energy In Domain');
        plot_axes = axes;
        if src == MenuCompareResults,
            hold all;
        end;
    end;
    
    time_axis = v.e_step:v.e_step:v.Nt;
    plot(plot_axes,time_axis,v.energy(time_axis));
    if isequal(eventdata,'add'), return; end;
    xlabel('time steps');
    ylabel('energy [Ws]');
    title('Total Energy In Domain');
    
end



%% displaying radiation efficiency
function EfficiencyDisplay(src,eventdata)
    
    TableData = cell(11,1+length(v.freq_ff));
    
    TableData(:,1) = {'  Frequency [Hz]:'; ...
                      '  Input power [W]:'; ...
                      '  Radiated power, total [W]:'; ...
                      '  Radiated power, theta polarized [W]:'; ...
                      '  Radiated power, phi polarized [W]:'; ...
                      '  Radiation efficiency, total [ratio]:'; ...
                      '  Radiation efficiency, theta polarized [ratio]:'; ...
                      '  Radiation efficiency, phi polarized [ratio]:'; ...
                      '  Radiation efficiency, total [dB]:'; ...
                      '  Radiation efficiency, theta polarized [dB]:'; ...
                      '  Radiation efficiency, phi polarized [dB]:'};
    
    for p = 1:length(v.freq_ff),
        
        f_sample = round(v.freq_ff(p)./v.Df) + 1;

        %% radiation efficiency
        REff = v.P_rad_ff(p)./v.P_in(f_sample);
        REff_dB = 10*log10(REff);

        TableData{ 1,1+p} = num2str(v.freq_ff(p),'%.4g');
        TableData{ 2,1+p} = num2str(v.P_in(f_sample),'%.4g');
        TableData{ 3,1+p} = num2str(v.P_rad_ff(p),'%.4g');
        TableData{ 6,1+p} = num2str(REff,'%.4g');
        TableData{ 9,1+p} = num2str(REff_dB,'%.4g');

        %% if the polarized powers exist, add the polarizations
        if ~isempty(v.P_rad_ff_theta) && ~isempty(v.P_rad_ff_phi),
            
            REff_theta = v.P_rad_ff_theta(p)./v.P_in(f_sample);
            REff_phi   = v.P_rad_ff_phi(p)  ./v.P_in(f_sample);
            
            REff_theta_dB = 10*log10(REff_theta);
            REff_phi_dB   = 10*log10(REff_phi);
            
            TableData{ 4,1+p} = num2str(v.P_rad_ff_theta(p),'%.4g');
            TableData{ 5,1+p} = num2str(v.P_rad_ff_phi(p),'%.4g');
            
            TableData{ 7,1+p} = num2str(REff_theta,'%.4g');
            TableData{ 8,1+p} = num2str(REff_phi,'%.4g');
            
            TableData{10,1+p} = num2str(REff_theta_dB,'%.4g');
            TableData{11,1+p} = num2str(REff_phi_dB,'%.4g');
            
        else
            
            TableData([4,5,7,8,10,11],1+p) = {'N/A'};
            
        end;

    end;
    
    EffDialog = dialog('Visible','off', ...
        'Name','Radiation Efficiency', ...
        'WindowStyle','modal', ...
        'Units','pixels', ...
        'Position',[ScreenCenter-[320 170] 640 335], ...
        'KeyPressFcn',{@RetEscPressed,@CloseFigure,@CloseFigure});
    
    %% table display
    EffTable = uitable(EffDialog, ...
        'Units','pixels', ...
        'Position',[10 45 620 280], ...
        'BackgroundColor',[1 1 1], ...
        'FontName',pref.Font.FontName, ...
        'FontSize',pref.Font.FontSize, ...
        'FontWeight',pref.Font.FontWeight, ...
        'FontUnits',pref.Font.FontUnits, ...
        'FontAngle',pref.Font.FontAngle, ...
        'RowName',[], ...
        'ColumnName',[], ...
        'ColumnWidth',{230 'auto'}, ...
        'Data',TableData);

    %% buttons
    OKButton = Button(EffDialog,'OK',1,@CloseFigure);

    set(EffDialog,'Visible','on'); %% display figure

end



%% displaying Q at the center frequency and bandwidth based on the Q
%% For reference see Rodney Vaughan and Jorgen Bach Andersen's book 
%% "Channels, Propagation and Antennas for Mobile Communication", section 8.2
function QandBWDisplay(src,eventdata,type)

    %% Number of elements fitting the frequency resolution
    Nt_fft = round( 1/(v.Dt*v.Df) );
    
    %% "real" Df
    Df_real = 1/(v.Dt*Nt_fft);
    
    %% frequency axis up to the Nyquist frequency
    f_axis = (0:floor(Nt_fft/2)).'.*Df_real;
    
    Z_disp = v.Z_fft; %% local variable -- possible to edit
    Z_disp(1) = 0; %% impedance at zero frequency is usually nonsense
           
    %% s11 parameter
    s11 = (Z_disp - v.Rs)./(Z_disp + v.Rs);
    s11db = 20.*log10(abs(s11));

    if strcmp(v.ExcitationPulse,'gauss_sine') && ...
       strcmp(v.GaussSineShape,'auto'),
        %% mask only the frequency range of interest
        mask = ( f_axis > v.fmin ) & ( f_axis < v.fmax );
    else
        %% mask only where the excitation signal is higher than 10 %
        mask = v.Vs_fft > (0.1*max(v.Vs_fft));
    end;
    
    %% zero s11 everywhere else
    s11db(~mask) = 0;
    
    %% resonance where minimum of s11 occurs
    [s11db_min,ind] = min(s11db);
    
    %% Q-factor
    imZ = imag(Z_disp);
    Q = abs( ( imZ([2:end end]) - imZ([1 1:end-1]) ) .* f_axis ...
        ./ ((4*v.Df) .* real(Z_disp)) );
    
    %% warning about not well defined resonance
    if s11db_min > -10, str = 'Warning: resonance not well defined!';
    else                str = '';
    end;
    
    message = { ['Resonance found at frequency: ',num2str(f_axis(ind)),' Hz'], ...
                ['S11 at resonance: ',num2str(s11db_min),' dB'], ...
                str, ...
                [''], ...
                ['Q at resonance: ',num2str(Q(ind))], ...
                [''], ...
                ['Bandwidth at resonance: ',num2str(f_axis(ind)/Q(ind)),' Hz'], ...
                [''], ...
                ['Note: Q and bandwidth are valid only for lossless case.'], ...
                [''] };
    
    msgbox(message,'Q and Bandwidth');

end



%% Displaying the elapsed CPU time
function CPUTime(src,eventdata)
    
    if ~isfield(v,'overall') || isempty(v.overall),
        message = 'CPU time not available for this simulation.';
    else
        message = ['CPU time: ',num2str(v.overall),' sec.'];
    end;
    
    msgbox(message,'CPU Time');
    
end



%% displaying near fields
function ArrayDisplay(src,eventdata)

    %% select frequency if more than one
    if length(v.freq_nf) == 1,
        p = 1;
    else
        [p,ok] = listdlg('ListString',num2str(v.freq_nf.','%.4g'), ...
                         'SelectionMode','single', ...
                         'ListSize',[160 160], ...
                         'Name','Frequency', ...
                         'PromptString','Select frequency:');
        if not(ok) || isempty(p), return; end;
    end;
    TitleFreq = [', f = ',num2str(v.freq_nf(p),'%.4g'),' Hz'];
    
    if v.FieldNormalizing, TitleNorm = ' (normalized to 1W input power)';
    else                   TitleNorm = '';
    end;

    switch src
        
        case MenuExField,
            slicer(v.Ex_field(:,:,:,p),'V/m',['Ex-field',TitleFreq]);
            
        case MenuEyField,
            slicer(v.Ey_field(:,:,:,p),'V/m',['Ey-field',TitleFreq]);
            
        case MenuEzField,
            slicer(v.Ez_field(:,:,:,p),'V/m',['Ez-field',TitleFreq]);

        case MenuHxField,
            % new constant - p.184
            C = v.Dt/(2j*v.mu0*v.Dx*sin(pi*v.freq_nf(p)*v.Dt));
            Hx_field = C.*( v.Ey_field(:,:,[2:end end],p) - v.Ey_field(:,:,:,p) ...
                          - v.Ez_field(:,[2:end end],:,p) + v.Ez_field(:,:,:,p) );
            slicer(Hx_field,'A/m',['Hx-field',TitleFreq]);
            clear Hx_field;
            
        case MenuHyField,
            % new constant - p.184
            C = v.Dt/(2j*v.mu0*v.Dx*sin(pi*v.freq_nf(p)*v.Dt));
            Hy_field = C.*( v.Ez_field([2:end end],:,:,p) - v.Ez_field(:,:,:,p) ...
                          - v.Ex_field(:,:,[2:end end],p) + v.Ex_field(:,:,:,p) );
            slicer(Hy_field,'A/m',['Hy-field',TitleFreq]);
            clear Hy_field;
            
        case MenuHzField,
            % new constant - p.184
            C = v.Dt/(2j*v.mu0*v.Dx*sin(pi*v.freq_nf(p)*v.Dt));
            Hz_field = C.*( v.Ex_field(:,[2:end end],:,p) - v.Ex_field(:,:,:,p) ...
                          - v.Ey_field([2:end end],:,:,p) + v.Ey_field(:,:,:,p) );
            slicer(Hz_field,'A/m',['Hz-field',TitleFreq]);
            clear Hz_field;

        case MenuEmagField,
            slicer(v.E_field(:,:,:,p),'V/m',['E-field Magnitude',TitleNorm,TitleFreq]);
            
        case MenuHmagField,
            slicer(v.H_field(:,:,:,p),'A/m',['H-field Magnitude',TitleNorm,TitleFreq]);
            
        case MenuJmagField,
            slicer(v.J_field(:,:,:,p),'S/m',['J-field Magnitude',TitleNorm,TitleFreq]);

        case MenuSARDisp,
            slicer(v.SAR(:,:,:,p),'W/kg',['SAR Local',TitleNorm,TitleFreq]);
            
        case MenuSAR1gDisp,
            slicer(v.SAR_Aver_1g(:,:,:,p),'W/kg',['SAR 1g-averaged',TitleNorm,TitleFreq]);
            
        case MenuSAR10gDisp,
            slicer(v.SAR_Aver_10g(:,:,:,p),'W/kg',['SAR 10g-averaged',TitleNorm,TitleFreq]);

    end;
    
end



%% displaying wire currents
function WireCurrents(src,eventdata)
    
    %% selecting only wires for the list string
    TypeWires = false(length(v.Geometry),1);
    for i = 1:length(v.Geometry),
        if isequal(v.Geometry{i}.type,'wire'),
            TypeWires(i) = true;
        end;
    end;
    WireString = v.GeometryString(TypeWires,:);
    
    if isempty(WireString),
        msgbox('No wire objects found.','Info','warn','modal');
        uiwait;
        return;
    end;
    
    %% list dialog box for selecting a single wire
    [WirePointer,ok] = listdlg('ListString',WireString, ...
                               'ListSize',[160 160], ...
                               'Name','Wire Currents', ...
                               'PromptString','Select Wire:', ...
                               'SelectionMode','single');
    if not(ok), return; end;
    
    TypeWires = find(TypeWires); % indices into the Geometry
    GeomPointer = TypeWires(WirePointer); % now pointers into the Geometry
    
    if isempty(GeomPointer), return; end;
    
    %% select frequency if more than one
    if length(v.freq_nf) == 1,
        p = 1;
    else
        [p,ok] = listdlg('ListString',num2str(v.freq_nf.','%.4g'), ...
                         'SelectionMode','single', ...
                         'ListSize',[160 160], ...
                         'Name','Frequency', ...
                         'PromptString','Select frequency:');
        if not(ok) || isempty(p), return; end;
    end;
    TitleFreq = [', f = ',num2str(v.freq_nf(p),'%.4g'),' Hz'];
    
    size_E = size(v.Ex_field);
    
    Hx_field = zeros(size_E(1:3));
    Hy_field = zeros(size_E(1:3));
    Hz_field = zeros(size_E(1:3));

    Ix_field = zeros(size_E(1:3));
    Iy_field = zeros(size_E(1:3));
    Iz_field = zeros(size_E(1:3));
    
    % new constant - p.184
    C = v.Dt/(2j*v.mu0*v.Dx*sin(pi*v.freq_nf(p)*v.Dt));

    Hx_field = C.*( v.Ey_field(:,:,[2:end end],p) - v.Ey_field(:,:,:,p) ...
                  - v.Ez_field(:,[2:end end],:,p) + v.Ez_field(:,:,:,p) );
    Hy_field = C.*( v.Ez_field([2:end end],:,:,p) - v.Ez_field(:,:,:,p) ...
                  - v.Ex_field(:,:,[2:end end],p) + v.Ex_field(:,:,:,p) );
    Hz_field = C.*( v.Ex_field(:,[2:end end],:,p) - v.Ex_field(:,:,:,p) ...
                  - v.Ey_field([2:end end],:,:,p) + v.Ey_field(:,:,:,p) );
              
    Ix_field = ( Hz_field - Hz_field(:,[1 1:end-1],:) ...
               - Hy_field + Hy_field(:,:,[1 1:end-1]) ).*v.Dx;
    Iy_field = ( Hx_field - Hx_field(:,:,[1 1:end-1]) ...
               - Hz_field + Hz_field([1 1:end-1],:,:) ).*v.Dx;
    Iz_field = ( Hy_field - Hy_field([1 1:end-1],:,:) ...
               - Hx_field + Hx_field(:,[1 1:end-1],:) ).*v.Dx;
           
    clear Hx_field Hy_field Hz_field;
    
    %% origin of the FDTD coordinate system
    origin = v.limits(1,:);
              
    handle = v.Geometry{GeomPointer}.handle;

    %% retrieving vertices from the axes
    vertices = [ get(handle,'XData').' , ...
                 get(handle,'YData').' , ...
                 get(handle,'ZData').' ];

    %% shifting vertices into FDTD grid coordinates
    vertices = vertices - origin(ones(size(vertices,1),1),:);

    %% lengths of wire straight segments
    diff_vertices = diff(vertices);

    %% total wire length
    wire_length = sum(abs(diff_vertices(:)));

    %% current vector
    Iw = zeros(wire_length,1);
    ind = 1;

    for i = 1:size(vertices,1)-1,
        
        switch find(diff_vertices(i,:)) %% which direction
            
            case 1, %% x-dir
                sub = vertices(i,1):sign(diff_vertices(i,1)):vertices(i+1,1);
                if diff_vertices(i,1)>0,
                    sub(1) = [];
                else
                    sub(end) = [];
                end;
                for j = sub,
                    Iw(ind) = Ix_field(j,vertices(i,2)+1,vertices(i,3)+1) ...
                        *sign(diff_vertices(i,1));
                    ind = ind + 1;
                end;
                
            case 2, %% y-dir
                sub = vertices(i,2):sign(diff_vertices(i,2)):vertices(i+1,2);
                if diff_vertices(i,2)>0,
                    sub(1) = [];
                else
                    sub(end) = [];
                end;
                for j = sub,
                    Iw(ind) = Iy_field(vertices(i,1)+1,j,vertices(i,3)+1) ...
                        *sign(diff_vertices(i,2));
                    ind = ind + 1;
                end;
                
            case 3, %% z-dir
                sub = vertices(i,3):sign(diff_vertices(i,3)):vertices(i+1,3);
                if diff_vertices(i,3)>0,
                    sub(1) = [];
                else
                    sub(end) = [];
                end;
                for j = sub,
                    Iw(ind) = Iz_field(vertices(i,1)+1,vertices(i,2)+1,j) ...
                        *sign(diff_vertices(i,3));
                    ind = ind + 1;
                end;
                
        end;
        
    end;
    
    figure('NumberTitle','off', ...
           'Name',['Wire Currents - Magnitude',TitleFreq]);
    plot(0.5:length(Iw),abs(Iw));
    xlabel('cells');
    ylabel('Total current [A]');
    YLimits = ylim;
    YLimits(1) = 0;
    ylim(YLimits);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%% GEOMETRY DISPLAYS %%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Display wire
function DisplayWire(pointer,wirecolor)

    %% approximate the wire by staircasing
    vertices = StaircaseWire(v.Geometry{pointer}.point1, ...
                             v.Geometry{pointer}.point2);

    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;
                         
    %% display
    v.Geometry{pointer}.handle = ...
        line(vertices(:,1),vertices(:,2),vertices(:,3), ...
            'LineWidth',2, ...
            'Color',wirecolor, ...
            'Visible',visible);
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;

end %% end DisplayWire



%% Display PEC plate
function DisplayPlate(pointer,faces,vertices,patch_count)
    
    %% aproximate the plate by staircasing,
    %% if the f,v is not available from the Edit checking
    if nargin==1,
        [faces,vertices,patch_count] = StaircasePlane( ...
            v.Geometry{pointer}.point1, ...
            v.Geometry{pointer}.point2, ...
            v.Geometry{pointer}.point3, ...
            v.Geometry{pointer}.type);
    end;

    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;

    %% drawing patch object
    v.Geometry{pointer}.handle = patch('Vertices',vertices,'Faces',faces, ...
                                       'FaceColor','y', ...
                                       'EdgeColor','k', ...
                                       'LineStyle',PatchLineStyle, ...
                                       'Visible',visible, ...
                                       'UserData',patch_count);
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;
    
end %% end DisplayPlate



%% Display brick
function DisplayBrick(pointer)

    %% brick corners
    corners = round( [ v.Geometry{pointer}.point1 ; ...
                       v.Geometry{pointer}.point2 ] );
    
    %% extremes of the coordinates = boundaries of the brick
    minima = min( corners );
    maxima = max( corners );
    
    %% auxiliary variables
    x_vec = minima(1):maxima(1);
    y_vec = minima(2):maxima(2);
    z_vec = minima(3):maxima(3);
    
    size_x = maxima(1) - minima(1) + 1;
    size_y = maxima(2) - minima(2) + 1;
    size_z = maxima(3) - minima(3) + 1;
    
    size_xy = size_x*size_y;
    size_yz = size_y*size_z;
    size_zx = size_z*size_x;
    
    %% patch vertices coordinates
    [Y_x,Z_x,X_x] = ndgrid(y_vec,z_vec,x_vec([1,end]));
    [Z_y,X_y,Y_y] = ndgrid(z_vec,x_vec,y_vec([1,end]));
    [X_z,Y_z,Z_z] = ndgrid(x_vec,y_vec,z_vec([1,end]));
  
    vertices = [ X_x(:) Y_x(:) Z_x(:) ; ...
                 X_y(:) Y_y(:) Z_y(:) ; ...
                 X_z(:) Y_z(:) Z_z(:) ];
             
    %% patch faces indices into list of vertices
    ind_x = zeros(size_y,size_z);
    ind_y = zeros(size_z,size_x);
    ind_z = zeros(size_x,size_y);
    
    ind_x(:) = 1:size_yz;
    ind_y(:) = 1:size_zx;
    ind_z(:) = 1:size_xy;

    ind_x = ind_x(1:end-1,1:end-1);
    ind_y = ind_y(1:end-1,1:end-1);
    ind_z = ind_z(1:end-1,1:end-1);
    
    %% x-oriented faces
    ind_x = ind_x(:);
    ind_x(:,4) = ind_x(:,1) + size_y;
    ind_x(:,2) = ind_x(:,1) + 1;
    ind_x(:,3) = ind_x(:,4) + 1;
    
    %% y-oriented faces
    ind_y = ind_y(:) + 2*size_yz;
    ind_y(:,4) = ind_y(:,1) + size_z;
    ind_y(:,2) = ind_y(:,1) + 1;
    ind_y(:,3) = ind_y(:,4) + 1;
    
    %% z-oriented faces
    ind_z = ind_z(:) + 2*(size_yz + size_zx);
    ind_z(:,4) = ind_z(:,1) + size_x;
    ind_z(:,2) = ind_z(:,1) + 1;
    ind_z(:,3) = ind_z(:,4) + 1;
    
    faces = [ ind_x ; ind_x + size_yz ; ...
              ind_y ; ind_y + size_zx ; ...
              ind_z ; ind_z + size_xy ];
          
    %% number of patches in particular directions - simplifies preprocessing
    patch_count = 2.*[ size(ind_x,1) size(ind_y,1) size(ind_z,1) ];
    
    %% reduction of vertices by eliminating the repeated ones
    [vertices,ind1,ind2] = unique(vertices,'rows');
    faces = ind2(faces);

    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;

    %% drawing patch object
    v.Geometry{pointer}.handle = patch('Vertices',vertices, ...
                                       'Faces',faces, ...
                                       'FaceColor','g', ...
                                       'EdgeColor','k', ...
                                       'LineStyle',PatchLineStyle, ...
                                       'Visible',visible, ...
                                       'UserData',patch_count);
    
    %% bricks need to be below plates in stacking
    uistack(v.Geometry{pointer}.handle,'bottom');
    
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;

end %% end DisplayBrick



%% Display source
function DisplaySource(pointer)
        
    %% source line
    vertices = [ v.Geometry{pointer}.gndpoint ; ...
                 v.Geometry{pointer}.point2 ];

    %% visibility
    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;
             
    %% display
    v.Geometry{pointer}.handle = ...
        line(vertices(:,1),vertices(:,2),vertices(:,3), ...
        'LineWidth',2,'Color','r','Visible',visible);
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;
        
end %% end DisplaySource



%% Display field source
function DisplayFieldSource(pointer)
    
    %% array corners
    corners = round( [ v.Geometry{pointer}.point1 ; ...
                       v.Geometry{pointer}.point2 ] );
    
    %% extremes of the coordinates = boundaries of the source array
    minima = min( corners );
    maxima = max( corners );
    
    %% auxiliary variables
    x_vec = minima(1):maxima(1);
    y_vec = minima(2):maxima(2);
    z_vec = minima(3):maxima(3);

    %% coordinate arrays
    [X,Y,Z] = ndgrid(x_vec,y_vec,z_vec);
    
    %% each column represents one line
    X = [X(:)';X(:)'];
    Y = [Y(:)';Y(:)'];
    Z = [Z(:)';Z(:)'];
    
    switch v.Geometry{pointer}.component % H-grid shifted
        case {4,5,6}, X = X + 0.5; Y = Y + 0.5; Z = Z + 0.5;
    end;
    
    switch v.Geometry{pointer}.component % Ex, Ey, Ez, Hx, Hy, Hz
        case 1, X(1,:) = X(1,:) + 0.25; X(2,:) = X(2,:) + 0.75;
        case 2, Y(1,:) = Y(1,:) + 0.25; Y(2,:) = Y(2,:) + 0.75;
        case 3, Z(1,:) = Z(1,:) + 0.25; Z(2,:) = Z(2,:) + 0.75;
        case 4, X(1,:) = X(1,:) - 0.25; X(2,:) = X(2,:) - 0.75;
        case 5, Y(1,:) = Y(1,:) - 0.25; Y(2,:) = Y(2,:) - 0.75;
        case 6, Z(1,:) = Z(1,:) - 0.25; Z(2,:) = Z(2,:) - 0.75;
    end;
    
    switch v.Geometry{pointer}.component % E red, H blue
        case {1,2,3}, field_color = 'r';
        case {4,5,6}, field_color = 'b';
    end;
    
    %% visibility
    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;
             
    %% display
    v.Geometry{pointer}.handle = line(X,Y,Z, ...
        'LineWidth',2,'Color',field_color,'Visible',visible);
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;
    
end %% end DisplayFieldSource



%% Display probe
function DisplayProbe(pointer)
    
    point1 = v.Geometry{pointer}.point1;
    
    switch v.Geometry{pointer}.component
        case 1 % Ex
            vertices = [ point1 ; point1 + [1 0 0] ];
        case 2 % Ey
            vertices = [ point1 ; point1 + [0 1 0] ];
        case 3 % Ez
            vertices = [ point1 ; point1 + [0 0 1] ];
        case 4 % Hx
            vertices = [ point1 + [-0.5 0.5 0.5] ; point1 + 0.5 ];
        case 5 % Hy
            vertices = [ point1 + [0.5 -0.5 0.5] ; point1 + 0.5 ];
        case 6 % Hz
            vertices = [ point1 + [0.5 0.5 -0.5] ; point1 + 0.5 ];
    end;
    
    %% visibility
    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;
             
    %% display
    v.Geometry{pointer}.handle = ...
        line(vertices(:,1),vertices(:,2),vertices(:,3), ...
        'LineWidth',2,'Color','m','Visible',visible);
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;
        
end %% end DisplayProbe



%% Display external object (optimized for speed on large objects)
function DisplayExtObj(pointer)

    %% reading from structure
    point1 = v.Geometry{pointer}.point1;
    VolArray = v.Geometry{pointer}.VolArray;
    VolFaceColor = v.Geometry{pointer}.VolFaceColor;
    VolEdgeColor = v.Geometry{pointer}.VolEdgeColor;
    PlxArray    = v.Geometry{pointer}.PlxArray;
    PlyArray    = v.Geometry{pointer}.PlyArray;
    PlzArray    = v.Geometry{pointer}.PlzArray;
    PlFaceColor = v.Geometry{pointer}.PlFaceColor;
    PlEdgeColor = v.Geometry{pointer}.PlEdgeColor;
    GNDPoint = v.Geometry{pointer}.GNDPoint;
    Direction = v.Geometry{pointer}.Direction;

    %% estimated number of materials
    NumberOfMaterials = max(VolArray(:));
    
    %% visibility
    if v.Geometry{pointer}.visible, visible = 'on';
                               else visible = 'off'; end;
    
    
    
    %% (1) Volume display %%
    
    handle_vol = zeros(NumberOfMaterials,1);

    VolArray(end+1,end+1,end+1) = 0; %% extending domain by 1 in all directions
    data_all = VolArray~=0; %% where some material
    ExtArraySize = size(VolArray); %% extended array dimensions

    %% auxiliary variables
    sz1 = ExtArraySize(1);
    sz2 = ExtArraySize(2);
    sz12 = sz1*sz2;

    %% face matrices: 1s where faces
    face_x = xor(data_all,data_all([end 1:end-1],:,:));
    face_y = xor(data_all,data_all(:,[end 1:end-1],:));
    face_z = xor(data_all,data_all(:,:,[end 1:end-1]));

    %% vertex matrix: 1s where vertices
    vertex = false(ExtArraySize);

    %% indices of faces and vertices
    indx = find(face_x);
    indx(:,4) = indx(:,1) + sz12;
    indx(:,2) = indx(:,1) + sz1;
    indx(:,3) = indx(:,2) + sz12;

    indy = find(face_y);
    indy(:,4) = indy(:,1) + sz12;
    indy(:,2) = indy(:,1) + 1;
    indy(:,3) = indy(:,2) + sz12;

    indz = find(face_z);
    indz(:,4) = indz(:,1) + sz1;
    indz(:,2) = indz(:,1) + 1;
    indz(:,3) = indz(:,2) + sz1;

    %% material data at faces
    data_x = VolArray(indx(:,1) - ~data_all(indx(:,1)));
    data_y = VolArray(indy(:,1) - ~data_all(indy(:,1))*sz1);
    data_z = VolArray(indz(:,1) - ~data_all(indz(:,1))*sz12);

    %% index matrix
    im = zeros(ExtArraySize,'uint32');

    for i = 1:NumberOfMaterials,

        %% mask where material
        mask_x = data_x == i;
        mask_y = data_y == i;
        mask_z = data_z == i;

        %% if the material is not on the surface, try another
        if ~any(mask_x) && ~any(mask_y) && ~any(mask_z),
            continue;
        end;

        %% filling the vertex matrix
        vertex(indx(mask_x,:)) = 1;
        vertex(indy(mask_y,:)) = 1;
        vertex(indz(mask_z,:)) = 1;

        %% linear indices into vertex matrix
        vertex_ind = find(vertex);
        im(vertex_ind) = (1:length(vertex_ind)).';

        %% face definition in terms of PATCH command
        faces = [ im(indx(mask_x,:)) ; ...
                  im(indy(mask_y,:)) ; ...
                  im(indz(mask_z,:)) ];

        %% erasing the arrays for later use (faster when indexed)
        vertex(vertex_ind) = 0;
        im(vertex_ind) = 0;

        %% vertices definition in terms of PATCH command
        ver = zeros(length(vertex_ind),3);
        [ver(:,1),ver(:,2),ver(:,3)] = ind2sub(ExtArraySize,vertex_ind);
        ver = ver - 1; %% actual coordinates start at zero

        %% shifting by position vector
        ver(:,1) = ver(:,1) + point1(1);
        ver(:,2) = ver(:,2) + point1(2);
        ver(:,3) = ver(:,3) + point1(3);

        %% drawing patch
        handle_vol(i) = patch('Vertices',ver,'Faces',faces, ...
                              'FaceColor',VolFaceColor(i,:), ...
                              'EdgeColor',VolEdgeColor{i}, ...
                              'LineStyle',PatchLineStyle, ...
                              'Visible',visible);

    end;

    %% removing blanks where no materials occurred
    handle_vol(handle_vol==0) = [];

    
    %% (2) Planes display %%

    handle_pl = [];

    if ~isempty(PlxArray) || ~isempty(PlyArray) || ~isempty(PlzArray),

        NumberOfMaterials = max( [ max(PlxArray(:)), ...
                                   max(PlyArray(:)), ...
                                   max(PlzArray(:)) ] );
        handle_pl = zeros(NumberOfMaterials,1);

        %% extending domain to the same dimension as VolArray
        if ~isempty(PlxArray), PlxArray(end,end+1,end+1) = 0; end;
        if ~isempty(PlyArray), PlyArray(end+1,end,end+1) = 0; end;
        if ~isempty(PlzArray), PlzArray(end+1,end+1,end) = 0; end;

        %% auxiliary variables
        sz1 = ExtArraySize(1);
        sz2 = ExtArraySize(2);
        sz12 = sz1*sz2;

        %% face matrices: 1s where faces
        if ~isempty(PlxArray),
            face_x = ~and(data_all,data_all([end 1:end-1],:,:)) & PlxArray~=0;
        else
            face_x = false(ExtArraySize);
        end;

        if ~isempty(PlyArray),
            face_y = ~and(data_all,data_all(:,[end 1:end-1],:)) & PlyArray~=0;
        else
            face_y = false(ExtArraySize);
        end;

        if ~isempty(PlzArray),
            face_z = ~and(data_all,data_all(:,:,[end 1:end-1])) & PlzArray~=0;
        else
            face_z = false(ExtArraySize);
        end;

        %% vertex matrix: 1s where vertices
        vertex = false(ExtArraySize);

        %% indices of faces and vertices
        indx = find(face_x);
        indx(:,4) = indx(:,1) + sz12;
        indx(:,2) = indx(:,1) + sz1;
        indx(:,3) = indx(:,2) + sz12;

        indy = find(face_y);
        indy(:,4) = indy(:,1) + sz12;
        indy(:,2) = indy(:,1) + 1;
        indy(:,3) = indy(:,2) + sz12;

        indz = find(face_z);
        indz(:,4) = indz(:,1) + sz1;
        indz(:,2) = indz(:,1) + 1;
        indz(:,3) = indz(:,2) + sz1;

        %% material data at faces
        if ~isempty(PlxArray),
            data_x = PlxArray(indx(:,1));
        else
            data_x = zeros(size(indx,1),1);
        end;

        if ~isempty(PlyArray),
            data_y = PlyArray(indy(:,1));
        else
            data_y = zeros(size(indy,1),1);
        end;

        if ~isempty(PlzArray),
            data_z = PlzArray(indz(:,1));
        else
            data_z = zeros(size(indz,1),1);
        end;

        %% index matrix
        im = zeros(ExtArraySize,'uint32');

        for i = 1:NumberOfMaterials,

            %% mask where material
            mask_x = data_x == i;
            mask_y = data_y == i;
            mask_z = data_z == i;

            %% if the material is not on the surface, try another
            if ~any(mask_x) && ~any(mask_y) && ~any(mask_z),
                continue;
            end;

            %% filling the vertex matrix
            vertex(indx(mask_x,:)) = 1;
            vertex(indy(mask_y,:)) = 1;
            vertex(indz(mask_z,:)) = 1;

            %% linear indices into vertex matrix
            vertex_ind = find(vertex);
            im(vertex_ind) = (1:length(vertex_ind)).';

            %% face definition in terms of PATCH command
            faces = [ im(indx(mask_x,:)) ; ...
                      im(indy(mask_y,:)) ; ...
                      im(indz(mask_z,:)) ];

            %% erasing the arrays for later use (faster when indexed)
            vertex(vertex_ind) = 0;
            im(vertex_ind) = 0;

            %% vertices definition in terms of PATCH command
            ver = zeros(length(vertex_ind),3);
            [ver(:,1),ver(:,2),ver(:,3)] = ind2sub(ExtArraySize,vertex_ind);
            ver = ver - 1; %% actual coordinates start at zero

            %% shifting by position vector
            ver(:,1) = ver(:,1) + point1(1);
            ver(:,2) = ver(:,2) + point1(2);
            ver(:,3) = ver(:,3) + point1(3);

            %% drawing patch
            handle_pl(i) = patch('Vertices',ver,'Faces',faces, ...
                                 'FaceColor',PlFaceColor(i,:), ...
                                 'EdgeColor',PlEdgeColor{i}, ...
                                 'LineStyle',PatchLineStyle, ...
                                 'Visible',visible);

        end;

        %% removing blanks where no materials occurred
        handle_pl(handle_pl==0) = [];

    end;


    %% (3) Source display %%
    
    handle_src = [];
    
    if ~isempty(GNDPoint) && ~isempty(Direction),

        GNDPoint = GNDPoint + point1;

        %% vertices of the source
        point2 = GNDPoint;
        switch Direction
            case 1, point2(1) = point2(1) + 1;
            case 2, point2(1) = point2(1) - 1;
            case 3, point2(2) = point2(2) + 1;
            case 4, point2(2) = point2(2) - 1;
            case 5, point2(3) = point2(3) + 1;
            case 6, point2(3) = point2(3) - 1;
        end;
        vertices = [ GNDPoint ; point2 ];

        %% display source
        handle_src = line(vertices(:,1),vertices(:,2),vertices(:,3), ...
            'LineWidth',2,'Color','r','Visible',visible);

    end;

    %% all together
    v.Geometry{pointer}.handle = [ handle_vol ; handle_pl ; handle_src ];
    if Highlight,
        set(v.Geometry{pointer}.handle,'SelectionHighlight','on');
    else
        set(v.Geometry{pointer}.handle,'SelectionHighlight','off');
    end;

    axis(GeomAxes,v.limits(:));

end %% end DisplayExtObj



%% Staircasing of a wire 
function vertices = StaircaseWire(point1,point2)

    %% new version capable of handling non-aligned endpoints
    %% but should be compatible with the old integer version

    vector = point2 - point1; %% vector pointing from 1 to 2

    round_point1 = round(point1);
    round_point2 = round(point2);
    steps = abs( round_point2 - round_point1 ); %% number of steps in directions
    
    %% if the segment will be parallel to the axis
    if sum(steps==0)>1, vertices = [ round_point1 ; round_point2 ]; return; end;
    
    dir = sign(vector); %% directions relative to the coordinates
    dir( dir==0 ) = 1;
    c = abs( sort([1 2 3].*dir) ); %% order of coordinates to proceed

    offset = dir.*( round_point1 - point1 ); %% ...when points are not in grid
    span = abs( vector );

    %% parameters of the line where transitions between cells occur
    parameters{1} = ( (1:steps(1)).' - 0.5 + offset(1) )./span(1);
    parameters{2} = ( (1:steps(2)).' - 0.5 + offset(2) )./span(2);
    parameters{3} = ( (1:steps(3)).' - 0.5 + offset(3) )./span(3);

    %% coordinates of the corresponding transitions
    coords = [ repmat(c(1),steps(c(1)),1) ; ...
               repmat(c(2),steps(c(2)),1) ; ...
               repmat(c(3),steps(c(3)),1) ];
    
    %% order of the transitions
    [sorted,order] = sort( [ parameters{c(1)} ; ...
                             parameters{c(2)} ; ...
                             parameters{c(3)} ] );
    
    %% coordinates ordered
    coords = coords(order);
    
    %% indices where transitions change direction (a point will come here)
    ind = find( diff(coords)~=0 );
    length_ind = length(ind);
    
    %% ...and picking only corresponding coordinates
    coords = coords(ind);
    
    %% cell steps in particular directions
    steps = diff([0 ; ind]).*dir(coords).';
    
    %% allocation of the vertices
    vertices = zeros(length_ind,3);
    
    %% casting differences between vertices
    vertices( sub2ind([length_ind 3],(1:length_ind).',coords) ) = steps;
    
    %% adding the starting point
    vertices = [ round_point1 ; vertices ];
    
    %% computing the vertices incrementally
    vertices = cumsum(vertices);
    
    %% adding the endpoint
    vertices(end+1,:) = round_point2;
    
end %% end StaircaseWire



%% Staircasing of a plane
function [faces,vertices,patch_count] = StaircasePlane(point1,point2,point3,shape)
    
    if strcmp(shape,'rectangle'),
        
        if isequal(point1,point2), %% singular case
            point3 = point1;
            point4 = point1;
        else
            %% constructing the rectangle in 3D
            p12 = point1 - point2;
            p32 = point3 - point2;
            point3 = point3 - p12.*( (p12*p32.') / (p12*p12.') );
            point4 = point3 + p12;
        end;
        
        %% approximating the edges by staircased wires
        vert = [ StaircaseWire(point1,point2) ; ...
                 StaircaseWire(point2,point3) ; ...
                 StaircaseWire(point3,point4) ; ...
                 StaircaseWire(point4,point1) ];

        points = [ point1 ; point2 ; point3 ; point4 ];

    elseif strcmp(shape,'triangle'),
        
        %% approximating the edges by staircased wires
        vert = [ StaircaseWire(point1,point2) ; ...
                 StaircaseWire(point2,point3) ; ...
                 StaircaseWire(point3,point1) ];
        
        points = [ point1 ; point2 ; point3 ];
        
    else error('Impossible branch');

    end;
        
    %% extremes of the vertex coordinates = volume affected by the triangle
    minima = round( min( points ) );
    maxima = round( max( points ) );

    %% size of the auxiliary array and the masks
    array_size = maxima - minima + 1;
    
    %% projection masks of the patches
    mask_xy = false( array_size([1,2]) );
    mask_yz = false( array_size([2,3]) );
    mask_zx = false( array_size([3,1]) );
    
    %% normalizing the vertex coordinates to subscripts to the aux. array
    vert(:,1) = vert(:,1) - minima(1) + 1;
    vert(:,2) = vert(:,2) - minima(2) + 1;
    vert(:,3) = vert(:,3) - minima(3) + 1;
    
    %% coordinate increments forming the segments of the staircased edges
    diff_vert = diff(vert);

    %% directions of the segments
    dir = sign(diff_vert);
    dir( dir == 0 ) = 1;

    %% fill the three masks with points of all the segments
    for i = 1:size(diff_vert,1),
        
        %% x-oriented segment
        if     diff_vert(i,1),
            ind_x = vert(i,1):dir(i,1):vert(i,1)+diff_vert(i,1);
            mask_xy(ind_x,vert(i,2)) = 1;
            mask_zx(vert(i,3),ind_x) = 1;
            
        %% y-oriented segment
        elseif diff_vert(i,2),
            ind_y = vert(i,2):dir(i,2):vert(i,2)+diff_vert(i,2);
            mask_yz(ind_y,vert(i,3)) = 1;
            mask_xy(vert(i,1),ind_y) = 1;
            
        %% z-oriented segment
        elseif diff_vert(i,3),
            ind_z = vert(i,3):dir(i,3):vert(i,3)+diff_vert(i,3);
            mask_zx(ind_z,vert(i,1)) = 1;
            mask_yz(vert(i,2),ind_z) = 1;
            
        end;
        
    end;

    %% filling the mask inner areas
    mask_xy = cumsum(mask_xy) & flipud(cumsum(flipud(mask_xy)));
    mask_yz = cumsum(mask_yz) & flipud(cumsum(flipud(mask_yz)));
    mask_zx = cumsum(mask_zx) & flipud(cumsum(flipud(mask_zx)));
    
    %% marking the patches where there are four neighboring points in square
    mask_xy = mask_xy(1:end-1,1:end-1) & ...
              mask_xy(1:end-1,2:end  ) & ...
              mask_xy(2:end  ,1:end-1) & ...
              mask_xy(2:end  ,2:end  );
    
    mask_yz = mask_yz(1:end-1,1:end-1) & ...
              mask_yz(1:end-1,2:end  ) & ...
              mask_yz(2:end  ,1:end-1) & ...
              mask_yz(2:end  ,2:end  );
    
    mask_zx = mask_zx(1:end-1,1:end-1) & ...
              mask_zx(1:end-1,2:end  ) & ...
              mask_zx(2:end  ,1:end-1) & ...
              mask_zx(2:end  ,2:end  );
    
    %% if the plane will not generate any patch, return
    if ~any(mask_xy(:)) && ~any(mask_yz(:)) && ~any(mask_zx(:)),
        faces = []; vertices = []; patch_count = [];
        return;
    end;

    %% plane description:
    %% normal vector - decimated to single precision for compatibility
    %% between plates of different point order
    n = cross(point1,point2) + cross(point2,point3) + cross(point3,point1);
    n = double(single(n));
    if n(2)>0, n = -n; end; %% proper order of the patches
    if n(3)>0, n = -n; end; %% for compatibility with wires
    d = n*(point1 - minima + 0.5).'; %% distance from origin - shifted above zero
         
    %% auxiliary coordinate vectors of the auxiliary volume
    x_vec = 1:array_size(1)-1;
    y_vec = 1:array_size(2)-1;
    z_vec = 1:array_size(3)-1;
    
    %% auxiliary variables
    sz1 = array_size(1);
    sz2 = array_size(2);
    sz12 = sz1*sz2;

    %% initialization of indices
    indx = zeros(0,4);
    indy = zeros(0,4);
    indz = zeros(0,4);
    
    %% here, intersections of the lines connecting the centroids of the
    %% cells with the plane are found; rounding is oriented with respect to
    %% orientation of the "corners" - compatible with wires
    
    if n(3)~=0, %% z-faces

        [X,Y] = ndgrid(x_vec,y_vec);

        Z = ceil( ( d - n(1).*X - n(2).*Y )./n(3) );

        %% checking for errors due to limited precision
        error_mask = ( n(1).*X + n(2).*Y + n(3).*Z - d ) > 0;
        Z( error_mask ) = Z( error_mask ) + 1;

        %% indices of patch vertices
        indz = sub2ind(array_size,X(mask_xy),Y(mask_xy),Z(mask_xy));
        indz = indz(:);
        indz(:,4) = indz(:,1) + sz1;
        indz(:,2) = indz(:,1) + 1;
        indz(:,3) = indz(:,2) + sz1;
        
    end;
    
    if n(1)~=0, %% x-faces

        [Y,Z] = ndgrid(y_vec,z_vec);

        if n(1)<0,
            X = ceil(  ( d - n(2).*Y - n(3).*Z )./n(1) );
            error_mask = ( n(1).*X + n(2).*Y + n(3).*Z - d ) > 0;
        else % n(1)>0,
            X = floor( ( d - n(2).*Y - n(3).*Z )./n(1) ) + 1;
            error_mask = ( n(1).*X + n(2).*Y + n(3).*Z - d ) < 0;
        end;
        X( error_mask ) = X( error_mask ) + 1;
        
        %% indices of patch vertices
        indx = sub2ind(array_size,X(mask_yz),Y(mask_yz),Z(mask_yz));
        indx = indx(:);
        indx(:,4) = indx(:,1) + sz12;
        indx(:,2) = indx(:,1) + sz1;
        indx(:,3) = indx(:,2) + sz12;
        
    end;

    if n(2)~=0, %% y-faces

        [Z,X] = ndgrid(z_vec,x_vec);

        if n(2)<0,
            Y = ceil(  ( d - n(3).*Z - n(1).*X )./n(2) );
            error_mask = ( n(1).*X + n(2).*Y + n(3).*Z - d ) > 0;
        else % n(2)>0,
            Y = floor( ( d - n(3).*Z - n(1).*X )./n(2) ) + 1;
            error_mask = ( n(1).*X + n(2).*Y + n(3).*Z - d ) < 0;
        end;
        Y( error_mask ) = Y( error_mask ) + 1;
        
        %% indices of patch vertices
        indy = sub2ind(array_size,X(mask_zx),Y(mask_zx),Z(mask_zx));
        indy = indy(:);
        indy(:,4) = indy(:,1) + sz12;
        indy(:,2) = indy(:,1) + 1;
        indy(:,3) = indy(:,2) + sz12;
        
    end;

    %% index matrix of the patch vertices
    im = zeros(array_size,'uint32');

    %% number of patches in particular directions - simplifies preprocessing
    patch_count = [ size(indx,1), size(indy,1), size(indz,1) ];
    
    ind = [ indx ; indy ; indz ];

    %% linear indices into vertex matrix
    vertex_ind = unique(ind(:));
    im(vertex_ind) = (1:length(vertex_ind)).';

    %% face definition in terms of PATCH command
    faces = im(ind);

    %% vertices definition in terms of PATCH command
    vertices = zeros(length(vertex_ind),3);
    [vertices(:,1),vertices(:,2),vertices(:,3)] = ...
        ind2sub(array_size,vertex_ind);

    %% shifting from normalized to actual position
    vertices(:,1) = vertices(:,1) + minima(1) - 1;
    vertices(:,2) = vertices(:,2) + minima(2) - 1;
    vertices(:,3) = vertices(:,3) + minima(3) - 1;
    
end %% end StaircasePlane



%% Locking the controls when finished simulation - main figure
function LockMainFig;

    %% menus
    set(get(MenuInsert,    'Children'),'Enable','off');
    set(get(MenuSimulation,'Children'),'Enable','off');
    
    %% uicontrols
    h = findobj(MainFig,'Type','uicontrol','-not','Style','pushbutton');
    set(h,'Enable','inactive');
    set([GeomPushInsert,GeomPushDelete],'Enable','off');
    set(GeomListBox,'Enable','on');
    set([GeomMenuInsert,GeomMenuDuplicate,GeomMenuTranslate, ...
         GeomMenuDelete,GeomMenuMoveUp,GeomMenuMoveDown],'Enable','off');
    
    LockImageAxes = axes( ...
        'Parent',MainFig, ...
        'DataAspectRatio',[1 1 1], ...
        'Units','pixels', ...
        'Position',[FigWidth-30 FigHeight-15 15 15], ...
        'Visible','off', ...
        'YDir','reverse', ...
        'XLim',[0.5 15.5], ...
        'YLim',[0.5 15.5]);
    LockImage = image( ...
        'Parent',LockImageAxes, ...
        'CData',LockImageData);

end;



%% Locking particular geometry dialog windows
function LockDialog(fig,ok_button,cancel_button);
    
    h = findobj(fig,'Type','uicontrol');
    set(h,'Enable','inactive');
    set(ok_button,'Enable','off');
    set(cancel_button,'Enable','on');
    
end;



%% Download file with PuTTY
function PuttyDownload(source_file,target_dir)
    
    %% make sure that the remote path is in unix format (with slashes)
    source_file = strrep(source_file,'\','/');

    %% call PSCP
    [null,null] = system(['"', ...
            fullfile(Engine.PuTTyPath,'pscp.exe'), ...
            '" -q -C -l ',Engine.Username, ...
                  ' -pw ',Engine.Password, ...
            ' ',Engine.HostName, ...
            ':"',source_file,'" "',target_dir,'"']);

end

%% Upload file with PuTTY
function PuttyUpload(source_file,target_dir)

    %% make sure that the remote path is in unix format (with slashes)
    target_dir = strrep(target_dir,'\','/');

    %% call PSCP
    system(['"', ...
            fullfile(Engine.PuTTyPath,'pscp.exe'), ...
            '" -q -C -l ',Engine.Username, ...
                  ' -pw ',Engine.Password, ...
            ' "',source_file,'" ', ...
            Engine.HostName,':"',target_dir,'"']);

end

%% Run command with PuTTY
function PuttyRun(command)

    %% make sure that the remote path is in unix format (with slashes)
    command = strrep(command,'\','/');

    %% call PLINK
    system(['"', ...
            fullfile(Engine.PuTTyPath,'plink.exe'), ...
            '" -X -ssh -C -l ',Engine.Username, ...
                       ' -pw ',Engine.Password, ...
            ' ',Engine.HostName,' "',command,'"']);

end



%% compute parameters of the Gauss-Sine pulse
function [fc,fs,offset] = ComputeGaussSinePars(fmin,fmax,td_thres,fd_thres)
    
    fc = ( fmax + fmin )/2;

    fs = ( fmax - fmin ) * pi * sqrt(5/log(10)/-fd_thres);

    offset = 1/(2*fc) * ...
        round( sqrt( fc^2*(-td_thres)*log(10)/5/fs^2 - 0.25 ) );

end



%% definition of the Gauss-Sine pulse
function out = GaussSinePulse(mag,fc,fs,offset)
    
    %% offset of the center of the pulse in time steps
    offset_Dt = round( offset/v.Dt );
    
    %% span of the Gauss-Sine pulse in time steps
    n = (-offset_Dt:offset_Dt).';
    
    %% Gauss-weighted Sine
    out = mag .* sin((2*pi*fc*v.Dt).*n) .* exp( -((fs*v.Dt).*n).^2 );
    
end



%% setting of SAR averaging limits (but can be used for other limits)
function [new_limits,ok] = limitdlg(title,def_limits)

    %% default answer
    new_limits = [];
    ok = 0;

    %% open window
    LimitsDialog = dialog('Visible','off', ...
        'Name',title, ...
        'Units','pixels', ...
        'Position',[(ScreenSize(3:4)-[268 420])./2 268 170], ...
        'KeyPressFcn',{@RetEscPressed,@LimitsCallback,@CloseFigure});

    LimitsPanel = uipanel('Parent',LimitsDialog, ...
                          'ForegroundColor','b', ...
                          'Title','Limits', ...
                          'Units','pixels', ...
                          'Position',[7 45 254 120]);

    uitext(LimitsPanel,'x',[ 50 86 60 16],'HorizontalAlignment','center');
    uitext(LimitsPanel,'y',[115 86 60 16],'HorizontalAlignment','center');
    uitext(LimitsPanel,'z',[180 86 60 16],'HorizontalAlignment','center');

    uitext(LimitsPanel,'lower:',[10 65 40 17]);
    uitext(LimitsPanel,'upper:',[10 40 40 17]);

    xlim1 = uiedit(LimitsPanel,int2str(def_limits(1,1)),[ 50 65 60 21]);
    ylim1 = uiedit(LimitsPanel,int2str(def_limits(1,2)),[115 65 60 21]);
    zlim1 = uiedit(LimitsPanel,int2str(def_limits(1,3)),[180 65 60 21]);
    xlim2 = uiedit(LimitsPanel,int2str(def_limits(2,1)),[ 50 40 60 21]);
    ylim2 = uiedit(LimitsPanel,int2str(def_limits(2,2)),[115 40 60 21]);
    zlim2 = uiedit(LimitsPanel,int2str(def_limits(2,3)),[180 40 60 21]);

    %% edit box handles in the same shape as the corresponding limits
    Editboxes = [ xlim1 ylim1 zlim1 ; xlim2 ylim2 zlim2 ];

    %% pushbutton resets the limits to defaults
    ResetLimitsPush = uipushbutton(LimitsPanel,'Reset',[165 10 75 23], ...
                                   @ResetLimitsCallback);
    function ResetLimitsCallback(src,eventdata)
        for i = 1:6,
            set(Editboxes(i),'String',int2str(def_limits(i)));
        end;
    end

    %% buttons
    OKButton     = Button(LimitsDialog,'OK'    ,2,@LimitsCallback);
    CancelButton = Button(LimitsDialog,'Cancel',1,@CloseFigure);

    %% display figure
    set(LimitsDialog,'Visible','on');

    %% wait until closed
    uiwait(LimitsDialog);

    function LimitsCallback(src,eventdata)

        %% temporary limits
        temp_limits = zeros(size(def_limits));

        %% check the values in the edit boxes for correctness
        for i = 1:6,

            value = str2num( get(Editboxes(i),'String') );
            if isempty(value) || ~isfinite(value) || ...
               imag(value)~=0 || value~=round(value),
                errordlg('Limits must be integer numeric values', ...
                         'Bad input','modal');
                pause(1);
                uicontrol(Editboxes(i)); %% give focus to the field
                return;
            end;
            temp_limits(i) = value;

        end;

        %% check if the upper limits are higher than the lower limits
        if any( temp_limits(1,:) >= temp_limits(2,:) ),
            errordlg('Upper limit must be higher than lower limit', ...
                     'Bad input','modal');
            pause(1);
            uicontrol(Editboxes(1)); %% give focus to the first field
            return;
        end;

        %% set the new limits and ok
        new_limits = temp_limits;
        ok = 1;

        %% close the window
        close(LimitsDialog);

    end

end



%% button for dialogs
%% order ... number of button from the right
function handle = Button(parent,string,order,callback);
    old_units = get(parent,'Units');
    set(parent,'Units','pixels');
    FigurePos = get(parent,'Position');
    handle = uicontrol(parent,'Style','pushbutton', ...
                              'String',string, ...
                              'Units','pixels', ...
                              'Position',[FigurePos(3)-order*84 12 75 23], ...
                              'Callback',callback, ...
                              'KeyPressFcn',@ReturnPressed);
    set(parent,'Units',old_units);
end

%% object callback when Return pressed
function ReturnPressed(src,eventdata)
    if strcmp(eventdata.Key,'return'), 
        feval(get(src,'Callback'));
    end;
end

%% window callback when return or escape pressed
function RetEscPressed(src,eventdata,return_cb,escape_cb)
    if strcmp(eventdata.Key,'return'),
        drawnow;
        feval(return_cb);
    elseif strcmp(eventdata.Key,'escape'),
        feval(escape_cb);
    end;
end

%% callback for closing a figure
function CloseFigure(src,eventdata); close(gcbf); end

%% moving the focus from multiline editboxes on Tab and Shift-Tab key
function TabFocus(src,eventdata);
    if strcmp(eventdata.Key,'tab'),
        siblings = get(get(src,'Parent'),'Children');
        focus = siblings==src;
        if strcmp(eventdata.Modifier,'shift'),
            uicontrol(siblings(circshift(focus(:),1))); %% move back
        else
            uicontrol(siblings(circshift(focus(:),-1))); %% move forward
            drawnow;
            %% remove unintentional Tab character
            set(src,'String',strrep(cellstr(get(src,'String')),char(9),''));
        end;
    end;
end

%% check the editbox for valid numeric value
function status = CheckEditboxValue(handle,sign)
    
    value = str2num( get(handle,'String') );

    if isempty(value) || ~isfinite(value) || imag(value) ~= 0,
        errordlg('Input value must be real number','Bad input','modal');
    elseif strcmp(sign,'negative') && value >= 0,
        errordlg('Input value must be negative number','Bad input','modal');
    elseif isempty(sign) && value <= 0,
        errordlg('Input value must be positive number','Bad input','modal');
    else
        status = 1; %% status 1 if OK
        return;
    end;

    uiwait;
    uicontrol(handle); %% bring the error into focus
    status = 0; %% status 0 if error
    
end

%% callback for the [...] pushbutton :: directory
function GetDir(src,eventdata,EditHandle,TitleString)
    
    PathName = uigetdir(get(EditHandle,'String'),TitleString);
    if isequal(PathName,0), return; end;
    set(EditHandle,'String',PathName);
    
end

%% callback for the [...] pushbutton :: file
function GetFile(src,eventdata,EditHandle,TitleString)

    [FileName,PathName] = uigetfile({'*','All Files (*.*)'}, ...
        TitleString,[WorkDir,filesep]);
    if isequal(FileName,0) || isequal(PathName,0), return; end;
    set(EditHandle,'String',fullfile(PathName,FileName));

end



%% simplified uicontrols follow

%% uitext = uicontrol('Style','text')
function handle = uitext(parent,string,position,varargin)
    handle = uicontrol(parent,'Style','text', ...
                              'String',string, ...
                              'HorizontalAlignment','left', ...
                              'Units','pixels', ...
                              'Position',position, ...
                              varargin{:});
end

%% uiedit = uicontrol('Style','edit')
function handle = uiedit(parent,string,position,varargin)
    handle = uicontrol(parent,'Style','edit', ...
                              'String',string, ...
                              'HorizontalAlignment','right', ...
                              'BackgroundColor','w', ...
                              'Units','pixels', ...
                              'Position',position, ...
                              'Max',1,'Min',0, ...
                              varargin{:});
end

%% uimultiedit = uicontrol('Style','edit','Max',2,'Min',0)
function handle = uimultiedit(parent,string,position,varargin)
    handle = uicontrol(parent,'Style','edit', ...
                              'String',string, ...
                              'HorizontalAlignment','left', ...
                              'BackgroundColor','w', ...
                              'Units','pixels', ...
                              'Position',position, ...
                              'Max',2,'Min',0, ...
                              varargin{:});
end

%% uipopupmenu = uicontrol('Style','popupmenu')
function handle = uipopupmenu(parent,string,value,position,varargin)
    handle = uicontrol(parent,'Style','popupmenu', ...
                              'String',string, ...
                              'Value',value, ...
                              'HorizontalAlignment','left', ...
                              'BackgroundColor','w', ...
                              'Units','pixels', ...
                              'Position',position, ...
                              varargin{:});
end

%% uiradiobutton = uicontrol('Style','radiobutton')
function handle = uiradiobutton(parent,value,position,varargin)
    handle = uicontrol(parent,'Style','radiobutton', ...
                              'Max',1,'Min',0, ...
                              'Value',value, ...
                              'Units','pixels', ...
                              'Position',position, ...
                              varargin{:});
end

%% uicheckbox = uicontrol('Style','checkbox')
function handle = uicheckbox(parent,string,value,position,callback,varargin)
    handle = uicontrol(parent,'Style','checkbox', ...
                              'String',string, ...
                              'Max',1,'Min',0, ...
                              'Value',value, ...
                              'Units','pixels', ...
                              'Position',position, ...
                              'Callback',callback, ...
                              varargin{:});
end

%% uipushbutton = uicontrol('Style','pushbutton')
function handle = uipushbutton(parent,string,position,callback,varargin)
    handle = uicontrol(parent,'Style','pushbutton', ...
                              'String',string, ...
                              'Units','pixels', ...
                              'Position',position, ...
                              'Callback',callback, ...
                              varargin{:});
end

%% uitogglebutton = uicontrol('Style','togglebutton')
function handle = uitogglebutton(parent,string,value,position,callback,varargin)
    handle = uicontrol(parent,'Style','togglebutton', ...
                              'String',string, ...
                              'Max',1,'Min',0, ...
                              'Value',value, ...
                              'Units','pixels', ...
                              'Position',position, ...
                              'Callback',callback, ...
                              varargin{:});
end

%% uiradiogroup = uibuttongroup
function handle = uiradiogroup(parent,position,callback,varargin)
    handle = uibuttongroup('Parent',parent, ...
                           'BorderType','none', ...
                           'Units','pixels', ...
                           'Position',position, ...
                           'SelectionChangeFcn',callback, ...
                           varargin{:});
end

%% simple horizontal or vertical line in GUI -- using uipanel
function handle = uiline(parent,color,position,varargin)
    handle = uipanel('Parent',parent, ...
                     'BorderType','line', ...
                     'HighlightColor',color, ...
                     'Units','pixels', ...
                     'Position',position, ...
                     varargin{:});
end

end % end afc
