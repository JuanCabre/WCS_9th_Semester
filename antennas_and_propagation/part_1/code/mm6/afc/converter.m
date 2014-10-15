function converter
% BDF file converter for AFC
% converts BDF files with old features into new AAF files

% (c) 2007-- Ondrej Franek
% $Revision v1.06$ $06-01-2011$



if isunix, %% nicer open dialog
    DefaultDialogs = getappdata(0,'UseNativeSystemDialogs');
    setappdata(0,'UseNativeSystemDialogs',0);
end;



%% Specify files to convert
[FileNames,PathName] = uigetfile( ...
    {'*.bdf','Binary Definition Files (*.bdf)'}, ...
    'Specify Files to Convert', ...
    'Multiselect','on');

if isequal(FileNames,0) || isequal(PathName,0), return; end;

if ~iscellstr(FileNames), FileNames = { FileNames }; end;



%% Specify whether to add results from .RES files
button = questdlg('Add results from .RES files?','Results','Cancel');
switch button
    case 'Yes'
        ResPathName = uigetdir(PathName,'Select folder with .RES files');
        if isequal(ResPathName,0), return; end;
    case 'No'
        ResPathName = [];
    case 'Cancel'
        return;
end;



SuccessCount = 0;
for i = 1:length(FileNames),

    try
    
        %% Open
        v = load('-mat',fullfile(PathName,FileNames{i}));

        %% Correct compatibility issues

        %% domain
        if all(isfield(v,{'corner1','corner2'})),
            v.limits = [ min(v.corner1,v.corner2) ; max(v.corner1,v.corner2) ];
            v = rmfield(v,{'corner1','corner2'});
        end;
        
        %% FDTD
        if ~isfield(v,'Drop'), v.Drop = 1; end;
        if ~isfield(v,'Drop_dB'), v.Drop_dB = 10.*log10(v.Drop/100); end %% conv. to dB
        v = rmfield(v,'Drop');
        if ~isfield(v,'drop_flag'), v.drop_flag = false; end;
        if  isfield(v,'ParamDepSwitch'), v = rmfield(v,'ParamDepSwitch'); end;
        if  isfield(v,'f'), v.fc = v.f; v = rmfield(v,'f'); end;
        if ~isfield(v,'offset'), v.offset = 370; end
        if  isfield(v,'delay'), v.offset = v.delay*v.Dt;
                                v = rmfield(v,'delay'); end;
        if  isfield(v,'t0'),    v.offset = v.t0;
                                v = rmfield(v,'t0'); end;
        if ~isfield(v,'Dt_mode'), v.Dt_mode = 'manual'; end;
        if ~isfield(v,'CFL_coef'), v.CFL_coef = 0.999; end;
        if ~isfield(v,'mag'), v.mag = 1; end;
        if ~isfield(v,'fmin'), v.fmin =  600e6; end;
        if ~isfield(v,'fmax'), v.fmax = 1200e6; end;
        if ~isfield(v,'fd_thres'), v.fd_thres = -20; end;
        if ~isfield(v,'td_thres'), v.td_thres = -60; end;
        if ~isfield(v,'GaussSineShape'), v.GaussSineShape = 'manual'; end;
        if ~isfield(v,'ExcitationPulse'), v.ExcitationPulse = 'gauss_sine'; end;
        
        %% PML
        if ~isfield(v,'m_s'), v.m_s = 4; end;
        if  isfield(v,'m'), v.m_s = v.m; v = rmfield(v,'m'); end;
        if ~isfield(v,'m_k'), v.m_k = v.m_s; end;
        if ~isfield(v,'m_a'), v.m_a = 1; end;
        if ~isfield(v,'kappa_max'), v.kappa_max = 15; end;
        if ~isfield(v,'a_max'), v.a_max = 0.2; end;
        
        %% Misc
        if ~isfield(v,'field_flag'), v.field_flag = true; end;
        if ~isfield(v,'nfff_flag'), v.nfff_flag = true; end;
        if ~isfield(v,'energy_flag'), v.energy_flag = true; end;
        if ~isfield(v,'e_step'), v.e_step = 100; end;
        if ~isfield(v,'locked_flag'), v.locked_flag = false; end;
        
        %% nf
        if ~isfield(v,'freq_nf'), v.freq_nf = 900e6; end;
        if  isfield(v,'fc'), v.freq_nf = v.fc; end;
        
        %% ff
        if ~isfield(v,'freq_ff'), v.freq_ff = 900e6; end;
        if  isfield(v,'fc'), v.freq_ff = v.fc; end;
        if ~isfield(v,'x_nfff'), v.x_nfff = [1 1]; end;
        if ~isfield(v,'y_nfff'), v.y_nfff = [1 1]; end;
        if ~isfield(v,'z_nfff'), v.z_nfff = [1 1]; end;
        if ~isfield(v,'phcenter'), v.phcenter = [0 0 0]; end;
        if ~isfield(v,'resolution'), v.resolution = [90 180]; end;
        
        %% Geometry
        for GeomPointer = 1:length(v.Geometry),
            
            e = v.Geometry{GeomPointer};
            
            if strcmp(e.type,'volobj'),
                e.type = 'extobj';
                if ~iscell(e.VolEdgeColor),
                    e.VolEdgeColor = num2cell(e.VolEdgeColor,2);
                end;
                if ~iscell(e.PlEdgeColor),
                    e.PlEdgeColor = num2cell(e.PlEdgeColor,2);
                end;
                v.Geometry{GeomPointer} = e;
                v.GeometryString(GeomPointer,:) = ...
                    strrep(v.GeometryString(GeomPointer,:),'volobj','extobj');
                
            elseif strcmp(e.type,'brick') && ~isfield(e,'PEC'),
                e.point2 = e.point1 + e.point2;
                if isinf(e.permittivity), e.PEC = true;
                                          e.permittivity = 1;
                                          e.conductivity = 0;
                else e.PEC = false;
                end;
                v.Geometry{GeomPointer} = e;
                
            elseif strcmp(e.type,'plane'),
                e.type = 'rectangle';
                
                mod   = e.Extension;
                alpha = e.Angle*pi/180;
                
                %% coincident points of a base
                if e.point1==e.point2, e.point2(2) = e.point1(2) + 1; end;
                if e.point1(2)==e.point2(2), %% x-oriented base
                    e.point3 = e.point1 + [ 0 mod*cos(alpha) mod*sin(alpha) ];
                else %% otherwise
                    e.point3 = e.point1 + [ mod*cos(alpha) 0 mod*sin(alpha) ];
                end;
                e = rmfield(e,{'Extension','Angle'});
                
                %% replacing string in the geometry list
                k = strfind(v.GeometryString(GeomPointer,:),'plane');
                v.GeometryString(GeomPointer,k(1)+9:end+4) = ...
                    v.GeometryString(GeomPointer,k(1)+5:end); %% shifting comment
                v.GeometryString(GeomPointer,k(1):k(1)+8) = 'rectangle';
                
                v.Geometry{GeomPointer} = e;
                
            elseif strcmp(e.type,'resistor') && isfield(e,'Value'),
                e.resistance = e.Value;
                e = rmfield(e,'Value');
                v.Geometry{GeomPointer} = e;

            elseif strcmp(e.type,'capacitor') && isfield(e,'Value'),
                e.capacitance = e.Value;
                e = rmfield(e,'Value');
                v.Geometry{GeomPointer} = e;

            end;
            
            if ~isfield(e,'visible'),
                e.visible = true;
                v.Geometry{GeomPointer} = e;
            end;
            
        end;
        
        v.locked_flag = false;
        
        offset_Dt = round( v.offset/v.Dt );
        n = (-offset_Dt:offset_Dt).';
        v.Vs = v.mag .* sin((2*pi*v.fc*v.Dt).*n) .* exp( -((v.fs*v.Dt).*n).^2 );
        
        %% default result values - empty
        v.V = []; v.I = [];
        v.Vs_fft = []; v.V_fft = []; v.I_fft = []; v.Z_fft = [];
        v.P_in = []; v.P_rad_nf = []; v.P_rad_ff = [];
        v.Eth = []; v.Eph = [];
        v.Ex_field = []; v.Ey_field = []; v.Ez_field = [];
        v.energy = []; v.SAR = []; v.SAR_Aver_1g = []; v.SAR_Aver_10g = [];
        v.E_field = []; v.H_field = []; v.J_field = [];
        
        [null,FileNameBody,null,null] = fileparts(FileNames{i});
        
        if ~isempty(ResPathName),
            r = load('-mat',fullfile(ResPathName,[FileNameBody '.res']));
            VarNames = fieldnames(r);
            for j = 1:length(VarNames),
                if strcmp(VarNames{j},'out'), continue; end;
                if strcmp(VarNames{j},'Z'), v.Z_fft = r.Z; continue; end;
                v.(VarNames{j}) = r.(VarNames{j});
            end;
            clear r;
            
            Nt_fft = round( 1/(v.Dt*v.Df) );
            v.Vs_fft = 2.*fft(v.Vs,Nt_fft)./Nt_fft;
            v.V_fft  = 2.*fft(v.V ,Nt_fft)./Nt_fft;
            v.I_fft  = 2.*fft(v.I ,Nt_fft)./Nt_fft;
            v.Vs_fft = v.Vs_fft(1:floor(Nt_fft/2)+1);
            v.V_fft  = v.V_fft( 1:floor(Nt_fft/2)+1);
            v.I_fft  = v.I_fft( 1:floor(Nt_fft/2)+1);
            v.Z_fft  = v.Z_fft( 1:floor(Nt_fft/2)+1);
            
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

                Pi = (abs(v.Eth).^2 + abs(v.Eph).^2)./(120*pi)./2; %% Poynting vector

                %% radiated power in far field
                v.P_rad_ff = sum(sum(Pi.*S));
                
            end;
            
            if v.field_flag,
                
                v.Ex_field = 2.*v.Ex_field./Nt_fft;
                v.Ey_field = 2.*v.Ey_field./Nt_fft;
                v.Ez_field = 2.*v.Ez_field./Nt_fft;
                
            end;
            
            %% power entering the domain
            v.P_in = real( abs(v.V_fft).^2 ./ conj(v.Z_fft) )./2;
            
            v.locked_flag = true;
        end;
        
        %% Save
        save(fullfile(PathName,[FileNameBody '.aaf']),'-struct','v');

    catch
        
        errordlg(['Conversion operation failed at ',FileNames{i}],'Error');
        uiwait;
        
        rethrow(lasterror);
        
        break;
        
    end;

    SuccessCount = SuccessCount + 1;

end;

msgbox([num2str(SuccessCount),' file(s) successfully converted.'],'Info','modal');

if isunix,
    setappdata(0,'UseNativeSystemDialogs',DefaultDialogs);
end;

