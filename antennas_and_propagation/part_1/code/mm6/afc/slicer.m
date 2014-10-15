function slicer(Data,DataDescription,FigureTitle)
%SLICER Displays 3-D arrays in slices
%
% SLICER(DATA) where DATA is a 3-D array
%
% SLICER(DATA,DESC) adds DESC string as a description of colorbar
%
% SLICER(DATA,DESC,TITLE) uses TITLE string as a figure title

% Copyright 2007-- Ondrej Franek
% $Revision v1.06$ $(14-05-2008)$

AvailableColormaps = ...
    {'autumn','bone','colorcube','cool','copper','flag','gray','hot', ...
     'hsv','jet','lines','pink','prism','spring','summer','white','winter'};

if any(imag(Data(:))), Data = abs(Data); end;
[SizeX,SizeY,SizeZ] = size(Data);

FiniteMask = isfinite( Data(:) ); % selecting only the finite elements

% in case the array does not have finite elements
if ~any(FiniteMask),
    MinData = -1; MaxData = 1;
else
    MinData = double( min(Data(FiniteMask)) ); % limits of the data
    MaxData = double( max(Data(FiniteMask)) );
end;

% in case the array is homogeneous
if MinData==MaxData, MaxData = MaxData + 1; MinData = MinData - 1; end;

DispRange = [MinData MaxData]; % display range - initialized to the limits
Eps = ( MaxData - MinData )/10000; % minimum display range (zero doesn't work)

Slice = [1 1 1]; % intersection of the slices
Point = [1 1 1]; % current point

% flags declaration - shared among nested functions
MovingIntersection  = false;
MovingDisplayWindow = false;
MovingLowerLimit    = false;
MovingUpperLimit    = false;
Zooming             = false;
Panning             = false;
InAxes              = false;

% shared variables declaration
hCurrentAxes         = [];
hColorbar            = [];
OldPointerLocation   = [];
AxesOriginalPosition = [];

% creating the figure
if nargin < 3, FigureTitle = 'Slicer'; end; % default figure title
hFigure = figure( ...
    'Units','normalized', ...
    'Position',[.2 .2 .6 .6], ...
    'WindowButtonMotionFcn',@PointerMotion, ...
    'WindowButtonUpFcn',@ButtonUp, ...
    'NumberTitle','off', ...
    'Name',FigureTitle, ...
    'Toolbar','figure', ...
    'Colormap',jet(256)); % default colormap
Colormap = get(hFigure,'Colormap');

ZoomMode = zoom(hFigure); % camera modes
PanMode  = pan(hFigure);

% slice in ZX plane
hAxesZX  = axes('Units','normalized','OuterPosition',[0.1 .51 .43 .47]);
hImageZX = image('CData',permute(Data(:,1,:),[3 1 2]), ...
                 'CDataMapping','scaled', ...
                 'XData',[1 SizeX]-0.5, ...
                 'YData',[1 SizeZ]-0.5);
xlabel x;
ylabel z;
title('ZX plane');

% slice in YZ plane
hAxesYZ  = axes('Units','normalized','OuterPosition',[.55 .5 .43 .47]);
hImageYZ = image('CData',permute(Data(1,:,:),[3 2 1]), ...
                 'CDataMapping','scaled', ...
                 'XData',[1 SizeY]-0.5, ...
                 'YData',[1 SizeZ]-0.5);
xlabel y;
ylabel z;
title('YZ plane');

% slice in XY plane
hAxesXY  = axes('Units','normalized','OuterPosition',[0.1 0.02 .43 .47]);
hImageXY = image('CData',permute(Data(:,:,1),[2 1 3]), ...
                 'CDataMapping','scaled', ...
                 'XData',[1 SizeX]-0.5, ...
                 'YData',[1 SizeY]-0.5);
xlabel x;
ylabel y;
title('XY plane');

% properties common for all axes
hAllAxes = [hAxesZX hAxesYZ hAxesXY];
set(hAllAxes,'Box','on', ...
             'TickDir','out', ...
             'CLim',DispRange);
axis(hAllAxes,'image');
set([hImageZX hImageYZ hImageXY],'ButtonDownFcn',@ButtonDownOnImages);

% draw slice lines for the first time
hLineZX_Z = line([0 SizeX],Slice([3 3])-.5,'Parent',hAxesZX);
hLineZX_X = line(Slice([1 1])-.5,[0 SizeZ],'Parent',hAxesZX);
hLineYZ_Z = line([0 SizeY],Slice([3 3])-.5,'Parent',hAxesYZ);
hLineYZ_Y = line(Slice([2 2])-.5,[0 SizeZ],'Parent',hAxesYZ);
hLineXY_Y = line([0 SizeX],Slice([2 2])-.5,'Parent',hAxesXY);
hLineXY_X = line(Slice([1 1])-.5,[0 SizeY],'Parent',hAxesXY);
set([hLineZX_Z hLineZX_X hLineYZ_Z hLineYZ_Y hLineXY_Y hLineXY_X], ...
    'LineWidth',2, ...
    'LineStyle',':', ...
    'Color','w', ...
    'ButtonDownFcn',@ButtonDownOnImages);

% colorbar
hColorbar = axes( ...
    'Units','normalized', ...
    'OuterPosition',[.02 .02 .06 .96], ...
    'Box','on', ...
    'XTick',[], ...
    'TickDir','out', ...
    'CLim',[0 1], ...
    'XLim',[.4 1.5], ... % 0.4 so that the left borderline is visible
    'YLim',DispRange, ...
    'ButtonDownFcn',@ButtonDownOnColorbar);
hColormapImage = image( ...
    'CData',(0:0.001:1)', ... % 1000 points for smoothness
    'CDataMapping','scaled', ...
    'YData',DispRange, ...
    'ButtonDownFcn',@ButtonDownOnColorbar);
if nargin > 1, % colorbar title
    hColorbarTitle = title(hColorbar,DataDescription);
    %% title must be elevated to avoid collision with y-axis multiplier
    set(hColorbarTitle,'Units','pixels');
    set(hColorbarTitle,'Position',get(hColorbarTitle,'Position') + [0 20 0]);
    set(hColorbarTitle,'Units','normalized'); % ...because of resizing
end;

set(hColorbar,'Units','pixels');
ColorbarPosition = get(hColorbar,'Position');
set(hColorbar,'Units','normalized'); % ...because of resizing

% colorbar adjustable limits
hLowerLimit = uicontrol( ...
    'Style','frame', ...
    'Units','pixels', ...
    'Position',[ColorbarPosition(1:2)-[1 3] 18 5], ...
    'Enable','inactive', ...
    'ButtonDownFcn',@ButtonDownOnLimits);
hUpperLimit = uicontrol( ...
    'Style','frame', ...
    'Units','pixels', ...
    'Position',[ColorbarPosition(1:2)-[1 3]+[0 ColorbarPosition(4)] 18 5], ...
    'Enable','inactive', ...
    'ButtonDownFcn',@ButtonDownOnLimits);

% display panel
set(hFigure,'Units','pixels');
FigurePosition = get(hFigure,'Position');
PanelSize = [230 150];
% panel will be initially positioned in the lower right corner
PanelCurrentPosition = [FigurePosition(3)-PanelSize(1)-10 10 PanelSize];
PanelOriginalPosition = PanelCurrentPosition;
hPanel = uipanel(hFigure, ...
    'Units','pixels', ...
    'Position',PanelCurrentPosition, ...
    'BorderType','beveledout');

hTextDescribe(1) = uitext(hPanel,'Section Coordinates:',[10 100 120 17]);
hTextDescribe(2) = uitext(hPanel,'Cursor Coordinates:' ,[10  75 120 17]);
hTextDescribe(3) = uitext(hPanel,'Cursor Voxel Value:' ,[10  50 120 17]);
hTextDescribe(4) = uitext(hPanel,'Colormap:'           ,[10  10 120 17]);

hTextCoords(1) = uitext(hPanel,'x',[120 125 30 17]);
hTextCoords(2) = uitext(hPanel,'y',[155 125 30 17]);
hTextCoords(3) = uitext(hPanel,'z',[190 125 30 17]);

hDispSection(1) = uitext(hPanel,num2str(Slice(1)),[120 100 30 17]);
hDispSection(2) = uitext(hPanel,num2str(Slice(2)),[155 100 30 17]);
hDispSection(3) = uitext(hPanel,num2str(Slice(3)),[190 100 30 17]);

hDispCursor(1) = uitext(hPanel,'',[120 75 30 17]);
hDispCursor(2) = uitext(hPanel,'',[155 75 30 17]);
hDispCursor(3) = uitext(hPanel,'',[190 75 30 17]);

hDispVoxel = uitext(hPanel,'',[120 50 65 17]); % voxel coordinates
hDispColor = uicontrol(hPanel, ... % voxel color
    'Style','frame', ...
    'Units','pixels', ...
    'Position',[197 50 17 17], ...
    'BackgroundColor','default');

hDispDivisionLine = uipanel( ...
    'Parent',hPanel, ...
    'Units','pixels', ...
    'Position',[10 40 PanelSize(1)-20 2]);

hDispColormap = uicontrol(hPanel, ...
    'Style','popupmenu', ...
    'Units','pixels', ...
    'Position',[70 15 80 17], ...
    'String',AvailableColormaps, ...
    'BackgroundColor','w', ...
    'Value',10, ...
    'Callback',@ChangeColormap);
    function ChangeColormap(src,eventdata)
        val = get(hDispColormap,'Value');
        Colormap = eval([AvailableColormaps{val},'(256)']);
        set(hFigure,'Colormap',Colormap);
    end

set([hDispSection hDispCursor hDispVoxel], ...
    'FontWeight','bold');
set([hTextCoords hDispSection hDispCursor hDispVoxel], ...
    'HorizontalAlignment','center');
set([hTextDescribe hTextCoords hDispSection hDispCursor hDispVoxel hDispColor], ...
    'Enable','inactive');
set([hPanel hTextDescribe hTextCoords hDispSection hDispCursor hDispVoxel hDispColor], ...
    'ButtonDownFcn',@ButtonDownOnDisplay);
set(hPanel,'Units','normalized'); % to enable automatic placement on resize
set(hFigure,'ResizeFcn',@FigureResize);



%%%%%%%%%%%%%%%%%%%%%%%%
%% Callback Functions %%
%%%%%%%%%%%%%%%%%%%%%%%%

    function FigureResize(src,eventdata)
        
        % keep colorbar width at 15 pixels
        set(hColorbar,'Units','pixels');
        ColorbarPosition = get(hColorbar,'Position');
        ColorbarPosition(3) = 15;
        set(hColorbar,'Position',ColorbarPosition);
        set(hColorbar,'Units','normalized');
        
        % place limits to correct positions
        set([hLowerLimit hUpperLimit], ...
            'Units','pixels', ...
            'Position',[ColorbarPosition(1)-1 0 18 5]);
        RedrawLimit(hLowerLimit,DispRange(1));
        RedrawLimit(hUpperLimit,DispRange(2));
        
        % keep display panel size, but centered at new position
        set(hPanel,'Units','pixels');
        PanelCurrentPosition = get(hPanel,'Position');
        PanelCurrentPosition(1:2) = PanelCurrentPosition(1:2) + ...
            (PanelCurrentPosition(3:4) - PanelSize)./2;
        PanelCurrentPosition(3:4) = PanelSize;
        set(hPanel,'Position',PanelCurrentPosition);
        set(hPanel,'Units','normalized');
        
    end


    function ButtonDownOnDisplay(src,eventdata)

        % dragging display window with right mouse button
        if ~strcmp(get(hFigure,'SelectionType'),'alt'), return; end;
        setptr(hFigure,'closedhand'); % undocumented pointer
        set(hFigure,'Units','pixels');
        PanelOriginalPosition = PanelCurrentPosition - ...
                                [get(hFigure,'CurrentPoint') 0 0];
        MovingDisplayWindow = true;

    end


    function ButtonDownOnLimits(src,eventdata)

        switch src
            case hLowerLimit,
                set(hFigure,'Pointer','bottom');
                MovingLowerLimit = true;
            case hUpperLimit,
                set(hFigure,'Pointer','top');
                MovingUpperLimit = true;
        end;

    end


    function ButtonDownOnImages(src,eventdata)
        
        switch get(hFigure,'SelectionType')
            
            case 'normal' % normal click = moving intersection of slices
                RedrawImages;
                MovingIntersection = true;

            case 'extend' % center or left and right mouse button = zooming
                setptr(hFigure,'glass'); % undocumented pointer
                OldPointerLocation = get(0,'PointerLocation');
                Zooming = true;
                switch src % which axes will be zoomed ?
                    case {hImageZX,hLineZX_Z,hLineZX_X}, hCurrentAxes = hAxesZX;
                    case {hImageYZ,hLineYZ_Y,hLineYZ_Z}, hCurrentAxes = hAxesYZ;
                    case {hImageXY,hLineXY_X,hLineXY_Y}, hCurrentAxes = hAxesXY;
                end;

            case 'alt' % right mouse button = panning
                setptr(hFigure,'closedhand'); % undocumented pointer
                Panning = true;
                switch src % which axes will be panned ?
                    case {hImageZX,hLineZX_Z,hLineZX_X}, hCurrentAxes = hAxesZX;
                    case {hImageYZ,hLineYZ_Y,hLineYZ_Z}, hCurrentAxes = hAxesYZ;
                    case {hImageXY,hLineXY_X,hLineXY_Y}, hCurrentAxes = hAxesXY;
                end;
                set([hFigure,hCurrentAxes],'Units','pixels');
                AxesOriginalPosition = ...
                    get(hCurrentAxes,'Position') - [get(hFigure,'CurrentPoint') 0 0];

            case 'open' % doubleclick = return to original zoom
                switch src
                    case {hImageZX,hLineZX_Z,hLineZX_X},
                        set(hAxesZX,'CameraViewAngleMode','auto');
                    case {hImageYZ,hLineYZ_Y,hLineYZ_Z},
                        set(hAxesYZ,'CameraViewAngleMode','auto');
                    case {hImageXY,hLineXY_X,hLineXY_Y},
                        set(hAxesXY,'CameraViewAngleMode','auto');
                end;
            
        end;

    end


    function ButtonDownOnColorbar(src,eventdata)
        
        if strcmp(get(hFigure,'SelectionType'),'alt'), % right mouse button = panning
                
            setptr(hFigure,'closedhand'); % undocumented pointer
            Panning = true;
            hCurrentAxes = hColorbar; %% colorbar is the panned axes
            set([hFigure,hCurrentAxes],'Units','pixels');
            AxesOriginalPosition = ...
                get(hCurrentAxes,'Position') - [get(hFigure,'CurrentPoint') 0 0];

        end;

    end


    function PointerMotion(src,eventdata)
        
        % moving mouse pointer inside the figure can mean many things,
        % according to the flag which was set when pressing mouse button:
        
        if MovingDisplayWindow,
            set(hFigure,'Units','pixels');
            PanelCurrentPosition = PanelOriginalPosition + ...
                                   [get(hFigure,'CurrentPoint') 0 0];
            set(hPanel,'Units','pixels', ...
                       'Position',PanelCurrentPosition);
            set(hPanel,'Units','normalized');
            return;
            
        elseif MovingLowerLimit,
            CBarPoint = get(hColorbar,'CurrentPoint');
            if     CBarPoint(1,2) <  MinData,      DispRange(1) = MinData;
            elseif CBarPoint(1,2) >= DispRange(2), DispRange(1) = DispRange(2) - Eps;
            else                                   DispRange(1) = CBarPoint(1,2);
            end;
                set(hColormapImage,'YData',DispRange);
                RedrawLimit(hLowerLimit,DispRange(1));
                set(hAllAxes,'CLim',DispRange);
            return;
            
        elseif MovingUpperLimit,
            CBarPoint = get(hColorbar,'CurrentPoint');
            if     CBarPoint(1,2) >  MaxData,      DispRange(2) = MaxData;
            elseif CBarPoint(1,2) <= DispRange(1), DispRange(2) = DispRange(1) + Eps;
            else                                   DispRange(2) = CBarPoint(1,2);
            end;
                set(hColormapImage,'YData',DispRange);
                RedrawLimit(hUpperLimit,DispRange(2));
                set(hAllAxes,'CLim',DispRange);
            return;

        elseif Zooming,
            NewPointerLocation = get(0,'PointerLocation');
            % double the zoom on every 30 pixels drag
            ZoomFactor = 2^(( NewPointerLocation(2) - OldPointerLocation(2) )/30);
            OldPointerLocation = NewPointerLocation;
            camzoom(hCurrentAxes,ZoomFactor);
            return;
            
        elseif Panning,
            AxesCurrentPosition = AxesOriginalPosition + ...
                                  [get(hFigure,'CurrentPoint') 0 0];
            set(hCurrentAxes,'Units','pixels', ...
                             'Position',AxesCurrentPosition);
            if hCurrentAxes==hColorbar, %% panning the colorbar
                RedrawLimit(hLowerLimit,DispRange(1));
                RedrawLimit(hUpperLimit,DispRange(2));
            end;
            return;

        end;
        
        % to prevent interference with the built-in zoom, pan and camera
        if strcmp(get(ZoomMode,'Enable'),'on') || ...
           strcmp(get(PanMode ,'Enable'),'on'), return; end;
        camera_mode = cameratoolbar('GetMode');
        if ~isempty(camera_mode), return; end;
        
        % mouse over colorbar limits lights up the appropriate pointer icon
        set([hFigure,hLowerLimit,hUpperLimit],'Units','pixels');
        FigurePoint = get(hFigure,'CurrentPoint');
        LowerLimitPosition = get(hLowerLimit,'Position');
        UpperLimitPosition = get(hUpperLimit,'Position');
        if     ( FigurePoint >= LowerLimitPosition(1:2) - [5 5] ) & ...
               ( FigurePoint <= LowerLimitPosition(1:2) + [5 5] + [18 5] ),
            set(hFigure,'Pointer','bottom');
            return;
        elseif ( FigurePoint >= UpperLimitPosition(1:2) - [5 5] ) & ...
               ( FigurePoint <= UpperLimitPosition(1:2) + [5 5] + [18 5] ),
            set(hFigure,'Pointer','top');
            return;
        end;
                          

        InAxes = false;

        % mouse pointer is in axes zx
        AxesPoint = get(hAxesZX,'CurrentPoint');
        AxesPoint = AxesPoint(1,1:2);
        if ( AxesPoint(1) >= 0 ) && ( AxesPoint(1) < SizeX ) && ...
           ( AxesPoint(2) >= 0 ) && ( AxesPoint(2) < SizeZ )
            AxesPoint = round(AxesPoint+0.5);
            Point = [AxesPoint(1) Slice(2) AxesPoint(2)];
            InAxes = true;
        end;

        % mouse pointer is in axes yz
        AxesPoint = get(hAxesYZ,'CurrentPoint');
        AxesPoint = AxesPoint(1,1:2);
        if ( AxesPoint(1) >= 0 ) && ( AxesPoint(1) < SizeY ) && ...
           ( AxesPoint(2) >= 0 ) && ( AxesPoint(2) < SizeZ )
            AxesPoint = round(AxesPoint+0.5);
            Point = [Slice(1) AxesPoint(1) AxesPoint(2)];
            InAxes = true;
        end;

        % mouse pointer is in axes xy
        AxesPoint = get(hAxesXY,'CurrentPoint');
        AxesPoint = AxesPoint(1,1:2);
        if ( AxesPoint(1) >= 0 ) && ( AxesPoint(1) < SizeX ) && ...
           ( AxesPoint(2) >= 0 ) && ( AxesPoint(2) < SizeY )
            AxesPoint = round(AxesPoint+0.5);
            Point = [AxesPoint(1) AxesPoint(2) Slice(3)];
            InAxes = true;
        end;

        if InAxes,
            
            % update display of the cursor position
            value = double(Data(Point(1),Point(2),Point(3)));
            set(hFigure,'Pointer','crosshair');
            set(hDispCursor,{'String'},cellstr(num2str(Point(:))));
            set(hDispVoxel,'String',num2str(value));
            
            % scaled index into colormap
            CMapIndex = fix((value-DispRange(1))/(DispRange(2)-DispRange(1))*256)+1;
            
            % cropping if out of range
            if isnan(CMapIndex), CMapIndex = 1;
            elseif CMapIndex < 1  , CMapIndex = 1;
            elseif CMapIndex > 256, CMapIndex = 256;
            end;
            
            % setting the colored patch in display panel
            set(hDispColor,'BackgroundColor',Colormap(CMapIndex,:));
            if MovingIntersection, RedrawImages; end;
            return;
            
        else
            
            set(hFigure,'Pointer','default');
            set([hDispCursor hDispVoxel],'String','');
            set(hDispColor,'BackgroundColor','default');
            
        end;

    end


    function ButtonUp(src,eventdata)

        if     MovingIntersection,  MovingIntersection  = false; return;
        elseif MovingDisplayWindow, MovingDisplayWindow = false;
        elseif MovingLowerLimit,    MovingLowerLimit    = false;
        elseif MovingUpperLimit,    MovingUpperLimit    = false;
        elseif Zooming,             Zooming             = false;
        elseif Panning,             Panning             = false;
            set(hCurrentAxes,'Units','normalized');
            if hCurrentAxes==hColorbar, %% panning the colorbar
                RedrawLimit(hLowerLimit,DispRange(1));
                RedrawLimit(hUpperLimit,DispRange(2));
            end;
        end;
        set(hFigure,'Pointer','default');

    end


    function RedrawImages

        Slice = Point;
        % update display of the section coordinates
        set(hDispSection,{'String'},cellstr(num2str(Slice(:))));
        set(hImageZX,'CData',permute(Data(:,Slice(2),:),[3 1 2]));
        set(hImageYZ,'CData',permute(Data(Slice(1),:,:),[3 2 1]));
        set(hImageXY,'CData',permute(Data(:,:,Slice(3)),[2 1 3]));
        set(hLineZX_Z,'XData',[0 SizeX],'YData',Slice([3 3])-.5);
        set(hLineZX_X,'XData',Slice([1 1])-.5,'YData',[0 SizeZ]);
        set(hLineYZ_Z,'XData',[0 SizeY],'YData',Slice([3 3])-.5);
        set(hLineYZ_Y,'XData',Slice([2 2])-.5,'YData',[0 SizeZ]);
        set(hLineXY_Y,'XData',[0 SizeX],'YData',Slice([2 2])-.5);
        set(hLineXY_X,'XData',Slice([1 1])-.5,'YData',[0 SizeY]);

    end


    function RedrawLimit(handle,Level)

        set(hColorbar,'Units','pixels');
        ColorbarPosition = get(hColorbar,'Position');
        set(hColorbar,'Units','normalized');

        set(handle,'Units','pixels');
        LimitPosition = get(handle,'Position');
        LimitPosition(1) = ColorbarPosition(1) - 1;
        LimitPosition(2) = ColorbarPosition(2) + ...
            ( Level - MinData )./( MaxData - MinData ) ...
            .*ColorbarPosition(4) - 3;
        set(handle,'Position',LimitPosition);
        
    end


end % end of main function



%% uitext = uicontrol('Style','text')
function handle = uitext(parent,string,position,varargin)
    handle = uicontrol(parent, ...
        'Style','text', ...
        'String',string, ...
        'HorizontalAlignment','left', ...
        'Units','pixels', ...
        'Position',position, ...
        varargin{:});
end
