% template exercise 1A
%=============================================================================
% Copyright (c) by Patrick Eggers 2008. All rigts reserved.
% This copyright notice may NOT be removed from the program
% and must accompany any part of the program being copied .
%==============================================================================
clear all;
clf;
close all

%% Transmitted Power in dB

Tx_power = 30; 


j=sqrt(-1);

test_power=0;

del_ang=1;%2.5; % angular resolution
%check what happens at -180
lo_deg=1;
hi_deg=360;
phi=lo_deg:del_ang:hi_deg; % full circle
deg2rad=pi/180;
rad2deg=180/pi;
rad=phi*deg2rad; % in radians
phi_0=0;

% 1st part of pattern
a=zeros(size(rad)); % why do we do this, what happen w this statement?
a1=sin(rad).^30;  % 53 to 128 deg

% 2nd part of pattern
a2=0.25*cos(rad-pi/4).^80; % 24 to 66 deg

% combine thge two in to one common pattern 'a' .. how ??

a(53:128)=a1(53:128);
a(24:66)=a(24:66)+a2(24:66);

%% dB and plot (polar and cartesian plot of the antenna pattern)
a_dB=10.*log10(a + 0.001);
a_dB_polar=[ transpose(rad) transpose(a_dB+Tx_power*ones(size(a_dB,1),size(a_dB,2)))];

figure
plot(rad,a_dB);
title('Antenna Pattern Cartesian Plot')
xlabel('Angle [rad]')
axis([0, pi, -30, 5])
ylabel('Power [dB]')
grid on

figure
polar(a_dB_polar(:,1),a_dB_polar(:,2))  
title('Antenna Pattern Polar Plot')


%% Correlation with dirac and distribution by using the Fourier Transform

dirac90=zeros(size(rad));
dirac60=zeros(size(rad));
dirac90(90)=1;
dirac60(60)=3/4;

rect1=zeros(size(rad));
rect1(80:100)=1;

rect2=zeros(size(rad));
rect2(54:67)=3/4;

Rect1=fft(rect1);
Rect2=fft(rect2);

Dirac= fft(dirac90);
Dirac2= fft(dirac60);


A=fft(a);
RESULT=ifft((A).*Dirac);
RESULT2=ifft((A).*Dirac2);
RESULT_dB=10.*log10(RESULT + 0.001);
RESULT2_dB=10.*log10(RESULT2 + 0.001);


RESULT3=ifft((A).*Rect1);
RESULT4=ifft((A).*Rect2);
RESULT3_dB=10.*log10(RESULT3 + 0.001);
RESULT4_dB=10.*log10(RESULT4 + 0.001);

RESULT_dB_polar=[ transpose(rad) transpose(RESULT_dB+Tx_power*ones(size(RESULT_dB,1),size(RESULT_dB,2)))];
RESULT2_dB_polar=[ transpose(rad) transpose(RESULT2_dB+Tx_power*ones(size(RESULT2_dB,1),size(RESULT2_dB,2)))];

figure
polar(RESULT_dB_polar(:,1),RESULT_dB_polar(:,2))
hold on
polar(RESULT2_dB_polar(:,1),RESULT2_dB_polar(:,2),'r')
hold off


RESULT3_dB_polar=[ transpose(rad) transpose(RESULT3_dB+Tx_power*ones(size(RESULT_dB,1),size(RESULT_dB,2)))];
RESULT4_dB_polar=[ transpose(rad) transpose(RESULT4_dB+Tx_power*ones(size(RESULT2_dB,1),size(RESULT2_dB,2)))];
figure

polar(RESULT3_dB_polar(:,1),RESULT3_dB_polar(:,2))
hold on
polar(RESULT4_dB_polar(:,1),RESULT4_dB_polar(:,2),'r')
hold off


%% SDMA Exercise: Same Antenna Pattern as before, now repeated multiple times with a certain offset angle

    a_SDMA=zeros(size(rad));
    a_SDMA(53:128)=transpose(a1(53:128));
    a_SDMA(24:66)=a_SDMA(24:66)+(a2(24:66));

    a_SDMAdB=transpose(10.*log10(a_SDMA + 0.001));
    
    
%% Plot of the antenna pattern, to see for which offset we reach which interference level

figure
polar(a_dB_polar(:,1),a_dB_polar(:,2))  
hold on
h9=plot([0 cosd(24.5+90)*40],[0 sind(24.5+90)*40],'color',[0.85 0.33 0.1]);
h10=plot([0 cosd(-24.5+90)*40],[0 sind(-24.5+90)*40],'color',[0.85 0.33 0.1]);
h1=plot([0 cosd(29+90)*40],[0 sind(29+90)*40],'color',[0 0.45 0.74]);
h2=plot([0 cosd(-29+90)*40],[0 sind(-29+90)*40],'color',[0 0.45 0.74]);
h3=plot([0 cosd(52.8+90)*40],[0 sind(52.8+90)*40],'color',[0.49 0.18 0.56]);
h4=plot([0 cosd(-52.8+90)*40],[0 sind(-52.8+90)*40],'color',[0.49 0.18 0.56]);
h5=plot([0 cosd(56+90)*40],[0 sind(56+90)*40],'color',[0.47 0.67 0.19]);
h6=plot([0 cosd(-56+90)*40],[0 sind(-56+90)*40],'color',[0.47 0.67 0.19]);
h7=plot([0 cosd(58+90)*40],[0 sind(58+90)*40],'color',[0.64 0.08 0.18]);
h8=plot([0 cosd(-58+90)*40],[0 sind(-58+90)*40],'color',[0.64 0.08 0.18]);
h11=plot([0 cosd(60+90)*40],[0 sind(60+90)*40],'color',[0.3 0.75 0.93]);
h12=plot([0 cosd(-60+90)*40],[0 sind(-60+90)*40],'color',[0.3 0.75 0.93]);
hold off

legend([h9,h1,h3,h5,h7,h9,h11],'24.5°',...
    '29°',...
    '52.8°',...
    '56°',...
    '58°',...
    '60°')
    

title('Powers at different Offset-Angles')    
   

%% SDMA plots with different offsets
    %% offset 29°  -  target CIR: 9dB
    
    offset=29;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    

    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  

    title('Offset: 29°, CIR>9dB')
    
 %% offset 24.5°  -  target CIR: 9dB
    
    offset=24.5;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    

    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  

    title('Offset: 24.5°, CIR<9dB')
    
    
%% Offset 56°   -  target CIR: 9dB

    offset=56;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    

    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
 
    title('Offset: 56°, CIR>9dB')
    
    
%% Offset 52.76°   -  target CIR: 9dB

    offset=52.76;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    
    
    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    title('Offset: 52.8°, CIR>9dB')
    
    
    %% Offset 60.09°   -  target CIR: 15dB

    offset=60;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    
    
    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2)) 
    
    
    title('Offset: 60°, CIR>15dB')
    
    
    %% Offset 58.09°   -  target CIR: 15dB

    offset=58;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    
    
    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi) (a_SDMAdB+Tx_power*ones(size(a_SDMAdB,1),size(a_SDMAdB,2)))];
    polar(a_SDMAdB_polar(:,1),a_SDMAdB_polar(:,2)) 
    
    
    title('Offset: 58°, CIR>15dB')
    

    %%
    
    %%
    
    %%
    
test=transpose(RESULT3_dB);
    
    
RESULT3_dB_polar=[ transpose(rad-pi/2) transpose(RESULT3_dB+Tx_power*ones(size(RESULT3_dB,1),size(RESULT3_dB,2)))];

figure
polar(RESULT3_dB_polar(:,1),RESULT3_dB_polar(:,2))

title('Equivalent radiation pattern')
    

figure
plot(rad-pi/2,RESULT3_dB-max(RESULT3_dB));
title('Antenna Pattern Cartesian Plot')
xlabel('Angle [rad]')
axis([0, pi, -15, 0])
ylabel('Power [dB]')
grid on

%%

    offset=59;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    
    
    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2)) 
    
    
    title('Offset: 59°, CIR>9dB')

     
    
    %%
    
    %%

    offset=63.3;
    o1=offset;
    o2=0;
    o3=-offset;
    o4=2*offset;
    o5=-2*offset;
    
    figure
    
    
    
    a_SDMAdB_polar=[ transpose(rad+o1/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    hold on
    
    a_SDMAdB_polar=[ transpose(rad+o2/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o3/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o4/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2))  
    
    
    a_SDMAdB_polar=[ transpose(rad+o5/180*pi-pi/2) (RESULT3_dB_polar+Tx_power*ones(size(RESULT3_dB_polar,1),size(RESULT3_dB_polar,2)))];
    polar(a_SDMAdB_polar(:,1),RESULT3_dB_polar(:,2)) 
    
    
    title('Offset: 63.3°, CIR>15dB')