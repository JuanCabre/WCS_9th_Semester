% templat exercise 1A
%=============================================================================
% Copyright (c) by Patrick Eggers 2008. All rigts reserved.
% This copyright notice may NOT be removed from the program
% and must accompany any part of the program being copied .
%==============================================================================
clear;
clf;

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
a1=sin(rad).^30;

% 2nd part of pattern
a2=0.25*cos(rad-pi/4).^80;

% combine thge two in to one common pattern 'a' .. how ??

% dB and plot
a_dB=10.*log10(a + 0.001);

%plot(a1_dB);
% or polar(mmmmmm,nnnnnn) what should arguments be?
%axis([0 2*pi -60 0])

%-------------- question 1  ---------------------
%calculate C0/I0

%CS/IS

%2)------------- question 2 ------------------------
% make angular correlation/convolution of patterna and new environment
% how ??? .. think of convolution --- vs frequency domain
% NOTE dB is only for plotting/visualisation .. all manipulations are in liear power
% (OR du you have an altrenativeto solving the question .. then do both)

%initiate Cs and Is

% perform correlation of a and Cs

% perform correlation of a and Is

% plot result of correlations dBs _> polar plot, remeber offset and negative
% value trunking .. why??

% read out the values asked for in question .. how /where do you read and
% compare??