function PO_MM9
% JE Berg PO propagation - FFT
% (c) Patrick Eggers 2008

close all
clear all

%%%%% Scenario setting %%%%%%%%%%%%%%%%%
f=900*10^6;
c=3*10^8;
lambda=c/f;
k=2*pi/lambda;
delta_y=lambda/3 %screen aperture resolution
delta_x=10; % screen spacing
eps=0.01; % acceptable field error pertubation.. look out , too small blows up array sizes
%% define your height of screen
h_screen=5;
%% choose angle of incidence, start w 0
alpha=(0/360)*2*pi; % in radians
r_tx=1000; % transmitter/source speration to 1st screen edge
E0=10; %start field strength

N=2; %last screen

%%%%%%%%%% procedure %%%%%%%%%%%%%%

%1 determine array size  ny . remember Fresnel zone, terminator/neutraliser ..
%what else wrt FFT to think of??
%walfisch & bertoni truncation paramters
yc=N*delta_x*tan(alpha) + sqrt(N*lambda*delta_x/(2*pi*eps))*sec(alpha)
w=sqrt(lambda*delta_x)
ny=round(100) %modifiy .. how many lambda or meters??/delta_y /// yc + 3W + how deep times 2
ny= round(2*yc/delta_y + 6*w/delta_y + h_screen/delta_y);

%1a set ky , i..e space frequency variable, what range does it have and
%why?, what units?
nyy=52; %try to vary this .. what happens why? .. what should it be?
ky=(((1:ny)-nyy)/ny).*2*pi/delta_y;

red_zone= yc+3*w

%2 initiate arrays
E(1,:)=zeros(size(ky));
%% now fill in the start field over 1st screen following JE Berg
% remember field array above screen is in a plane.. map the cyl. field to
% this plane. Also remember = 0 at/below screen

% Height of the Source
H= (h_screen + tan(alpha)*r_tx) + red_zone;
% Distance from source
r(1:ny)= sqrt( (H -[1:ny].*delta_y).^2 + (r_tx)^2);
% Initial Field
 E(1,:)=1./(r.^0.5) .* exp(-i*k.*r);
%E(1,:)=1;

plot(abs(E(1,:)),'-b');
hold on; grid on

%3 condition array before FFT.. what do you do here?
% We truncate and slowly kill (Neutralizer)
E(1,1:round(h_screen/delta_y + red_zone/delta_y))=0;
y= 0:delta_y:red_zone;
length(y)
length(E(1,end-length(y)+1:end))

E(1,end+1-length(y):end)=E(1,end+1-length(y):end).*neutraliser(y,yc,w);
plot(abs(E(1,:)),'--r')



%4 FFT of the E field 1st screen 
E(1,:)=fft(E(1,:),ny);


%5 propagator
% Propagator difined by Patrick
w_prop_old=exp(-((ky.^2-k.^2).^0.5) * delta_x); % why not as in paper : ky^2<=k^2 vs ky^2>k^2 ???

% Propagator defined by the paper
for l=1:length(ky)
    if ky(l)^2<=k^2
        w_prop(l)=exp(-1i*((-ky(l)^2+k.^2).^0.5) * delta_x);
    else
        w_prop(l)=exp(-1*((-k.^2+ky(l)^2).^0.5) * delta_x);
    end
end
% Both propagators provide the same results. it means that we need to
% correct the fliping introduced by Matlab

%6 multiply and inverse FFT
% Due to the shift of the FFT implemented by Matlab, we need to conjugate
% the propagator
w_prop=circshift(w_prop,[1 -52])
E(2,:)=ifft((w_prop).*((E(1,:))),ny);


%7 extract and plot results at position of last screen
%% identify index corresponiong to you angle of incidence/grazing incidence
% plot(abs(E(2,:)),'-m')

%% map to total field dependency for 'true life' launched spherical  wave

plot(abs(E(2,:)),'-k')
legend('Transmited field','Truncated and neutralized field','Resulting field after the screen')
title('Fields When rtx=100')

1/delta_y
% plot([red_zone/delta_y,red_zone/delta_y],[0,0.4],'k')

%what do observe wrt expected diffraction field?? - what does it indicate
%we might have to do?

hold off

% figure
% % plot(phase(flip(conj(w_prop))))
% plot(phase(w_prop))
% 
% hold on
% plot(phase(E(1,:)),'r')
% 
% figure
% % plot(abs(flip(conj(w_prop))))
% plot(abs(w_prop))
% hold on
% plot(abs(E(1,:)),'r')


%8 when this works (try differnt alpha and r_tx) .. try to to one more
%screen (and try to vary screen heights..1st all same , then different)

 
function [Neutral] = neutraliser(y,yc,w)
%from walfisch & bertoni
%y is array storing y coordinates
%yc is where soft transition starts to kick in
%factor 1/2 due to JE Berg
   Neutral=exp(-(y-yc).^2/w^2); % for yc<y<yc+3*w
   Neutral=(y<=(yc+3*w)).*Neutral; % 0 above neutralising zone
  Neutral=(y<yc)+(y>=yc).*Neutral; %(1->)1/2 below neutralising zone
  disp(length(Neutral))
 


