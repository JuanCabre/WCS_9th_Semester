% street.m; reflections from ground and street buildings;
% width of street a meters;
% transmit antenna h1 m, receive ant. h2 ;
% frequency 2400 MHz;
% use ray tracing and image theory;
% Rg is reflection from the ground, assume -0.85;
% Rw is reflection from the walls, assume Rw=-0.5;
% N is number of reflections;
%
% J�rgen Bach Andersen 2002
% modified Patrick Eggers 2003-2008

clear all
close all
clf
f=2400;
lambda=300/f;
k=2*pi/lambda;
N=10;   
a=10;   % what is this ? = Distance between the walls
h1=10;  % height terminal 1
h2=2;   % height terminal 2
Rg=-0.85; % ground reflection
Rw=-0.5;  % wall reflection
d=10:1:10000; %distance along street;
E=d-d;

for i=1:N+1
R1=sqrt(d.^2+(h1-h2)^2+((i-1)*a)^2);  % distance between antenna and side images;
R2=sqrt(d.^2+(h1+h2)^2+((i-1)*a)^2);  % distance between ground image and side images;

%insert here the contribution from the various images and sum them to a total E;
Floss1 = lambda./(4.*pi.*R1); 
Floss2 = lambda./(4.*pi.*R2);
if i == 1
E = E + Floss1.*((Rw^(i-1))*exp(-j*k*R1))+ (Rg).*Floss2.*((Rw^(i-1))*exp(-j*k*R2));
else
E = E + 2*Floss1.*((Rw^(i-1))*exp(-j*k*R1))+ 2*(Rg).*Floss2.*((Rw^(i-1))*exp(-j*k*R2));
end
end

P=20*log10(abs(E));
Plos=(lambda./(4*pi*d)).^2; %direct field;
P2=10*log10(abs(Plos));

plot(log10(d),P,'r',log10(d),P2,'g','linewidth',2)
xlabel('log(d)','FontSize',15);
ylabel('Field strength dB','FontSize',15);
title('Field strength and path loss vs distance','FontSize',15);
handle=legend('Field strength', 'Path loss');
set(handle,'FontSize',12);
grid on
% axis([1 3 -70 -10]);

% discuss the result. Try also w low N=1 or 2

% construct the time response
% How will you plot the time response (power versus delay)?


% % Impulse response
% figure()
dref = 100;
ddis = sqrt(dref.^2+(h1-h2)^2);

R1=sqrt(dref.^2+(h1-h2)^2);  % distance between antenna and side images;
R2=sqrt(dref.^2+(h1+h2)^2);  % distance between ground image and side images;

Floss1 = lambda./(4.*pi.*R1); 
Floss2 = lambda./(4.*pi.*R2);

E1 = Floss1.*exp(-j*k*R1);
E2 = (Rg).*Floss2.*exp(-j*k*R2);

dis = [R1;R2];
amp = [E1;E2];

if(N>0)
for i=2:N+1
%insert here the contribution from the various images and sum them to a total E;
R1=sqrt(dref.^2+(h1-h2)^2+((i-1)*a)^2);  % distance between antenna and side images;
R2=sqrt(dref.^2+(h1+h2)^2+((i-1)*a)^2);  % distance between ground image and side images;
Floss1 = lambda./(4.*pi.*R1); 
Floss2 = lambda./(4.*pi.*R2);

E1 = 2*Floss1.*((Rw^(i-1))*exp(-j*k*R1));
E2 = 2*(Rg).*Floss2.*((Rw^(i-1))*exp(-j*k*R2));

dis = [dis [R1;R2]];
amp = [amp [E1;E2]];

end
end

delay = (dis - ddis)./(3e8);
P=20*log10(abs(amp));


delay_v = reshape(delay,[1 numel(delay)]);
P = P - P(1,1);
P_v = reshape(P,[1 numel(P)]);
P_ti = 10.^(P_v/20);

% figure;
% stem(delay_v,P_v,'r','linewidth',2)
% xlabel('Delay (s)','FontSize',15);
% ylabel('Power (dB)','FontSize',15);
% title('Impulse Response','FontSize',15);


figure;
stem(delay_v,P_ti,'b','linewidth',2)
xlabel('Delay (s)','FontSize',15);
ylabel('Power','FontSize',15);
title('Impulse Response','FontSize',15);


% impulse response for 200m

dref = 200;
ddis = sqrt(dref.^2+(h1-h2)^2);

R1=sqrt(dref.^2+(h1-h2)^2);  % distance between antenna and side images;
R2=sqrt(dref.^2+(h1+h2)^2);  % distance between ground image and side images;

Floss1 = lambda./(4.*pi.*R1); 
Floss2 = lambda./(4.*pi.*R2);

E1 = Floss1.*exp(-j*k*R1);
E2 = (Rg).*Floss2.*exp(-j*k*R2);

dis = [R1;R2];
amp = [E1;E2];

if(N>0)
for i=2:N+1
%insert here the contribution from the various images and sum them to a total E;
R1=sqrt(dref.^2+(h1-h2)^2+((i-1)*a)^2);  % distance between antenna and side images;
R2=sqrt(dref.^2+(h1+h2)^2+((i-1)*a)^2);  % distance between ground image and side images;
Floss1 = lambda./(4.*pi.*R1); 
Floss2 = lambda./(4.*pi.*R2);

E1 = 2*Floss1.*((Rw^(i-1))*exp(-j*k*R1));
E2 = 2*(Rg).*Floss2.*((Rw^(i-1))*exp(-j*k*R2));

dis = [dis [R1;R2]];
amp = [amp [E1;E2]];

end
end

delay = (dis-ddis)./(3e8);
P=20*log10(abs(amp));


delay_v = reshape(delay,[1 numel(delay)]);
P = P - P(1,1);
P_v = reshape(P,[1 numel(P)]);
P_ti = 10.^(P_v/20);

% figure;
% stem(delay_v,P_v,'r','linewidth',2)
% xlabel('Delay (s)','FontSize',15);
% ylabel('Power (dB)','FontSize',15);
% title('Impulse Response','FontSize',15);


hold all;
stem(delay_v,P_ti,'r','linewidth',2)
xlabel('Delay (s)','FontSize',15);
ylabel('Power','FontSize',15);
title('Impulse Response','FontSize',15);

