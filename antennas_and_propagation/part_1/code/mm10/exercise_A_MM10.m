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
d=10:1:1000; %distance along street;
E=d-d;

for i=1:N+1
R1=sqrt(d.^2+(h1-h2)^2+((i-1)*a)^2);  % distance between antenna and side images;
R2=sqrt(d.^2+(h1+h2)^2+((i-1)*a)^2);  % distance between ground image and side images;

%insert here the contribution from the various images and sum them to a total E;
if i == 1
E = E + (Rw^(i-1))*exp(-j*k*R1)+ (Rg)*(Rw^(i-1))*exp(-j*k*R2);
else
E = E + 2*(Rw^(i-1))*exp(-j*k*R1)+ 2*(Rg)*(Rw^(i-1))*exp(-j*k*R2);
end
end

P=20*log10(abs(E));
Plos=1 ./d.^2; %direct field;
P2=10*log10(abs(Plos));

plot(log10(d),P+P2,'r',log10(d),P2,'g')
xlabel('log(d)');
grid
axis([1 3 -70 -10]);

% discuss the result. Try also w low N=1 or 2

% construct the time response
% How will you plot the time response (power versus delay)?

fd = 1/2*linspace(-1,1,length(E)); %Sampling frequency divided by the scale

Et=exp(-j*k.*d);
H = fftshift(fft(E))./fftshift(fft(Et));

% Impulse response
figure()
plot(fd,abs(H).^2);

% This is the dopler shift
figure()
plot(fd,abs(fftshift(fft(E))).^2)
