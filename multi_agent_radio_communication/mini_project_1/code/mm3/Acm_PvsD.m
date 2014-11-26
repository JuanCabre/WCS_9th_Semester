close all
clear all

P=100;
Pdc=10;

d=1:0.1:30;
P1=P*d.^(-4);
P2=0.5*P1;
area=pi*d.^2;
figure
h1=plot(d,P1);
hold on
h2=plot(d,P2);
plot([1,30],[Pdc, Pdc],'-')
xlabel('Distance d [m]')
ylabel('Power at distance d')
legend([h1, h2],'Power for one antenna','P for two antennas')


figure
h1=plot(area,P1);
hold on
h2=plot(area,P2);
plot([1,900],[Pdc, Pdc],'-')
xlabel('Area [m^2]')
ylabel('Power')
legend([h1, h2],'Power for one antenna','P for two antennas')