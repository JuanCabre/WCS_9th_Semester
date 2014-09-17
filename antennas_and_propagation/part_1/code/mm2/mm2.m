%% Template script for mini module
%%Excercise 2.43 e and f
clear all
close all
% freq = 1*10^9;
% lambda = 3*10^8/freq;
freq_range = [500:10:2000].*10^6;
lambda_range = 3*10^8./freq_range;

sigma = 5.7*10^7;
mu = 1.257*10^-6;
l = lambda_range./60;
a = lambda_range./200;
k = 2*pi./lambda_range;
R_loss = 4.415*10^-3;
R_rad = 0.21932;
Z_0 = 50;

%impedance does not change
for i = 1:length(freq_range)
    X_in(i) = -120 * ( (log(l(i)/a(i))-1) / (tan(k(i)*l(i))));
    Z_in(i) = (R_rad + R_loss) + j*X_in(i);
    Gamma(i) = (Z_in(i) - Z_0) / (Z_in(i) + Z_0); 
    SWR_1(i) = (1+abs(Gamma(i))) / (1-abs(Gamma(i)));
end
 
 % impedance changes according to formula
 for i = 1:length(freq_range)
     R_rad(i) = 80*(pi^2)*(l(i)/lambda_range(i))^2;
     R_loss(i) = l(i)/(2*pi*a(i)) * sqrt( (2*pi*freq_range(i)*mu) / (2*sigma) );
     X_in(i) = -120 * ( (log(l(i)/a(i))-1) / (tan(k(i)*l(i))));
     Z_in(i) = (R_rad(i) + R_loss(i)) + j*X_in(i);
     Gamma(i) = (Z_in(i) - Z_0) / (Z_in(i) + Z_0); 
     SWR_2(i) = (1+abs(Gamma(i))) / (1-abs(Gamma(i)));
 end
 
 % plot the figure
 figure()
 plot(freq_range,SWR_1,'r'); hold on; grid on;
 plot(freq_range,SWR_2);
 xlabel('Frequency (GHz)')
 ylabel('VSWR')
 legend('constant impedance (1GHz)','frequency dependant impedance');