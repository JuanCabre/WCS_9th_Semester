close all
clear all

gamma = 2;
A_c = 100; % just picked a random number
C = 1; % dont know this shit yet
M = 1:1:10;

A_cm = (A_c .* M.^-(2/gamma));
A_cmtot = A_cm .* M;

pl_max = C*(A_c / pi)^(gamma/2);

figure()
plot(M,A_cm); grid on;
title('A_{cm} [m^2] for a given number of antennas');
xlabel('M');
ylabel('Area [M^2]');

figure()
plot(M,A_cmtot); grid on;
title('A_{cm,total} for a given number of antennas');
xlabel('M');
ylabel('Area [M^2]');


%% 
gamma = 1:1:10;
M = 2;
x = M.^-(gamma./2).*M*A_c;

figure()
plot(gamma,x); grid on;
xlabel('\gamma');
ylabel('Area [M^2]');