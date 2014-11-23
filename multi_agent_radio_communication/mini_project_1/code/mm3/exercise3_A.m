close all
clear all

gamma = 2;
A_c = 100; % just picked a random number
C = 1; % dont know this shit yet
M = 1:1:10;

A_cm = (A_c .* M.^-(2/gamma));
A_cmtot = A_cm .* M;

figure()
plot(M,A_cm); grid on;
title('A_{cm} [m^2] for a given number of antennas, \gamma=2');
xlabel('M');
ylabel('Area [m^2]');

figure()
plot(M,A_cmtot); grid on;
title('A_{cm,total} for a given number of antennas, \gamma=2');
xlabel('M');
ylabel('Area [m^2]');


%% 
gamma = 1:1:10;
M = 2;
x = M.^(1-2./gamma)*A_c;

figure()
plot(gamma,x); grid on;
title('A_{cm,total} for a given path loss coefficient, M=2');
xlabel('\gamma');
ylabel('Area [m^2]');

%%
gamma = 1:1:10;
M = 2;
x = M.^-(2./gamma)*A_c;

figure()
plot(gamma,x); grid on;
title('A_{cm} for a given path loss coefficient, M=2');
xlabel('\gamma');
ylabel('Area [m^2]');


%%
for gamma=5:1:25
    for M=1:1:10
        A_cmtot2(gamma-4,M)=M.^(-2./(gamma/5)).*M*A_c;
    end
end

gamma=1:0.2:5;
M=1:1:10;
figure(); grid on;
mesh(M,gamma,A_cmtot2)
title('A_{cm,total} for a given path loss coeff. and number of antennas');
ylabel('\gamma');
xlabel('M');
zlabel('Area [m^2]');