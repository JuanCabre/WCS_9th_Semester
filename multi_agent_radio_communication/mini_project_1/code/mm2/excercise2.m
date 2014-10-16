% MM2 excercise 2 a,b
clear all
close all

% excercise 1.2 a
N = 1000;
% correlated channels
p = .8;             % correlation coeffiecient
x = randn(1,N);     % random channel 1 (in-phase part)
y = randn(1,N);     % random channel 2 (in-phase part)
x1 = randn(1,N);
y1 = randn(1,N);

a = x;      % channel 1
a1 = x1;
for i = 1:N
    b(i) = (sqrt(1-p^2))*y(i)+p*x(i);   % channel 2, with correlation p to channel 1
    b1(i) = (sqrt(1-p^2))*y1(i)+p*x1(i);
end
u = a .* exp(-j*a1);
v = b .* exp(-j*b1);

% check correlation coefficient, can use both functions for same result,
% corrcoeff best for matrix coefficient
check_correlation = (corr(abs(u'),abs(v')))          
abs(corrcoef(u',v'))

corr_value = char(sprintf('Correlation coefficient %.3f',check_correlation));

figure()  % envelope of channels
plot(abs(u),'color','r'); hold on; grid on;
plot(abs(v),'color','b');
text(25,3,corr_value);
title('Envelope, two rayleigh fading channels');
xlabel('sample number');
ylabel('amplitude');
legend('channel 1','channel 2')


% Make the CDF plot of the two channels
u_db = 10*log10(abs(u));    % convert to dB
v_db = 10*log10(abs(v));
H1 = sort(u_db);            % sort for CDF
H2 = sort(v_db);

figure()
Percent_Axis = linspace (0 ,100 , N);
semilogy(H1,Percent_Axis,'b'); hold on; grid on;
semilogy(H2,Percent_Axis,'r');
xlabel('Power [dB]')
ylabel('Log (%)')
title('CDF data plot')
legend('channel 1','channel 2')


% excercise 1.2 b

%H = [u(1,1), v(1,1)];

% Maximum Ratio Combining
% Summing the power of the two signals, we assume a constant noise floor,
% so the noise can and is placed outside the sum.
% we are summing powers because the weights are  r*/N , and since N are the
% same for both signals, we end up summing the power of H1 and H2.
Max_ratio = (abs(u))+(abs(v));
Max_ratio_dB = 10*log10(Max_ratio);
Max_ratio_sort = sort(Max_ratio_dB);


figure()
Percent_Axis = linspace (0 ,100 , N);
semilogy(H1,Percent_Axis,'b'); hold on; grid on;
semilogy(H2,Percent_Axis,'r');
semilogy(Max_ratio_sort,Percent_Axis,'g');
text(-55,50,corr_value);
xlabel('Power [dB]')
ylabel('Log (%)')
title('CDF data plot')
legend('channel 1','channel 2','MRC combinig')

