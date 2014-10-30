% MM2 excercise 2 a,b
clear all
close all

% excercise 1.2 a
N = 1000;
% correlated channels             
p = 0.4              % correlation coeffiecient

x = randn(1,N);     % random channel 1 (in-phase part)
y = randn(1,N);     % random channel 2 (in-phase part)
x1 = rand(1,N);     % random channel 1 (quadrature)
y1 = rand(1,N);     % random channel 2 (quadrature)

% channel 2, with correlation p to channel 1
for i = 1:N
    a(i) = (sqrt(1-p^2))*y(i)+p*x(i);    % correlated channel 2(in-phase part)
    a1(i) = (sqrt(1-p^2))*y1(i)+p*x1(i); % correlated channel 2(quadrature)
end

ch1 = x .* exp(-j*x1);
ch2 = a .* exp(-j*a1);



% check correlation coefficient, can use both functions for same result,
% corrcoeff best for matrix coefficient
check_correlation = (corr((ch1'),(ch2')))          
abs(corrcoef(ch1,ch2))

corr_value = char(sprintf('Correlation coefficient %.3f',check_correlation));

figure()  % envelope of channels
plot(abs(ch1),'color','r'); hold on; grid on;
plot(abs(ch2),'color','b');
text(25,3,corr_value);
title('Envelope, two rayleigh fading channels');
xlabel('sample number');
ylabel('amplitude');
legend('channel 1','channel 2')


% Make the CDF plot of the two channels
u_db = 10*log10(abs(ch1));    % convert to dB
v_db = 10*log10(abs(ch2));
H1 = sort(u_db);            % sort for CDF
H2 = sort(v_db);

figure()
Percent_Axis = linspace (0 ,100 , N);
semilogy(H1,Percent_Axis,'b'); hold on; grid on;
semilogy(H2,Percent_Axis,'r');
xlabel('Power [dB]')
ylabel('Log (%)')
title('CDF data plot')
legend('channel 1','channel 2','Location','East')


% excercise 1.2 b

%H = [u(1,1), v(1,1)];

% Maximum Ratio Combining
% Summing the power of the two signals, we assume a constant noise floor,
% so the noise can and is placed outside the sum.
% we are summing powers because the weights are  r*/N , and since N are the
% same for both signals, we end up summing the power of H1 and H2.
Max_ratio = (abs(ch1))+(abs(ch2));
Max_ratio_dB = 10*log10(Max_ratio);
Max_ratio_CDF = sort(Max_ratio_dB);


figure()
Percent_Axis = linspace (0 ,100 , N);
semilogy(H1,Percent_Axis,'b'); hold on; grid on;
semilogy(H2,Percent_Axis,'r');
semilogy(Max_ratio_CDF,Percent_Axis,'g');
text(-40,20,corr_value);
xlabel('Power [dB]')
ylabel('Log (%)')
title('CDF data plot')
legend('channel 1','channel 2','MRC combinig','Location','West')

