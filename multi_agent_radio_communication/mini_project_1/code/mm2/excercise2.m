% MM2 excercise 2 a,b
clear all
close all

%% For Djam

% make X number of uncorrelated channels
nr_channels = 8;                % dont use more than 8, error on plotting due to color function.
nr_realization = 1000;
for k = 1:nr_channels
    channel(k,:) = (randn(1,nr_realization) .* exp(-j*randn(1,nr_realization)));
end

c = get(gca,'colororder');  % get color vector
c(8,:) = [1,1,1];           % get color vector

% plot the envelope of the channels
for k = 1:nr_channels
plot(abs(channel(k,:)),'color',c(k,:)); hold on; grid on;
end
text(100,3.5,'test')
title('Envelope, two rayleigh fading channels');
xlabel('sample number');
ylabel('amplitude');


%% excercise 1.2 a

% correlated channels
p = 0;             % correlation coeffiecient
x = randn(1,1000) .* exp(-j*randn(1,1000));     % random channel 1
y = randn(1,1000) .* exp(-j*randn(1,1000));     % random channel 2
u = x;      % channel 1
for i = 1:1000
    v(i) = (sqrt(1-p^2))*y(i)+p*x(i);   % channel 2, with correlation p to channel 1
end



% check correlation coefficient, can use both functions for same result,
% corrcoeff best for matrix coefficient
check_correlation = abs(corr(u',v'))          
% abs(corrcoef(u',v'))


figure()
plot(abs(u),'color','r'); hold on; grid on;
plot(abs(v),'color','b');
text(100,3.5,'test')
title('Envelope, two rayleigh fading channels');
xlabel('sample number');
ylabel('amplitude');
legend('channel 1','channel 2')
