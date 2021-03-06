clear all
close all

j=sqrt(-1);

%% Excercise A

delays = [0, 0.1, 0.3, 0.5, 0.8, 1.1, 1.3, 1.7, 2.3, 3.1, 3.2, 5];
power_DB = [-4, -3, 0, -2.6, -3, -5, -7, -5, -6.5, -8.6, -11, -10];
linear = 10.^(power_DB./10);



%Experiment, ignore this
%%%%%%%%%%%%%%%%%%%%%%
% % % old_delays = [0, 0.1, 0.3, 0.5, 0.8, 1.1, 1.3, 1.7, 2.3, 3.1, 3.2, 5];
% % % old_power_DB = [-4, -3, 0, -2.6, -3, -5, -7, -5, -6.5, -8.6, -11, -10];
% % % old_linear = 10.^(old_power_DB./10);
% % % delays = [0:0.1:5];
% % % power_DB = ones(1,length(delays))*-1e17;
% % % m=1;
% % % thres= 1e-10; % Simple substractions do not result in 0 WTF!!!
% % % for k=1:length(delays)
% % %     if abs(delays(k) - old_delays(m)) < thres
% % %         power_DB(k) = old_power_DB(m);
% % %         m = m+1;
% % %     end
% % % end
% % % linear = 10.^(power_DB./10);
%%%%%%%%%%%%%%%%%%%%%%%%

delay_spread = spread(delays,linear);

%% Excercise B
for i=1:1000
    for k = 1:length(delays)
        h(i,k) = (randn() + j*randn()) * 10.^(power_DB(k)/10);
    end
end

%     stem(delays,abs(h(1,:)))
%     plot(delays,10.^(power_DB./10))

h_mean_abs = mean(abs(h));
stem(delays,h_mean_abs) %% Average Power delay profile
hold on
plot(delays,linear,'r')
delay_spread_h = spread(delays,h_mean_abs);

%% Exercise C

numbAntennas=1;

% numbLinks=size(h,1);
numbLinks=1;
channelResponse=size(h , 1);
% y=zeros(1,channelResponse)


% Aplying time reversal means convolving h(t) with h*(-t)
% where h*(-t) is the flipped conjugated of h(t)
for i=1:numbLinks
    for k=1:length(h)
        y(k,:)= conv(h(k,:),conj(flip(h(k,:))));
    end
end


% in y now we have stored the new equivalent channel impulse response which
% is the regular impulse response h(t) convolved with the time-reversal
% filter h*(-t) /// note that with "*" I mean conjugated ///

% however we need now to fix the x axis, the length of the convolved signal
% when we convolve a signal with itself should be double of the original
% signal
tr_axis = linspace(0,10,size(y,2));


% This also could be done with the following code (in frequency domain)
% however that leads to a problem due to how matlab finds the fft, we only
% end with half of the answer, uncoment the next lines triple commented to
% see the difference

% % % for i=1:numbLinks
% % %     for k=1:length(h)
% % %         Y_FF(k,:)= abs(fftshift(fft(h(k,:)))).^2;
% % %         y_ff(k,:)= ifft(Y_FF(k,:));
% % %     end
% % % end
% % % figure
% % % plot(abs(y(1,:)))
% % % hold on
% % % plot(abs(y_ff(1,:)),'r')


% Calculate the new PDP
pdp_y = mean(abs(y));
% Calculate new delay spread
delay_spread_y= spread(tr_axis,pdp_y);
% To answer Patrick's question of "How does it compare to the delay spread 
% calculated in (b)?" we can say that as expected, this signal is
% time-expanded so the delay spread will be greater.
figure
plot(tr_axis,abs(y(1,:)))

% Plotting real and imag part separately

% % symmetry_line = (length(y(1,:))+1)/2;
% % figure
% % hold on
% % plot(real(y(1,:)),'r')
% % plot(imag(y(1,:)),'b')
% % legend('real part', 'imag part')
% % y_axis_limit = get(gca,'YLim'); 
% % height = arrayfun(@(symmetry_line) line([symmetry_line symmetry_line],y_axis_limit),symmetry_line);
% % hold off

%% Exercise D

% Generate some other realizations h_2 independent of h following the same
% power delay profile this could be assumed to be impulse responses from
% another antena
for i=1:1000
    for k = 1:length(delays)
        h_2(i,k) = (randn() + j*randn()) * 10.^(power_DB(k)/10);
    end
end

%% Exercise E


% So, now, what does it means to apply time reversal in a MISO fashion? it
% means to sum the autocorrelation functions (convolutions) of each of the
% received antenas (see last slide in page 4 of mm1)

for k=1:size(h,1)
    y_MISO(k,:)= conv(h(k,:),conj(flip(h(k,:)))) + ...
                 conv(h_2(k,:),conj(flip(h_2(k,:))));
end

% Once again we need to fix the axis of the convolutions
tr_axis_MISO = linspace(0,10,size(y_MISO,2));

% Calculate the new PDP in MISO fashion
pdp_y_MISO = mean(abs(y_MISO));
% Calculate new delay spread in MISO fashion
delay_spread_y_MISO = spread(tr_axis_MISO,pdp_y_MISO);

figure
plot(tr_axis_MISO,abs(y_MISO(1,:)))

%% Extra

numb_links = 16;
% Create diferent set of impulse responses
for m=1:numb_links
    for k=1:size(h,1)
        y_MISO(k,:,m)= conv(h(k,:),conj(flip(h(k,:))));
    end
end

% Sum all the impulse responses

y_MISO=sum(y_MISO,3);

% Once again we need to fix the axis of the convolutions
tr_axis_MISO = linspace(0,10,size(y_MISO,2));

% Calculate the new PDP in MISO fashion with more than 2 antennas
pdp_y_MISO_multiple = mean(abs(y_MISO));
% Calculate new delay spread in MISO fashion with more than 2 antennas
delay_spread_y_MISO_multiple = spread(tr_axis_MISO,pdp_y_MISO_multiple);

figure
plot(tr_axis_MISO,abs(y_MISO(1,:)))
