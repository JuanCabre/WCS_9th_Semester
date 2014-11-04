% dual direction NB channels
% double ring model / variations
%
% 8/10-2006 ... 24/10-2012 Patrick Eggers

clf;
clear;
close all;

%%% radio paramter
c=3e8; %m/s
f=1e9; % Hz
lambda = c/f;
k=2*pi/lambda; %rad/m is wavenumber

%%%% antenna system
no_tx=4; %no antennas
no_rx=no_tx; %symmetric MIMO array
half_tx= round(no_tx/2);

d=0.3*lambda; %antenna spacing in ULA


%%% environment
radius1=100*lambda;
radius2=radius1;

scatter_separation=10000*lambda; % between two scatter rings

r1=(0:0.2:20)*lambda; % movement of one array
r2=r1; %zeros(1,101);

no_scatter=20;

%%%% NOTE THIS WILL CHANGE RICHNESS FROM 'nothing' to something noticable
%%%% I.e. directly related to eigen distribution
% YOU SKETCH AND EXPLAIN, WHAT IS RICHNESS?

section1=2*pi*rand(20,1); % try also split up into clusters
section2=2*pi*rand(20,1); % try also split up into clusters

scatter1_xy = radius1*exp(j*section1); % scatter ring localtion, xy plane expressed as complex number
scatter2_xy= scatter_separation + radius2*exp(j*section2);

scatter1 = 1*exp(j*rand(20,1)*2*pi); % random phase
scatter2 = 1*exp(j*rand(20,1)*2*pi); % random phase


%%% transfer function
% YOU CHECK
for x=1:length(r1);  % move one array
    for n=1:no_tx;   % tx anetnnas    centerd in scatter1
        for m=1:no_rx; % rx anetnnas, centerd in scatter2

           h(n,m,x)=sum(exp(j*(r1(x)+(n-half_tx)*d)*k.*cos(angle(scatter1_xy))).*scatter1.*...
                        exp(j*(r2(x)+(m-half_tx)*d)*k*cos(angle(scatter2_xy-scatter_separation))).*scatter2); 
                        
        % sum all scatteres to get total link signal
         % Doppler at Tx  
          % Doppler at Rx
        end;
    end;    
end;

% doppler
h11 = abs(fftshift(fft(squeeze(h(1,1,:))))).^2;
h12 = abs(fftshift(fft(squeeze(h(1,2,:))))).^2;
h44 = abs(fftshift(fft(squeeze(h(4,4,:))))).^2;
h32 = abs(fftshift(fft(squeeze(h(3,2,:))))).^2;


figure()
plot(h11); hold on;
plot(h12);
plot(h44);
plot(h32);
title('doppler')
legend('h_{11}','h_{12}','h_{44}','h_{32}');

% CDF
cdf_11 = sort(20*log10 (abs(squeeze(h(1,1,:)))));
cdf_12 = sort(20*log10 (abs(squeeze(h(1,2,:)))));
cdf_44 = sort(20*log10 (abs(squeeze(h(4,4,:)))));
cdf_32 = sort(20*log10 (abs(squeeze(h(3,2,:)))));
percent_axis = linspace(0,100,101)

figure()
semilogy(cdf_11,percent_axis); grid on; hold on;
semilogy(cdf_12,percent_axis);
semilogy(cdf_44,percent_axis);
semilogy(cdf_32,percent_axis);


figure()
h1=plot(r1/lambda,20*log10(squeeze(h(1,1,:))),'-r',r1/lambda,20*log10(squeeze(h(1,2,:))),'-b',...);
        r1/lambda,20*log10(squeeze(h(no_tx,no_tx,:))),':m',r1/lambda,20*log10(squeeze(h(half_tx+1,half_tx-1,:))),'-.k');

xlabel('x [\lambda]');
ylabel('h [dB]');
legend(h1 ,'h_{11}','h_{12}','h_{44}','h_{32}');



%--------------------------


%envelope correlation
% YOUR INPUT


%singular values
% YOU CHECK
for x=1:length(r1)
    S(:,x)=svd(squeeze(h(:,:,x)));
end;
figure()
h2=plot(1:length(r1),20*log10(sort(S(1,:))),'-r',1:length(r1),20*log10(sort(S(2,:))),'-b',...
     1:length(r1),20*log10(sort(S(3,:))),':m',1:length(r1),20*log10(sort(S(4,:))),'-.k');
xlabel('sorted distribution');
ylabel('eigen value [dB]');
legend(h2 ,'\lambda_1','\lambda_2','\lambda_3','\lambda_4');


%link Doppler
% YOUR INPUT



%dual virtual (array) Doppler   ??


for i = 1:length(scatter1_xy)
    for k = 1:length(scatter2_xy)
        links(i,k,:) = [scatter1_xy(i), scatter2_xy(k)];
    end
end

% plot scenario
figure()
scatter(real(scatter1_xy),imag(scatter1_xy)); hold on;
scatter(real(scatter2_xy),imag(scatter2_xy));
% for i = 1:length(scatter1_xy)
%     for k = 1:length(scatter2_xy)
%         plot([links(i,k,1)],[links(i,k,2)]);
%     end
% end
xlabel('lambda')
ylabel('lambda')

