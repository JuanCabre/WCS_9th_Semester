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

d=1*lambda; %antenna spacing in ULA


%%% environment
radius1=100*lambda;
radius2=radius1;

scatter_separation= 30*lambda; % between two scatter rings

x_step = 1000;
r1 = linspace(0,20,x_step)*lambda; % movement of one array
r2 = r1; %zeros(1,x_step);

no_scatter = 100;

%%%% NOTE THIS WILL CHANGE RICHNESS FROM 'nothing' to something noticable
%%%% I.e. directly related to eigen distribution
% YOU SKETCH AND EXPLAIN, WHAT IS RICHNESS?

for o = 1:1        % for loop to make an average doppler spectrum, can be kept at 1 for speed of execution
section1 = 2*pi*rand(no_scatter,1)  % try also split up into clusters
section2 = 2*pi*rand(no_scatter,1)  % try also split up into clusters

% % cluster
% t1 = 2*pi*rand(5,1);
% t2 = 2*pi*rand(5,1);
% section1 = [rand(no_scatter/4,1)/2+t1(1); rand(no_scatter/4,1)/2+t1(2); rand(no_scatter/4,1)/2+t1(3); rand(no_scatter/4,1)/2+t1(4)];
% section2 = [rand(no_scatter/4,1)/2+t2(1); rand(no_scatter/4,1)/2+t2(2); rand(no_scatter/4,1)/2+t2(3); rand(no_scatter/4,1)/2+t2(4)];

    scatter1_xy = radius1*exp(j*section1); % scatter ring localtion, xy plane expressed as complex number
    scatter2_xy= scatter_separation + radius2*exp(j*section2);
    
    scatter1 = 1*exp(j*rand(no_scatter,1)*2*pi); % random phase
    scatter2 = 1*exp(j*rand(no_scatter,1)*2*pi); % random phase
    
    
    %%% transfer function
    % YOU CHECK
    H_test = zeros(no_scatter,no_scatter);
    h_inter = zeros(1,50);
    for x=1:length(r1);  % move one array
        for n=1:no_tx;   % tx anetnnas    centerd in scatter1
            for m=1:no_rx; % rx anetnnas, centerd in scatter2
                
                %single sum 
                h_1(n,m,x)=sum(exp(j*(r1(x)+(n-half_tx)*d)*k.*cos(angle(scatter1_xy))).*scatter1.*...
                        exp(j*(r2(x)+(m-half_tx)*d)*k*cos(angle(scatter2_xy-scatter_separation))).*scatter2); 
                    
                % The product of two separate sums (Trick to get the
                % keyhole effect)
                sum1(n,m,x) = sum(exp(j*(r1(x)+(n-half_tx)*d)*k.*cos(angle(scatter1_xy))).*scatter1);
                sum2(n,m,x) = sum(exp(j*(r2(x)+(m-half_tx)*d)*k.*cos(angle(scatter2_xy-scatter_separation))).*scatter2);
                h_2 = sum1 .* sum2;
                
                
                %             %patricks + tracing angles between two rings
                %             h_2(n,m,x)=sum(exp(j*(r1(x)+(n-half_tx)*d)*k.*cos(angle(scatter1_xy))).*scatter1.*...
                %                 exp(j*(r2(x)+(m-half_tx)*d)*k*cos(angle(scatter2_xy-scatter_separation))).*scatter2.*...
                %                 exp(j*cos(angle(scatter2_xy-scatter1_xy))));
                %
                %             % sum of sum + tracing angles between
                %             h_3(n,m,x)=sum((sum(exp(j*(r2(x)+(m-half_tx)*d)*k*cos(angle(scatter2_xy-scatter_separation))).*scatter2)).*...
                %                 exp(j*(r1(x)+(n-half_tx)*d)*k.*cos(angle(scatter1_xy))).*scatter1.*exp(j*cos(angle(scatter2_xy-scatter1_xy)))); %
                
                %             % sum of sum
                %             for k = 1:no_scatter
                %                 for u = 1:no_scatter
                %                     H_test(k,u) = exp(j*(r1(x)+(n-half_tx)*d)*k.*cos(angle(scatter1_xy(k)))).*scatter1(k).*...
                %                         exp(j*(r2(x)+(m-half_tx)*d)*k*cos(angle(scatter2_xy(u)-scatter_separation))).*scatter2(u);
                %                 end
                %             end
                %             h_inter = sum(H_test);
                %             h_inter2 = sum(h_inter);
                %             h_4(n,m,x) = h_inter2;
                
                
            end;
        end;
    end;
    
%link Doppler
h11 = abs(fftshift(fft(squeeze(h_1(1,1,:,:))))).^2;
h12 = abs(fftshift(fft(squeeze(h_1(1,2,:,:))))).^2;
h44 = abs(fftshift(fft(squeeze(h_1(no_tx,no_tx,:,:))))).^2;
h32 = abs(fftshift(fft(squeeze(h_1(half_tx+1,half_tx-1,:,:))))).^2;


h11 = h11./sum(h11);
h12 = h12./sum(h12);
h44 = h44./sum(h44);
h32 = h32./sum(h32);

dop1(o,:) = h11;
dop2(o,:) = h12;
dop3(o,:) = h44;
dop4(o,:) = h32;
dopp1 = sum(dop1,1)/o;
dopp2 = sum(dop2,1)/o;
dopp3 = sum(dop3,1)/o;
dopp4 = sum(dop4,1)/o;
end
%extract normlized doppler from movement vector
max_d = max(r1);
f = 1/max_d*x_step;
norm_fd = linspace(-f/2,f/2,x_step);

figure()
plot(norm_fd,dopp1); hold on; grid on;
plot(norm_fd,dopp2);
plot(norm_fd,dopp3);
plot(norm_fd,dopp4);
xlabel( 'Doppler [c/m]' );
ylabel( 'Probability' );
title( 'Doppler Spectrum' );
legend('h_{11}','h_{12}','h_{44}','h_{32}');



% CDF
cdf_11 = sort(20*log10 (abs(squeeze(h_1(1,1,:)))));
cdf_12 = sort(20*log10 (abs(squeeze(h_1(1,2,:)))));
cdf_44 = sort(20*log10 (abs(squeeze(h_1(no_tx,no_tx,:)))));
cdf_32 = sort(20*log10 (abs(squeeze(h_1(half_tx+1,half_tx-1,:)))));
percent_axis = linspace(0,100,x_step);

figure()
semilogy(cdf_11,percent_axis); grid on; hold on;
semilogy(cdf_12,percent_axis);
semilogy(cdf_44,percent_axis);
semilogy(cdf_32,percent_axis);
title('CDF')
xlabel('power [DB]');
ylabel('Log %');
legend('h_{11}','h_{12}','h_{44}','h_{32}');


%envelope correlation
h11_envelope = 20*log10(squeeze(h_1(1,1,:)));
h12_envelope = 20*log10(squeeze(h_1(1,2,:)));
h44_envelope = 20*log10(squeeze(h_1(no_tx,no_tx,:)));
h32_envelope = 20*log10(squeeze(h_1(half_tx+1,half_tx-1,:)));

corr1_env = corr(abs(h11_envelope),abs(h32_envelope))
corr1_pow = corr(abs(h11_envelope).^2,abs(h32_envelope).^2)
corr1_complex = corr((h11_envelope),(h32_envelope))

corr_env = sprintf('corr_{env} = %.3f',corr1_env);
corr_pow = sprintf('corr_{pow} = %.3f',corr1_pow);
% corr_complex = sprintf('corr_{complex} = %.3f',corr1_complex);

figure()
h1=plot(r1/lambda,h11_envelope,'-r',r1/lambda,h12_envelope,'-b',...);
    r1/lambda,h44_envelope,':m',r1/lambda,h32_envelope,'-.k');
text(1,-3,1,corr_env);
text(1,-5,1,corr_pow);
% text(1,-7,1,corr_complex);
title('Envelope')
xlabel('x [\lambda]');
ylabel('h [dB]');
legend(h1 ,'h_{11}','h_{12}','h_{44}','h_{32}');



%--------------------------


%singular values
% YOU CHECK
for x=1:length(r1)
    S(:,x)=svd(squeeze(h_1(:,:,x)));
end;

figure()
h2=plot(1:length(r1),20*log10(sort(S(1,:))),'-r',1:length(r1),20*log10(sort(S(2,:))),'-b',...
    1:length(r1),20*log10(sort(S(3,:))),':m',1:length(r1),20*log10(sort(S(4,:))),'-.k');
title('Modes - mixed');
xlabel('sorted distribution');
ylabel('eigen value [dB]');
legend(h2 ,'\lambda_1','\lambda_2','\lambda_3','\lambda_4');


%dual virtual (array) Doppler   ??


% plot scenario
figure()
scatter(real(scatter1_xy),imag(scatter1_xy)); hold on;
scatter(real(scatter2_xy),imag(scatter2_xy));
title('2D visual of setup')
xlabel('X [m]')
ylabel('Y [m]')

