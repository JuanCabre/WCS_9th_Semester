clear all
close all

% make X number of uncorrelated channels
nr_antennas = 8;                % dont use more than 8, error on plotting due to color function.
nr_realizations = 1000;
SNR=100;
H=zeros(nr_antennas,nr_antennas,nr_realizations);

for r=1:nr_realizations
    for k = 1:nr_antennas
        H(k,:,r) = (randn(1,nr_antennas) .* exp(-1i*randn(1,nr_antennas)));
    end
end

lambda_vec=zeros(nr_realizations,1);
C_vec=zeros(nr_realizations,1);

for r=1:nr_realizations
    C=0;
    
    [U,S,V]=svd(H(:,:,r));

    lambda=diag(S);
    K=rank(S);

    for l=1:K
        C=C+log2(1+SNR*lambda(l)^2);
    end

    lambda_db=10*log(lambda(1));
    C_db=10*log(C);
    lambda_vec(r)= 2*lambda_db;
    C_vec(r)= C_db;
    
end

H1 = sort(C_vec);
H2 = sort(lambda_vec);
figure()
Percent_Axis = linspace (0 ,100 , 1000);
plot(H1,Percent_Axis,'b')
xlabel('Capacity [dB]')
ylabel('CDF (%)')
title('CDF data plot')

figure()
Percent_Axis = linspace (0 ,100 , 1000);
plot(H2,Percent_Axis,'b')
xlabel('\lambda^2 [dB]')
ylabel('CDF (%)')
title('CDF data plot')
