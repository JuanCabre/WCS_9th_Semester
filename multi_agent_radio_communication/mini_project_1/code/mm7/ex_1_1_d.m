clear all
close all

M=80; % Number of frequencies
K=10; % Number of networks

% Probability of network A colliding with a nother network asynchrounosly
Pc = (1/M) - (1/(M^2));

% Probability of collision with K networks

Pcol = 1 - (1-Pc)^(K-1);

throughput = (1-Pcol)*100/50e-6;