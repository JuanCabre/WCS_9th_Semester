function r = ray_doppler(fm,M,dt,N)
T = N*dt-dt;
t = 0:dt:T
c = sqrt(2/M);
w = 2*pi*fm;
x = 0;
y = 0;

for n = 1:M
    alpha = (2*pi*n-pi+(2*pi*rand-pi))/(4*M);
    ph1 = 2*pi*rand - pi;
    ph2 = 2*pi*rand - pi;
    x = x + c*cos(w*t*cos(alpha) + ph1);
    y = y + c*cos(w*t*sin(alpha) + ph2);
end
r = sqrt(x.^2 + y.^2)/sqrt(2);