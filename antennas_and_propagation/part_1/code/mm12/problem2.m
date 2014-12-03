%% problem 2 

% freqeuncies in [Hz]
f1 = 1E6; 
f2 = 4E6;
f3 = 1000E6;

fc = [f1,f2,f3];

c = 3E8;
lambda = c./fc;
k = 2*pi./(lambda);




epsR = 15;
eps0 = 8.8541878176E-12;

eps = epsR*eps0;

sigma = 0.12; % [mhos] or [S].

omega = 2*pi.*[f1,f2,f3];

epsCOMP = eps*(1+sigma./(j*eps.*omega));

eps1 = epsCOMP(3);

mu = 4*pi*10^(-7);
% mu = 1;
c2 = 1./sqrt(epsCOMP*1*mu);
lambda2 = c2./fc;
k2 = 2*pi./lambda2;



theta = 0:pi/100:pi/2;


for i = 1:3 
Rv(i,:) = (epsCOMP(i)/eps0.*cos(theta)-sqrt(epsCOMP(i)/eps0-sin(theta).^2))./(epsCOMP(i)/eps0.*cos(theta)+sqrt(epsCOMP(i)/eps0-sin(theta).^2));

theta_t = k(i)/k2(i).*sin(theta);
Tv(i,:) = (1-Rv(i,:)).*cos(theta)./cos(theta_t);

Rh(i,:)= (cos(theta)-sqrt(epsCOMP(i)/eps0-sin(theta).^2))./...
    (cos(theta)+sqrt(epsCOMP(i)/eps0-sin(theta).^2));

Th(i,:) = 2.*cos(theta)./(cos(theta)+sqrt(epsCOMP(i)/eps0-sin(theta).^2));
end

thetaDeg = theta*180/pi;


fig1 = figure;
plot(thetaDeg,real(Rv(1,:)),'b',thetaDeg,imag(Rv(1,:)),'--b',...
    thetaDeg,real(Rv(2,:)),'k',thetaDeg,imag(Rv(2,:)),'--k',...
    thetaDeg,real(Rv(3,:)),'r',thetaDeg,imag(Rv(3,:)),'--r')
legend('Magnitude 1M','Phase 1M','Magnitude  4M','Phase 4M',...
    'Magnitude 100M','Phase 100M','Location','SouthWest')
title('Reflection Coefficients for vertical','FontSize',18)
xlabel('Incidence angle in degree','FontSize',13)
ylabel('Reflection Coefficient','FontSize',13)



%% 
fig2 = figure;
plot(thetaDeg,real(Tv(1,:)),'b',thetaDeg,imag(Tv(1,:)),'--b',...
    thetaDeg,real(Tv(2,:)),'k',thetaDeg,imag(Tv(2,:)),'--k',...
    thetaDeg,real(Tv(3,:)),'r',thetaDeg,imag(Tv(3,:)),'--r')
legend('Magnitude 1M','Phase 1M','Magnitude 4M','Phase 4M',...
    'Magnitude 100M','Phase 100M','Location','West')
title('Transmission Coefficients for vertical','FontSize',18)
xlabel('Incidence angle in degree','FontSize',13)
ylabel('Transmission Coefficient','FontSize',13)

%% 
fig3 = figure;
plot(thetaDeg,real(Rh(1,:)),'b',thetaDeg,imag(Rh(1,:)),'--b',...
    thetaDeg,real(Rh(2,:)),'k',thetaDeg,imag(Rh(2,:)),'--k',...
    thetaDeg,real(Rh(3,:)),'r',thetaDeg,imag(Rh(3,:)),'--r')
legend('Magnitude 1M','Phase 1M','Magnitude  4M','Phase 4M',...
    'Magnitude 100M','Phase 100M','Location','West')
title('Reflection Coefficients for horizontal','FontSize',18)
xlabel('Incidence angle in degree','FontSize',13)
ylabel('Reflection Coefficient','FontSize',13)
%%
fig4 = figure;
plot(thetaDeg,real(Th(1,:)),'b',thetaDeg,imag(Th(1,:)),'--b',...
    thetaDeg,real(Th(2,:)),'k',thetaDeg,imag(Th(2,:)),'--k',...
    thetaDeg,real(Th(3,:)),'r',thetaDeg,imag(Th(3,:)),'--r')
legend('Magnitude 1M','Phase 1M','Magnitude 4M','Phase 4M',...
    'Magnitude 100M','Phase 100M','Location','NorthEast')
title('Transmission Coefficients for horizontal','FontSize',18)
xlabel('Incidence angle in degree','FontSize',13)
ylabel('Transmission Coefficient','FontSize',13)

name3 = 'HorizontalPol_1';
name4 = 'HorizontalPol_2';
name1 = 'VerticalPol_1';
name2 = 'VerticalPol_2';
