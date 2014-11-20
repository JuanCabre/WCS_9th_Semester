%% MM11 Exc 1_e

clear all
close all

%Et=(0:100);
Er=(0:100);

% Energy per Transmited packet per Throughput Rate
Et = 4/3 .* Er;
plot(Er,Et,'b','LineWidth',2); hold on; grid on
% Energy per Transmited packet
Et=8/3 .* Er;
plot(Er,Et,'r','LineWidth',2);

legend('Energy per Packet per TR','Energy per Packet','Location','northwest')
title('Below the lines, for each metric UC is better than NC','FontSize',15)
xlabel('Er','FontSize',15);
ylabel('Et','FontSize',15);

figure()


k=(0:5);
plot(k,2*k+2,'--b','LineWidth',2); hold on; grid on
plot(k,5/4*k+4,'b','LineWidth',2);
plot(k,4*k+1,'--r','LineWidth',2);
plot(k,k+5,'r','LineWidth',2);
% plot(k,8*k+8,'--k','LineWidth',2);
% plot(k,5*k+16,'k','LineWidth',2);

legend('Energy per Packet UC','Energy per Packet NC','Energy per Packet per TR UC',...
    'Energy per Packet per TR NC','Overall System Energy UC','Overall System Energy NC','Location','northwest')
title('Energy per Metric vs ratio Et/Er','FontSize',15)
xlabel('Et/Er','FontSize',15);
ylabel('Energy per Metric','FontSize',15);


% 
% % Overal Ennergy
% OE_uc= 8*(Et+Er);
% OE_nc= 5*Et + 16*Er;
% 
% % Energy per Transmited packet
% Eptp_uc= OE_uc/4;
% Eptp_nc= OE_nc/4;
% 
% % Energy per Transmited packet per Throughput Rate
% TR_uc = 1/2;
% TR_nc = 4/5;
% 
% EP_ptr_uc = ;
% EP_ptr_nc = ;
% %