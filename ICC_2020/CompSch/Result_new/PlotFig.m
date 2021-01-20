%% plot
G = [15,20:10:80];
% close all

figure
% load NMSE_SNR0_G15to80;
% plot(G,NMSE_dB,'-o','linewidth',1.6,'MarkerSize',10);

% set(gca,'fontsize',12);%'linewidth',4,'fontname','Times'

load NMSE_SNR10_G15to80;
plot(G,NMSE_dB,'-s','linewidth',1.6,'MarkerSize',8);
xlabel('\it G');ylabel('NMSE [dB]');
grid on;
hold on

load NMSE_SNR20_G15to80;
plot(G,NMSE_dB,'-v','linewidth',1.6,'MarkerSize',8);

load NMSE_SNR30_G15to80;
plot(G,NMSE_dB,'-d','linewidth',1.6,'MarkerSize',8);

legend('SNR = 10 dB',...
       'SNR = 20 dB',...
       'SNR = 30 dB')

axis tight
% \it {N}^{\rm pilot}_{T} = 28
% title('CR = 0.5, Lp = 3, N = 128');
