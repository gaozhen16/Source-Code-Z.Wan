%% plot
Lp = -5:5:40;
% close all

figure
load NMSE_DC;
plot(Lp,NMSE_dB,'--x','linewidth',1.6,'MarkerSize',10);
grid on;
% set(gca,'fontsize',12);%'linewidth',4,'fontname','Times'
hold on
xlabel('SNR [dB]');ylabel('NMSE [dB]');

load NMSE_R15;
plot(Lp,NMSE_dB,'-s','linewidth',1.6,'MarkerSize',8);
hold on

load NMSE_R20;
plot(Lp,NMSE_dB,'-v','linewidth',1.6,'MarkerSize',8);

load NMSE_R40;
plot(Lp,NMSE_dB,'-d','linewidth',1.6,'MarkerSize',8);

load NMSE_R60_8w;
plot(Lp,NMSE_dB,'-^','linewidth',1.6,'MarkerSize',8);

load NMSE_R80_8w;
plot(Lp,NMSE_dB,'-o','linewidth',1.6,'MarkerSize',8);

legend('DC-based Support Detection Scheme [7]',...
    'Proposed Scheme, \it G_{v} = \it G_{h} = 15',...
    'Proposed Scheme, \it G_{v} = \it G_{h} = 20',...
    'Proposed Scheme, \it G_{v} = \it G_{h} = 40',...
    'Proposed Scheme, \it G_{v} = \it G_{h} = 60',...
    'Proposed Scheme, \it G_{v} = \it G_{h} = 80')
% \it {N}^{\rm pilot}_{T} = 28
% title('CR = 0.5, Lp = 3, N = 128');
