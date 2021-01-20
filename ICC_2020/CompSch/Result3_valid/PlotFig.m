%% plot
Lp = 2:6;
% close all

figure
load NMSE_DCsche_128_half_SNR20;
plot(Lp,NMSE_dB,'b--x','linewidth',1.6,'MarkerSize',10);
grid on;
% set(gca,'fontsize',12);%'linewidth',4,'fontname','Times'
hold on
xlabel('Lp');ylabel('NMSE [dB]');
xlim([2, 6]);
set(gca,'XTick',[2:6]);

load NMSE_Redu_128_half_SNR20;
plot(Lp,NMSE_dB,'r-o','linewidth',1.6,'MarkerSize',8);
hold on

legend('DC-based Support Detection Scheme',...
    'Proposed Scheme, \it G_{v} = \it G_{h} = 20')

load NMSE_DCsche_128_half_SNR10;
p = plot(Lp,NMSE_dB,'b--x','linewidth',1.6,'MarkerSize',10);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
grid on;
% set(gca,'fontsize',12);%'linewidth',4,'fontname','Times'
hold on
xlabel('Lp');ylabel('NMSE [dB]');

load NMSE_Redu_128_half_SNR10;
p = plot(Lp,NMSE_dB,'r-o','linewidth',1.6,'MarkerSize',8);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on

load NMSE_DCsche_128_half_SNR0;
p = plot(Lp,NMSE_dB,'b--x ','linewidth',1.6,'MarkerSize',10);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
grid on;
% set(gca,'fontsize',12);%'linewidth',4,'fontname','Times'
hold on
xlabel('Lp');ylabel('NMSE [dB]');

load NMSE_Redu_128_half_SNR0;
p = plot(Lp,NMSE_dB,'r-o','linewidth',1.6,'MarkerSize',8);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on

