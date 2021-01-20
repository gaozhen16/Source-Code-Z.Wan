clear;

w_r = (-pi/2 : 0.001 : pi/2);
N = 16;  %天线数量
F = dftmtx(N) / sqrt(N);
select_en = 1;  % 1:not rand ; 0:rand

% 生成扫描矢量
base_vector = exp(1i*pi*(0:N-1).'*sin(w_r));
if (select_en)
    target_vector_1 = exp(1i*2*pi*rand(N,1));F(:,2);
    
    s_ad = F'*target_vector_1;
    s_ad(1:N/2,:) = zeros((N/2),1);
    target_vector_2 = F*s_ad/sqrt(N);
    target_vector_3 = exp(1i*angle(F*s_ad))/sqrt(N);
else
    target_vector = exp(1i*2*pi*rand(N,1));
end
f_r_1 = abs(base_vector'*target_vector_1) / N;
f_r_2 = abs(base_vector'*target_vector_2) / N;
f_r_3 = abs(base_vector'*target_vector_3) / N;

figure
polarplot(w_r+pi/2, f_r_1,'-','linewidth',3)   %加pi/2，将图像移动到上半圆
ax = gca;
axis tight
ax.ThetaLim = [0 180];
% ax.RLim = [0 1];
ax.ThetaTickLabel = {'-90','-60','-30','0','30','60','90'};

figure
polarplot(w_r+pi/2, f_r_2,'-','linewidth',3)   %加pi/2，将图像移动到上半圆
ax = gca;
axis tight
ax.ThetaLim = [0 180];
% ax.RLim = [0 1];
ax.ThetaTickLabel = {'-90','-60','-30','0','30','60','90'};

figure
polarplot(w_r+pi/2, f_r_3,'-','linewidth',3)   %加pi/2，将图像移动到上半圆
ax = gca;
axis tight
ax.ThetaLim = [0 180];
% ax.RLim = [0 1];
ax.ThetaTickLabel = {'-90','-60','-30','0','30','60','90'};


