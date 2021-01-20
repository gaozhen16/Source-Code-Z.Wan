tic
clear;
save_flag = 0;
%% parameter setting
Lp = 64; %多径数
D_MS = 6.4;4.7;10.1;   %假设x、y维度相同
D_BS = 0;   %假设x、y维度相同
ele_max_MS = pi;
azi_max_MS = pi;
ele_max_BS = pi;
azi_max_BS = pi;
N_MS = MY_number_of_antennas_FDLens(D_MS, ele_max_MS, azi_max_MS);
N_BS = 1;
N_RF_MS = 4;
N_RF_BS = 1;
N_MS_block = N_MS/N_RF_MS;
N_BS_block = N_BS/N_RF_BS; 
N_MS_beam = 64; %N_beam须同时为N_RF和N_block的倍数
N_BS_beam = N_BS; %N_beam须同时为N_RF和N_block的倍数
G_MS_r = 15;
G_BS_r = 1;
awgn_en = 1;

%% 冗余字典矩阵，设计时应保证正交性
A_MS_D_r = zeros(N_MS,G_MS_r*G_MS_r);
idx2 = 1;
for g1 = 1:G_MS_r
   sin_ele_r = -sin(ele_max_MS/2)+1/G_MS_r+(2*sin(ele_max_MS/2)-2/G_MS_r)*(g1-1)/(G_MS_r-1);
   for g2 = 1:G_MS_r
       sin_alpha_r = -sin(azi_max_MS/2)+1/G_MS_r+(2*sin(azi_max_MS/2)-2/G_MS_r)*(g2-1)/(G_MS_r-1);
       idx1 = 1;
       for me =  -floor(D_MS*sin(ele_max_MS/2))+1/2:floor(D_MS*sin(ele_max_MS/2))-1/2
            for ma = -floor(D_MS*sin(azi_max_MS/2)*cos(asin(me/D_MS))):floor(D_MS*sin(azi_max_MS/2)*cos(asin(me/D_MS)))
               A_MS_D_r(idx1,idx2) = sinc(me - D_MS*sin_ele_r)*sinc(ma - D_MS*cos(asin(sin_ele_r))*sin_alpha_r);
               idx1 = idx1 + 1;
            end
       end 
       A_MS_D_r(:,idx2) = A_MS_D_r(:,idx2)/norm(A_MS_D_r(:,idx2));
       idx2 = idx2 + 1;
   end
end
dim_of_dict_MS = size(A_MS_D_r,2);

% A_BS_D_r = zeros(N_BS,G_BS_r*G_BS_r);
% idx2 = 1;
% for g1 = 1:G_BS_r
%    sin_ele_r = -sin(ele_max_BS/2)+1/G_BS_r+(2*sin(ele_max_BS/2)-2/G_BS_r)*(g1-1)/(G_BS_r-1);
%    for g2 = 1:G_BS_r
%        sin_alpha_r = -sin(azi_max_BS/2)+1/G_BS_r+(2*sin(azi_max_BS/2)-2/G_BS_r)*(g2-1)/(G_BS_r-1);
%        idx1 = 1;
%        for me =  -floor(D_BS*sin(ele_max_BS/2))+1/2:floor(D_BS*sin(ele_max_BS/2))-1/2
%             for ma = -floor(D_BS*sin(azi_max_BS/2)*cos(asin(me/D_BS))):floor(D_BS*sin(azi_max_BS/2)*cos(asin(me/D_BS)))
%                A_BS_D_r(idx1,idx2) = sinc(me - D_BS*sin_ele_r)*sinc(ma - D_BS*cos(asin(sin_ele_r))*sin_alpha_r);
%                idx1 = idx1 + 1;
%             end
%        end 
%        A_BS_D_r(:,idx2) = A_BS_D_r(:,idx2)/norm(A_BS_D_r(:,idx2));
%        idx2 = idx2 + 1;
%    end
% end
% dim_of_dict_BS = size(A_BS_D_r,2);

%% 混合预编码码本设计，保证总相关度最小（MTC）
F_RF_temp = eye(N_MS);
W_RF_temp = eye(N_BS);
index_F = randperm(N_MS);    % random permutation
index_W = randperm(N_BS);
F_RF = F_RF_temp(:,index_F);
W_RF = W_RF_temp(:,index_W);

F_BB_q = dftmtx(N_RF_MS);
W_BB_q = dftmtx(N_RF_BS);
F_BB_q = F_BB_q(:,1:N_MS_beam/N_MS_block);
W_BB_q = W_BB_q(:,1:N_BS_beam/N_BS_block);
F_BB_MTC = kron(eye(N_MS_block),F_BB_q);
W_BB_MTC = kron(eye(N_BS_block),W_BB_q);

F_MTC = F_RF*F_BB_MTC;
F_MTC = sqrt(N_MS_beam)/norm(F_MTC,'fro')*F_MTC;
W_MTC = W_RF*W_BB_MTC;
W_MTC = sqrt(N_BS_beam)/norm(W_MTC,'fro')*W_MTC;

%% 用于加噪声
diag_W_MTC = zeros(N_BS*N_BS_beam/N_RF_BS,N_BS_beam);
for nn = 1:N_BS_beam/N_RF_BS
    diag_W_MTC((nn-1)*N_BS+1:nn*N_BS,(nn-1)*N_RF_BS+1:nn*N_RF_BS) = W_MTC(:,(nn-1)*N_RF_BS+1:nn*N_RF_BS);
end

%% 信噪比和预置空间
iterMax = 500;
PNR_dBs = -5:5:40;
NMSE = zeros(1,length(PNR_dBs));

%% Channel Estimation 
tic
% startmatlabpool(length(PNR_dBs))
for i_PNR_dB = 1:length(PNR_dBs)
%     tic
    sigma2 = 10^(-(PNR_dBs(i_PNR_dB)/10));
    sigma = sqrt(sigma2);
    for iter = 1:iterMax
%         tic
        [H_up,~,A_MS] = mmWave_uplink_channel_FDLens(D_MS, D_BS, ele_max_MS, azi_max_MS, ele_max_BS, azi_max_BS, Lp, 0);    
        noise = sigma*(normrnd(0,1,N_BS*N_BS_beam/N_RF_BS,N_MS_beam)+1i*normrnd(0,1,N_BS*N_BS_beam/N_RF_BS,N_MS_beam))/sqrt(2);
        Y_MTC_r = W_MTC'*H_up*F_MTC + awgn_en*diag_W_MTC'*noise;
        y_MTC_r = reshape(Y_MTC_r,N_BS_beam*N_MS_beam,1);         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Original below
        Q_CSM_MTC_r = kron(F_MTC.'*conj(A_MS),W_MTC');
        H_v_CS_MTC_r = Q_CSM_MTC_r \ y_MTC_r;
        H_up_OMP_MTC_r = H_v_CS_MTC_r.'*A_MS';
        
%         Q_CSM_MTC_r = kron(F_MTC.',W_MTC');
%         H_v_CS_MTC_r = OMP_delta(Q_CSM_MTC_r, y_MTC_r, N_BS, N_MS, sigma, Lp);
%         H_up_OMP_MTC_r = H_v_CS_MTC_r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Original above
        
        NMSE_temp = norm(H_up_OMP_MTC_r - H_up,'fro')^2/norm(H_up,'fro')^2;
        NMSE(i_PNR_dB) = NMSE(i_PNR_dB) + NMSE_temp;
        
%         toc
        disp(['  orLS, overhead = ' num2str(N_MS_beam) ', SNR = ' num2str(PNR_dBs(i_PNR_dB)) ', iter_max = ' num2str(iterMax) ', iter_now = ' num2str(iter) ...
            ', NMSE_now = ' num2str(10*log10(NMSE_temp)) 'dB, NMSE_avg = ' num2str(10*log10(NMSE(i_PNR_dB)/iter)) 'dB.']);
    end
    NMSE(i_PNR_dB) = NMSE(i_PNR_dB)/iterMax;
%     toc
    disp(['Finished ',num2str(i_PNR_dB),'/', num2str(length(PNR_dBs)) ' , NMSE_OMP = ' num2str(10*log10(NMSE(i_PNR_dB))) 'dB']);
end
toc
% closematlabpool
NMSE_dB = 10*log10(NMSE);
disp('Finished all');
if(save_flag)
    save NMSE_R15 NMSE_dB
end
%% Plot
figure
plot(PNR_dBs,NMSE_dB,'-m^','LineWidth',1.5); grid on;
xlabel('PNR [dB]'),ylabel('NMSE [dB]')