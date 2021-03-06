clear;
save_flag = 1;
%% parameter setting
Lp = 6; %多径数
D_MS = 4.7;   %假设x、y维度相同
D_BS = 4.7;   %假设x、y维度相同
ele_max_MS = pi;
azi_max_MS = pi;
ele_max_BS = pi;
azi_max_BS = pi;
N_MS = MY_number_of_antennas_FDLens(D_MS, ele_max_MS, azi_max_MS);
N_BS = MY_number_of_antennas_FDLens(D_BS, ele_max_BS, azi_max_BS);
N_RF_MS = 4;
N_RF_BS = 4;
N_MS_block = N_MS/N_RF_MS;
N_BS_block = N_BS/N_RF_BS; 
N_MS_beam = 32; %N_beam须同时为N_RF和N_block的倍数
N_BS_beam = N_BS; %N_beam须同时为N_RF和N_block的倍数
G_MS_r = 20;
G_BS_r = 20;
awgn_en = 1;

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
PNR_dBs = -15:5:-10;
NMSE = zeros(1,length(PNR_dBs));

%% Channel Estimation 
% startmatlabpool(length(PNR_dBs))
for i_PNR_dB = 1:length(PNR_dBs)
    tic
    sigma2 = 10^(-(PNR_dBs(i_PNR_dB)/10));
    sigma = sqrt(sigma2);
    for iter = 1:iterMax
        tic
        H_up = mmWave_uplink_channel_FDLens(D_MS, D_BS, ele_max_MS, azi_max_MS, ele_max_BS, azi_max_BS, Lp, 0);    
        noise = sigma*(normrnd(0,1,N_BS*N_BS_beam/N_RF_BS,N_MS_beam)+1i*normrnd(0,1,N_BS*N_BS_beam/N_RF_BS,N_MS_beam))/sqrt(2);
        Y_MTC_r = W_MTC'*H_up*F_MTC + awgn_en*diag_W_MTC'*noise;
        y_MTC_r = reshape(Y_MTC_r,N_BS_beam*N_MS_beam,1);         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Original below
        Q_CSM_MTC_r = kron(F_MTC.',W_MTC');
        [H_v_CS_MTC_r, iter_num, Atom_sel] = OMP_delta(Q_CSM_MTC_r, y_MTC_r, N_BS, N_MS, sigma, Lp, 300);
        H_up_OMP_MTC_r = H_v_CS_MTC_r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Original above
        
        NMSE_temp = norm(H_up_OMP_MTC_r - H_up,'fro')^2/norm(H_up,'fro')^2;
        NMSE(i_PNR_dB) = NMSE(i_PNR_dB) + NMSE_temp;
        
        toc
        disp(['  OMPu, overhead = ' num2str(N_MS_beam) ', SNR = ' num2str(PNR_dBs(i_PNR_dB)) ', iter_max = ' num2str(iterMax) ', iter_now = ' num2str(iter) ', OMP_iter = ' num2str(iter_num) ...
            ', NMSE_now = ' num2str(10*log10(NMSE_temp)) 'dB, NMSE_avg = ' num2str(10*log10(NMSE(i_PNR_dB)/iter)) 'dB.']);
    end
    NMSE(i_PNR_dB) = NMSE(i_PNR_dB)/iterMax;
    toc
    disp(['Finished ',num2str(i_PNR_dB),'/', num2str(length(PNR_dBs)) ' , NMSE_OMP = ' num2str(10*log10(NMSE(i_PNR_dB))) 'dB']);
end
% closematlabpool
NMSE_dB = 10*log10(NMSE);
disp('Finished all');
if(save_flag)
    switch N_MS_beam
        case 16
            save NMSE_OMPu_o16_m15tom10_L6 NMSE_dB
        case 32
            save NMSE_OMPu_o32_m15tom10_L6 NMSE_dB
        case 48
            save NMSE_OMPu_o48_m15tom10_L6 NMSE_dB
    end
end
%% Plot
% figure
% plot(PNR_dBs,ASE_perf,'-m^','LineWidth',1.5); grid on;
% xlabel('PNR [dB]'),ylabel('NMSE [dB]')