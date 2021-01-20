clear;
flag_save = 1;
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
Ns = 2;
N_RF_MS_used = N_RF_MS;
N_RF_BS_used = N_RF_BS;
N_MS_block = N_MS/N_RF_MS;
N_BS_block = N_BS/N_RF_BS; 
N_MS_beam = 32; %N_beam须同时为N_RF和N_block的倍数
N_BS_beam = N_BS; %N_beam须同时为N_RF和N_block的倍数
G_MS_r = 20;
G_BS_r = 20;
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
dimension_of_dict_MS = size(A_MS_D_r,2);

A_BS_D_r = zeros(N_BS,G_BS_r*G_BS_r);
idx2 = 1;
for g1 = 1:G_BS_r
   sin_ele_r = -sin(ele_max_BS/2)+1/G_BS_r+(2*sin(ele_max_BS/2)-2/G_BS_r)*(g1-1)/(G_BS_r-1);
   for g2 = 1:G_BS_r
       sin_alpha_r = -sin(azi_max_BS/2)+1/G_BS_r+(2*sin(azi_max_BS/2)-2/G_BS_r)*(g2-1)/(G_BS_r-1);
       idx1 = 1;
       for me =  -floor(D_BS*sin(ele_max_BS/2))+1/2:floor(D_BS*sin(ele_max_BS/2))-1/2
            for ma = -floor(D_BS*sin(azi_max_BS/2)*cos(asin(me/D_BS))):floor(D_BS*sin(azi_max_BS/2)*cos(asin(me/D_BS)))
               A_BS_D_r(idx1,idx2) = sinc(me - D_BS*sin_ele_r)*sinc(ma - D_BS*cos(asin(sin_ele_r))*sin_alpha_r);
               idx1 = idx1 + 1;
            end
       end 
       A_BS_D_r(:,idx2) = A_BS_D_r(:,idx2)/norm(A_BS_D_r(:,idx2));
       idx2 = idx2 + 1;
   end
end
dimension_of_dict_BS = size(A_BS_D_r,2);

%% 混合预编码码本设计，保证总相关度最小（MTC）
F_RF_temp = eye(N_MS);
W_RF_temp = eye(N_BS);
index_F = randperm(N_MS);    % random permutation
index_W = randperm(N_BS);
F_RF = F_RF_temp(:,index_F);
W_RF = W_RF_temp(:,index_W);

dftmtx_F = dftmtx(N_RF_MS);
dftmtx_W = dftmtx(N_RF_BS);
F_BB_q = dftmtx_F(:,1:N_MS_beam/N_MS_block);
W_BB_q = dftmtx_W;
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
PNR_dBs = -5;
NMSE = zeros(1,length(PNR_dBs));
BER = zeros(1,length(PNR_dBs));
N_bit_per_sym = 6;   % bits number of per symbol
M = 2^N_bit_per_sym; %调制阶数
N_frame = 250000;

%% Channel Estimation and BER
for iter = 1:iterMax
    H_up = mmWave_uplink_channel_FDLens(D_MS, D_BS, ele_max_MS, azi_max_MS, ele_max_BS, azi_max_BS, Lp, 1);  
    for i_PNR_dB = 1:length(PNR_dBs)
        tic
        sigma2 = 10^(-(PNR_dBs(i_PNR_dB)/10));
        sigma = sqrt(sigma2);
        N_error_bit = 0;

        %% Channel Estimation 
        noise = sigma*(normrnd(0,1,N_BS*N_BS_beam/N_RF_BS,N_MS_beam)+1i*normrnd(0,1,N_BS*N_BS_beam/N_RF_BS,N_MS_beam))/sqrt(2);
        Y_MTC_r = W_MTC'*H_up*F_MTC + awgn_en*diag_W_MTC'*noise;
        y_MTC_r = reshape(Y_MTC_r,N_BS_beam*N_MS_beam,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Original below
        Q_CSM_MTC_r = kron(F_MTC.'*conj(A_MS_D_r),W_MTC'*A_BS_D_r);
        [H_v_CS_MTC_r, iter_num, Atom_sel] = OMP_delta(Q_CSM_MTC_r, y_MTC_r, dimension_of_dict_BS, dimension_of_dict_MS, sigma, Lp);
        H_up_OMP_MTC_r = A_BS_D_r*H_v_CS_MTC_r *A_MS_D_r';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Original above
        
        NMSE_temp = norm(H_up_OMP_MTC_r - H_up,'fro')^2/norm(H_up,'fro')^2;
        NMSE(i_PNR_dB) = NMSE(i_PNR_dB) + NMSE_temp;
        
        %% find beam and select antenna
        h_est_temp = reshape(H_up_OMP_MTC_r,N_MS*N_BS,1);
        [~,idx_h] = sort(abs(h_est_temp),'descend');
        idx_sel_ant_MS = [];
        idx_sel_ant_BS = [];
        cnt = 1;
        while(length(idx_sel_ant_MS)~=N_RF_MS_used || length(idx_sel_ant_BS)~=N_RF_BS_used)
            [I,J]=ind2sub(size(H_up_OMP_MTC_r),idx_h(cnt));
            cnt = cnt + 1;
            if(length(idx_sel_ant_BS) < N_RF_BS_used), idx_sel_ant_BS = union(idx_sel_ant_BS,I); end
            if(length(idx_sel_ant_MS) < N_RF_MS_used), idx_sel_ant_MS = union(idx_sel_ant_MS,J); end
        end
        
       %% Reduced MIMO
        H_up_reduced = H_up_OMP_MTC_r(idx_sel_ant_BS,idx_sel_ant_MS);    
        
       %% precoding
        [U, ~, V] = svd(H_up_reduced);
        W = U(:,1:Ns); F = V(:,1:Ns);

       %% BER
        % Create a turbo encoder and decoder pair, where the interleaver indices
        % are supplied by an input argument to the |step| function.
        hTEnc = comm.TurboEncoder('InterleaverIndicesSource','Input port');
        hTDec = comm.TurboDecoder('InterleaverIndicesSource','Input port', ...
            'NumIterations',4);       
        hMod = comm.RectangularQAMModulator('ModulationOrder',M, ...
            'BitInput',true, ...
            'NormalizationMethod','Average power');
        hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M, ...
                'BitOutput',true, ...
                'NormalizationMethod','Average power', ...
                'DecisionMethod','Log-likelihood ratio', ...
                'Variance',sigma^2/N_bit_per_sym); 
        
        % Determine the interleaver indices given the frame length
        intrlvrIndices = randperm(N_frame*N_bit_per_sym);
    
        for ns = 1:Ns   
            % Generate random binary data
            data(ns,:) = randi([0 1],1,N_frame*N_bit_per_sym);
            
            % Turbo encode the data
            encodedData(ns,:) = (step(hTEnc,data(ns,:).',intrlvrIndices)).'; 
            
            % Modulate the encoded data
            modSignal(ns,:) = (step(hMod,(encodedData(ns,:)).')).';
        end
        
        y = W'*(H_up(idx_sel_ant_BS, idx_sel_ant_MS) * F * modSignal + ...
            sigma/sqrt(2)*(normrnd(0,1,N_RF_BS_used,size(modSignal,2)) + 1i*normrnd(0,1,N_RF_BS_used,size(modSignal,2))));

        receivedSignal = (W'*H_up_OMP_MTC_r(idx_sel_ant_BS, idx_sel_ant_MS) * F) \ y;
        
        for ns = 1:Ns  
            % Demodulate the received signal
            demodSignal(ns,:) = step(hDemod,receivedSignal(ns,:).').';
            
            % Turbo decode the demodulated signal. Because the bit mapping from the
            % demodulator is opposite that expected by the turbo decoder, the
            % decoder input must use the inverse of demodulated signal.
            receivedBits(ns,:) = step(hTDec,-demodSignal(ns,:).',intrlvrIndices).';
        end
        
        Tx_bit = reshape(data, N_frame*Ns*N_bit_per_sym, 1);
        Rx_bit = reshape(receivedBits, N_frame*Ns*N_bit_per_sym, 1);
        
        N_error_bit = sum(Tx_bit ~= Rx_bit) + N_error_bit;        
        BER(i_PNR_dB) = BER(i_PNR_dB) + N_error_bit;
        
        toc
        disp(['  MTC ' num2str(Ns) num2str(N_RF_MS_used) num2str(N_RF_BS_used) ', ' num2str(M) '-QAM, overhead = ' num2str(N_MS_beam) ', SNR = ' num2str(PNR_dBs(i_PNR_dB))...
            ', iter_max = ' num2str(iterMax) ', iter_now = ' num2str(iter) ', NMSE_now = ' num2str(10*log10(NMSE_temp)) 'dB, NMSE_avg = ' num2str(10*log10(NMSE(i_PNR_dB)/iter))...
            'dB, BER_temp = ' num2str(BER(i_PNR_dB)/(iter*N_frame*Ns*N_bit_per_sym)) ', bit_temp = ' num2str(iter*N_frame*Ns*N_bit_per_sym)...
            ', error_bit = ' num2str(N_error_bit)]);
    end
end
NMSE = NMSE/iterMax;
BER = BER/(iterMax*N_frame*Ns*N_bit_per_sym);
NMSE_dB = 10*log10(NMSE);
disp('Finished all');
if(flag_save)
    switch N_MS_beam
    case 16
        save BER_64QAM_MTC_o16_m5_L6_turbo BER
    case 32
        save BER_64QAM_MTC_o32_m5_L6_turbo BER
    case 48
        save BER_64QAM_MTC_o48_m5_L6_turbo BER
    end
end
%% Plot
% figure
% semilogy(PNR_dBs,BER,'-m^','LineWidth',1.5); grid on;
% xlabel('SNR [dB]'),ylabel('BER')