    clear;
    plot_en = 0;
    save_en = 1;
%% Para Setting
    c = 3e8;    % 光速
    Mx = 16;
    My = 16; 
    M = Mx * My;    % BS, UPA
    Nx = 16;
    Ny = 16;
    N = Nx * Ny;    % IRS, UPA

    fc = 30e9;  % 载频
    BW = 100e6;
    Ts = 1 / BW;
    K = 64;    % N_sub-carrier i.e. FFT-size
    delta_f = BW / K;
    N_cp = 64;
    
%% dictionary
    N_redu = 2;

    A_d_azi = dftmtx(N_redu * Mx) / sqrt(Mx);
    A_d_ele = dftmtx(N_redu * My) / sqrt(My);
    A_d = kron(A_d_azi(1:Mx,:), A_d_ele(1:My,:));

    A_r_azi = dftmtx(N_redu * Nx) / sqrt(Nx);
    A_r_ele = dftmtx(N_redu * Ny) / sqrt(Ny);
    A_r = kron(A_r_azi(1:Nx,:), A_r_ele(1:Ny,:));

    dim_d = size(A_d,2);
    dim_r = size(A_r,2);    
    
    Psi = [[A_d', zeros(N_redu^2*M,N)]; [zeros(N_redu^2*N,M), A_r']].';
    
%% NMSE vs. Ptx
    Np = 128;
    noise_en = 1;
    G_NLoS_en = 1;
    Ptx_dBms = -5:5:20;
    iterMax = 500;
    NMSE_d1 = zeros(1,length(Ptx_dBms));
    NMSE_r1 = zeros(1,length(Ptx_dBms));
    NMSE_d2 = zeros(1,length(Ptx_dBms));
    NMSE_r2 = zeros(1,length(Ptx_dBms));
    
%% Environment Para Setting, 参考3GPP毫米波信道模型 Umi-street
    N_IRS = 1;  % max 1
    N_UE = 1;
    Lc = 6;
    R = 100;    % 小区半径（基准量）
    height = 10;  % 基站高度（基准量）

    Loc_BS = [0;0;height];  % 基站位于0点，架高, UPA法方法为x正方向

    if(N_IRS == 1)
        Loc_IRS = [R;0;height];   % IRS位置，架高，UPA法方法为x负方向
    end

%     startmatlabpool(2)
    for iter = 1:iterMax
        tic
        
        Loc_UE = zeros(3,N_UE);   % UE位置，设高度为1.5m
        for i = 1:N_UE
           r = 10+(R-20)*rand;70+20*rand;20+20*rand;
           theta = -pi/3 + rand*2*pi/3;
           Loc_UE(:,i) = [r*cos(theta);r*sin(theta);1.5];
        end

        Loc_sctr = zeros(3,Lc); % 散射体位置，架高至用户高度-基站高度1/3
        for i = 1:Lc
            r = 10+(R-20)*rand;
            theta = -pi/3 + rand*2*pi/3;
            Loc_sctr(:,i) = [r*cos(theta);r*sin(theta);1.5 + rand * (height/3 - 1.5)];
        end

        if(plot_en)
           figure
           stem3(Loc_BS(1),Loc_BS(2),Loc_BS(3),'filled')
           hold on
           stem3(Loc_UE(1,:),Loc_UE(2,:),Loc_UE(3,:),'filled')
           hold on
           stem3(Loc_IRS(1,:),Loc_IRS(2,:),Loc_IRS(3,:),'filled')
           hold on
           stem3(Loc_sctr(1,:),Loc_sctr(2,:),Loc_sctr(3,:),'filled')
           legend('BS','UE','IRS','Scatters')
           xlabel('x');ylabel('y');zlabel('z');
        end

    %% 几何角度计算，采用“sin(方位角)”+“sin(俯仰角)cos(方位角)”的导向矢量形式
        ang_ele_BS_NLoS = zeros(Lc,1); % 基站处NLoS的俯仰角
        for i = 1:Lc
            ang_ele_BS_NLoS(i) = atan((Loc_BS(3)-Loc_sctr(3,i)) / Loc_sctr(1,i));
        end

        ang_azi_BS_NLoS = zeros(Lc,1); % 基站处NLoS的方位角
        for i = 1:Lc
            d_3D = norm(Loc_BS - Loc_sctr(:,i));
            ang_azi_BS_NLoS(i) = asin(Loc_sctr(2,i) / d_3D);
        end

        ang_ele_IRS_NLoS = zeros(Lc,1); % IRS处NLoS的俯仰角
        for i = 1:Lc
            ang_ele_IRS_NLoS(i) = atan((Loc_IRS(3)-Loc_sctr(3,i)) / (Loc_IRS(1) - Loc_sctr(1,i)));
        end

        ang_azi_IRS_NLoS = zeros(Lc,1); % IRS处NLoS的方位角
        for i = 1:Lc
            d_3D = norm(Loc_IRS - Loc_sctr(:,i));
            ang_azi_IRS_NLoS(i) = asin(Loc_sctr(2,i) / d_3D);
        end

        ang_ele_BS_LoS = zeros(N_UE,1); % 基站to用户LoS的俯仰角
        for i = 1:N_UE
            ang_ele_BS_LoS(i) = atan((Loc_BS(3)-Loc_UE(3,i)) / Loc_UE(1,i));
        end

        ang_azi_BS_LoS = zeros(N_UE,1); % 基站to用户LoS的方位角
        for i = 1:N_UE
            d_3D = norm(Loc_BS - Loc_UE(:,i));
            ang_azi_BS_LoS(i) = asin(Loc_UE(2,i) / d_3D);
        end

        ang_ele_IRS_LoS = zeros(N_UE,1); % IRS处LoS的俯仰角
        for i = 1:N_UE
            ang_ele_IRS_LoS(i) = atan((Loc_IRS(3)-Loc_UE(3,i)) / (Loc_IRS(1) - Loc_UE(1,i)));
        end

        ang_azi_IRS_LoS = zeros(N_UE,1); % IRS处LoS的方位角
        for i = 1:N_UE
            d_3D = norm(Loc_IRS - Loc_UE(:,i));
            ang_azi_IRS_LoS(i) = asin(Loc_UE(2,i) / d_3D);
        end

    %% 时延计算
        Delay_BS2UE_LoS = zeros(N_UE,1);
        for i = 1:N_UE
            Delay_BS2UE_LoS(i) = norm(Loc_BS - Loc_UE(:,i)) / c;
        end

        Delay_IRS2UE_LoS = zeros(N_UE,1);
        for i = 1:N_UE
            Delay_IRS2UE_LoS(i) = norm(Loc_IRS - Loc_UE(:,i)) / c;
        end

        Delay_BS2UE_NLoS = zeros(Lc,N_UE);
        for i = 1:Lc
            for j = 1:N_UE
                Delay_BS2UE_NLoS(i,j) = (norm(Loc_BS - Loc_sctr(:,i)) + norm(Loc_sctr(:,i) - Loc_UE(:,j))) / c;
            end
        end

        Delay_IRS2UE_NLoS = zeros(Lc, N_UE);
        for i = 1:Lc
            for j = 1:N_UE
                Delay_IRS2UE_NLoS(i,j) = (norm(Loc_IRS - Loc_sctr(:,i)) + norm(Loc_sctr(:,i) - Loc_UE(:,j))) / c;
            end
        end

        Delay_BS2IRS_NLoS = zeros(Lc, 1);
        for i = 1:Lc
           Delay_BS2IRS_NLoS(i) = (norm(Loc_IRS - Loc_sctr(:,i)) + norm(Loc_sctr(:,i) - Loc_BS)) / c;
        end

        Delay_BS2IRS_LoS = (norm(Loc_IRS - Loc_BS)) / c;

    %% LoS是否存在
        Hr_LoS_en = ones(N_UE, 1);
        Hd_LoS_en = zeros(N_UE, 1);
        
    %% 计算所有有效路径的绝对延时，并以此计算相对延迟
        Delay_min = zeros(N_UE,1);
        Delay_max = zeros(N_UE,1);
        
        for i = 1:N_UE
           Delay_G_valid = [Delay_BS2IRS_LoS;Delay_BS2IRS_NLoS];

           if(Hr_LoS_en(i))
               Delay_Hr_valid = [Delay_IRS2UE_LoS(i);Delay_IRS2UE_NLoS(:,i)];
           else
               Delay_Hr_valid = Delay_IRS2UE_NLoS(:,i);
           end

           if(Hd_LoS_en(i))
               Delay_Hd_valid = [Delay_BS2UE_LoS(i);Delay_BS2UE_NLoS(:,i)];
           else
               Delay_Hd_valid = Delay_BS2UE_NLoS(:,i);
           end

           % 计算级联信道时，不考虑NLoS + NLoS的级联
           Delay_cascade_min = Delay_BS2IRS_LoS + min(Delay_Hr_valid);
           Delay_cascade_max = Delay_BS2IRS_LoS + max(Delay_Hr_valid);

           Delay_min(i) = min(min(Delay_Hd_valid),Delay_cascade_min);
           Delay_max(i) = max(max(Delay_Hd_valid),Delay_cascade_max);
           disp(['  Tau_max = ' num2str((Delay_max(i)-Delay_min(i))/Ts) ' samples.']);
        end

    %% 路损(只考虑free space，不包括莱斯因子)
        PL_BS2UE_LoS = zeros(N_UE, 1);
        for i = 1:N_UE
            PL_BS2UE_LoS(i) =  32.44 + 20*log10(fc/1e6) + 20*log10(c*Delay_BS2UE_LoS(i) / 1000);
        end

        PL_BS2UE_NLoS = zeros(Lc, N_UE);
        for i = 1:Lc
           for j = 1:N_UE
               PL_BS2UE_NLoS(i,j) = 32.44 + 20*log10(fc/1e6) + 20*log10(c*Delay_BS2UE_NLoS(i,j)/1000);
           end
        end

        PL_BS2IRS_LoS = 10*log10(4*pi) + 20*log10(c*Delay_BS2IRS_LoS);32.44 + 20*log10(fc/1e6) + 20*log10(c*Delay_BS2IRS_LoS/1000);
        PL_BS2IRS_NLoS = zeros(Lc, 1);
        for i = 1:Lc
            PL_BS2IRS_NLoS(i) = 10*log10(4*pi) + 20*log10(c*Delay_BS2IRS_NLoS(i));32.44 + 20*log10(fc/1e6) + 20*log10(c*Delay_BS2IRS_NLoS(i)/1000);
        end

        PL_IRS2UE_LoS = zeros(N_UE, 1);
        for i = 1:N_UE
            PL_IRS2UE_LoS(i) = 32.44 + 20*log10(fc/1e6) + 20*log10(c*(Delay_IRS2UE_LoS(i))/1000);
        end

        PL_IRS2UE_NLoS = zeros(Lc, N_UE);
        for i = 1:Lc
            for j = 1:N_UE
                PL_IRS2UE_NLoS(i,j) = 32.44 + 20*log10(fc/1e6) + 20*log10(c*(Delay_IRS2UE_NLoS(i,j))/1000);
            end
        end

    %% 信道生成    
        Hd_f_LoS = zeros(K,N_UE,M);
        Hd_f_NLoS = zeros(K,N_UE,M);
        Hr_f_LoS = zeros(K,N_UE, N);
        Hr_f_NLoS = zeros(K, N_UE, N);
        G_f_LoS = zeros(K,N_UE,N,M);
        G_f_NLoS = zeros(K,N_UE,N,M);
        
        for i = 1:N_UE
            %% BS 2 UE
            AoD_ele = [ang_ele_BS_LoS(i);ang_ele_BS_NLoS];
            AoD_azi = [ang_azi_BS_LoS(i);ang_azi_BS_NLoS];
            tau = [Delay_BS2UE_LoS(i);Delay_BS2UE_NLoS(:,i)] - Delay_min(i);
            alpha = 10.^(-[PL_BS2UE_LoS(i);PL_BS2UE_NLoS(:,i)] ./ 20) .* exp(1i*2*pi*rand(Lc+1, 1));
            [Hd_f_LoS(:,i,:), Hd_f_NLoS(:,i,:)] = BS_2_UE_channel_3D(Mx, My, K, N_cp, Ts, Lc, 100, Hd_LoS_en(i), AoD_ele, AoD_azi, tau, alpha);

            %% IRS 2 UE
            AoD_ele = [ang_ele_IRS_LoS(i);ang_ele_IRS_NLoS];
            AoD_azi = [ang_azi_IRS_LoS(i);ang_azi_IRS_NLoS];
            tau = [Delay_IRS2UE_LoS(i);Delay_IRS2UE_NLoS(:,i)] - Delay_min(i);
            alpha = 10.^(-[PL_IRS2UE_LoS(i);PL_IRS2UE_NLoS(:,i)] ./ 20) .* exp(1i*2*pi*rand(Lc+1, 1));
            [Hr_f_LoS(:,i,:), Hr_f_NLoS(:,i,:)] = IRS_2_UE_channel_3D(Nx, Ny, K, N_cp, Ts, Lc, 100, Hr_LoS_en(i), AoD_ele, AoD_azi, tau, alpha);

            %% BS 2 IRS
            AoD_ele = [0;ang_ele_BS_NLoS];
            AoD_azi = [0;ang_azi_BS_NLoS];
            AoA_ele = [0;ang_ele_IRS_NLoS];
            AoA_azi = [0;ang_azi_IRS_NLoS];
            tau = [Delay_BS2IRS_LoS;Delay_BS2IRS_NLoS] - Delay_min(1);
            alpha = 10.^(-[PL_BS2IRS_LoS;PL_BS2IRS_NLoS] ./ 20) .* exp(1i*2*pi*rand(Lc+1, 1));
            [G_f_LoS(:,i,:,:), G_f_NLoS(:,i,:,:)] = BS_2_IRS_channel_3D(Mx, My, Nx, Ny, K, N_cp, Ts, Lc, 100, AoD_ele, AoD_azi, AoA_ele, AoA_azi, tau, alpha);
        end

        %% 总信道（为了简便，先考虑一个用户）
        G_f = zeros(K,N,M);
        Hd_f = zeros(K,M);
        Hr_f = zeros(K,N);
        for k = 1:K
            G_f(k,:,:) = squeeze(G_f_LoS(k,1,:,:)) + squeeze(G_NLoS_en * G_f_NLoS(k,1,:,:));
            Hr_f(k,:) = squeeze(Hr_f_LoS(k,1,:)).' + squeeze(Hr_f_NLoS(k,1,:)).';
            Hd_f(k,:) = squeeze(Hd_f_LoS(k,1,:)).' + squeeze(Hd_f_NLoS(k,1,:)).';
        end

    %% 导频信号      
        %% PA pilot design
        CB_BS_x_uni = -1+(0:Mx-1)*2/Mx;
        CB_BS_y_uni = -1+(0:My-1)*2/My;

        vec_temp =  CB_BS_x_uni - sin(0);
        [~,idx_azi] = min(abs(vec_temp));
        f_azi = exp(1i*pi*(0:Mx-1).' * CB_BS_x_uni(idx_azi))/sqrt(Mx);

        vec_temp =  CB_BS_y_uni - sin(0);
        [~,idx_ele] = min(abs(vec_temp));
        f_ele = exp(1i*pi*(0:My-1).' * CB_BS_y_uni(idx_ele) * cos(asin(CB_BS_x_uni(idx_azi))))/sqrt(My);

        f_PA = kron(f_azi,f_ele);
        
        %% S_PA
        S_RBPA = exp(1i*2*pi*rand(M,Np)) / sqrt(M);
        for t = 1:Np
            S_RBPA(:,t) = (S_RBPA(:,t) + f_PA) / sqrt(2);
        end
        
        %% S（可选）
        % 随机
        S_RB = exp(1i*2*pi*rand(M,Np)) / sqrt(M);
        
        S1 = S_RB;
        S2 = S_RBPA;
        
        %% IRS pilot（可选）
        % 随机
        Theta_rnd = exp(1i*2*pi*rand(N,Np));
        
        Theta1 = Theta_rnd;
        Theta2 = Theta_rnd;
        
        %% Rx signal
        y1_norm = zeros(K,Np);
        y2_norm = zeros(K,Np);
        for i = 1:Np
            for k = 1:K
                y1_norm(k,i) = [Hd_f(k,:), Hr_f(k,:)] * [S1(:,i) ; diag(Theta1(:,i)) * squeeze(G_f(1,:,:)) * S1(:,i)];
                y2_norm(k,i) = [Hd_f(k,:), Hr_f(k,:)] * [S2(:,i) ; diag(Theta2(:,i)) * squeeze(G_f(1,:,:)) * S2(:,i)];
            end
        end
        
        %% 感知矩阵
        Phi1_norm = zeros(Np, M+N);
        Phi2_norm = zeros(Np, M+N);

        for i = 1:Np
            Phi1_norm(i,:) = [S1(:,i).', (diag(Theta1(:,i)) * squeeze(G_f_LoS(1,1,:,:)) * S1(:,i)).'];
            Phi2_norm(i,:) = [S2(:,i).', (diag(Theta1(:,i)) * squeeze(G_f_LoS(1,1,:,:)) * S2(:,i)).'];
        end
        
        %% experiment iter
        for i_Ptx = 1:length(Ptx_dBms)
            %% Tx power
            Ptx_linear = 10^(Ptx_dBms(i_Ptx)/10) / 1000;  % 单位：瓦特
            Amp_linear = sqrt(Ptx_linear);
        
            %% 接收信号+观测矩阵
            y1 = Amp_linear * y1_norm;
            y2 = Amp_linear * y2_norm;
            Phi1 = Amp_linear * Phi1_norm;
            Phi2 = Amp_linear * Phi2_norm;

            %% Noise
            NSD = -174; % dBm/Hz
            sigma2 = 10^(NSD / 10) * delta_f / 1000;    % 单位瓦特
            sigma = sqrt(sigma2);
            noise = noise_en * sigma * (randn(K,Np) + 1i * randn(K,Np)) / sqrt(2);
            y1_n = y1 + noise;
            y2_n = y2 + noise;

            %% CE
            Th = sigma;            
            H_a_est1 = S_OMP_Algorithm(y1_n.', Phi1*Psi, Th^2, 100);
            H_a_est2 = S_OMP_Algorithm(y2_n.', Phi2*Psi, Th^2, 100);
            
            for k = 1:K
                Hd_a_est1 = H_a_est1(1:dim_d,k);
                Hr_a_est1 = H_a_est1(dim_d+1:end,k);
                Hd_est1 = Hd_a_est1.' * A_d';
                Hr_est1 = Hr_a_est1.' * A_r';

                NMSE_temp = norm(Hd_est1 - Hd_f(k,:),'fro')^2/norm(Hd_f(k,:),'fro')^2;
                NMSE_d1(i_Ptx) =  NMSE_d1(i_Ptx) + NMSE_temp;
                NMSE_temp = norm(Hr_est1 - Hr_f(k,:),'fro')^2/norm(Hr_f(k,:),'fro')^2;
                NMSE_r1(i_Ptx) =  NMSE_r1(i_Ptx) + NMSE_temp;
                
                Hd_a_est2 = H_a_est2(1:dim_d,k);
                Hr_a_est2 = H_a_est2(dim_d+1:end,k);
                Hd_est2 = Hd_a_est2.' * A_d';
                Hr_est2 = Hr_a_est2.' * A_r';

                NMSE_temp = norm(Hd_est2 - Hd_f(k,:),'fro')^2/norm(Hd_f(k,:),'fro')^2;
                NMSE_d2(i_Ptx) =  NMSE_d2(i_Ptx) + NMSE_temp;
                NMSE_temp = norm(Hr_est2 - Hr_f(k,:),'fro')^2/norm(Hr_f(k,:),'fro')^2;
                NMSE_r2(i_Ptx) =  NMSE_r2(i_Ptx) + NMSE_temp;
            end
            
            disp(['  r = ' num2str(Hr_LoS_en) ', d = ' num2str(Hd_LoS_en) ', CR = ' num2str((Np)/(N+M)) ', Ptx = ' num2str(Ptx_dBms(i_Ptx)) ' dBm, iter= ' num2str(iter) '/'  num2str(iterMax)...
            ', NMSE_d1 = ' num2str(10*log10(NMSE_d1(i_Ptx)/iter/K)) ' dB, NMSE_r1 = ' num2str(10*log10(NMSE_r1(i_Ptx)/iter/K))...
            'dB, NMSE_d2 = ' num2str(10*log10(NMSE_d2(i_Ptx)/iter/K)) ' dB, NMSE_r2 = ' num2str(10*log10(NMSE_r2(i_Ptx)/iter/K))]);
        end
        toc
    end
%     closematlabpool
    if(save_en)
        NMSE_dB = 10*log10(NMSE_d1/iterMax/K);
        save NMSE_d_BSrnd_500times NMSE_dB

        NMSE_dB = 10*log10(NMSE_r1/iterMax/K);
        save NMSE_r_BSrnd_500times NMSE_dB

        NMSE_dB = 10*log10(NMSE_d2/iterMax/K);
        save NMSE_d_BSPA_500times NMSE_dB

        NMSE_dB = 10*log10(NMSE_r2/iterMax/K);
        save NMSE_r_BSPA_500times NMSE_dB
    end
    

