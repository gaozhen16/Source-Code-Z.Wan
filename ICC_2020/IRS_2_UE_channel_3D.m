function [h_f_LoS, h_f_NLoS] = IRS_2_UE_channel_3D(Nx, Ny, K, N_cp, Ts, Lc, K_fac, LoS_en, AoD_ele, AoD_azi, tau, alpha)
    % {AoD_ele, AoD_azi, tau, alpha} 四个参数维度为Lc+1，第一项为LoS分量，其余为NLoS分量
    
    N = Nx * Ny;
    
    alpha(1) = alpha(1) * LoS_en;
    alpha(2:end) = alpha(2:end) * sqrt(1/K_fac);
        
    %% NLoS
    A_T = zeros(N,Lc+1);
    for i = 1:Lc+1
        A_T(:,i) = kron(exp(1i*pi*(0:Nx-1).'*sin(AoD_azi(i)))/sqrt(Nx),exp(1i*pi*(0:Ny-1).'*(sin(AoD_ele(i)).*cos(AoD_azi(i))))/sqrt(Ny));
    end
    
    %% frequency-domain channel
    h_f_LoS = zeros(K, N);
    h_f_NLoS = zeros(K, N);
    for k = 1:K
        for p = 1:Lc+1
            % beta_pk为第p条径在第k个子载波上的贡献
            beta_pk = 0;
            
            for d = 0:N_cp-1
                beta_pk = beta_pk + alpha(p) * raised_cosine(d*Ts - tau(p), 0.8, Ts) * exp(-1i*2*pi*(k-1)*d/K);
            end
            
            if(p == 1)
                h_f_LoS(k,:) = beta_pk * A_T(:,p)';
            else
                h_f_NLoS(k,:) = h_f_NLoS(k,:) + beta_pk * A_T(:,p)';
            end
        end
    end
    
    h_f_LoS = sqrt(N) * h_f_LoS;
    h_f_NLoS = sqrt(N) * h_f_NLoS;
end