function [H_f_LoS, H_f_NLoS] = BS_2_IRS_channel_3D(Mx, My, Nx, Ny, K, N_cp, Ts, Lc, K_fac, AoD_ele, AoD_azi, AoA_ele, AoA_azi, tau, alpha)
    % {AoD_ele, AoD_azi, AoA_ele, AoA_azi, tau, alpha} 六个参数维度为Lc+1，第一项为LoS分量，其余为NLoS分量
    
    M = Mx * My;
    N = Nx * Ny;
    
    alpha(1) = alpha(1);
    alpha(2:end) = alpha(2:end) * sqrt(1/K_fac);
        
          
    %% spacial channel
    A_T = zeros(M,Lc+1);
    for i = 1:Lc+1
        A_T(:,i) = kron(exp(1i*pi*(0:Mx-1).'*sin(AoD_azi(i)))/sqrt(Mx),exp(1i*pi*(0:My-1).'*(sin(AoD_ele(i)).*cos(AoD_azi(i))))/sqrt(My));
    end
    
    A_R = zeros(N,Lc+1);
    for i = 1:Lc+1
        A_R(:,i) = kron(exp(1i*pi*(0:Nx-1).'*sin(AoA_azi(i)))/sqrt(Nx),exp(1i*pi*(0:Ny-1).'*(sin(AoA_ele(i)).*cos(AoA_azi(i))))/sqrt(Ny));
    end
    
    %% frequency-domain channel
    H_f_LoS = zeros(K, N, M); 
    H_f_NLoS = zeros(K, N, M);
    
    for k = 1:K
        for p = 1:Lc+1
             % beta_pk为第p条径在第k个子载波上的贡献
            beta_pk = 0;
            
            for d = 0:N_cp-1
                beta_pk = beta_pk + alpha(p) * raised_cosine(d*Ts - tau(p), 0.8, Ts) * exp(-1i*2*pi*(k-1)*d/K);
            end
            
            if(p == 1)
                H_f_LoS(k,:,:) = sqrt(M*N)* beta_pk * A_R(:,p) * A_T(:,p)';
            else
                H_f_NLoS(k,:,:) = squeeze(H_f_NLoS(k,:,:)) + beta_pk * A_R(:,p) * A_T(:,p)';
            end
        end
    end
    
    H_f_NLoS = sqrt(M*N) * H_f_NLoS;
end