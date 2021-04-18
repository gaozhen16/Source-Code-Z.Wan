clear;clc;
plot_en = 0;
noise_en = 1;
save_en = 0;
N_redu_Ad = 2;redu_Dd_en = 0;
N_redu_Dd = 2;redu_Ad_en = 0;
IRS_type = 'C';    % C: CMS。U: ultra-dense(0.25)。D: critical(0.5)。
DSC_type = 'R';    % R: random。B: block。U: uniform。
%% Para Setting
c = 3e8;    % lightspeed [m/s]
fc = 0.15e12;  % carrier frequency [Hz]
lambda = c/fc; 

m = 64;
M = m^2;    % BS, UPA   

n = 8;
N = n^2;    % UE, UPA

An = 100;
g = 2*round(An);
A = An*lambda;  % IRS physical size

M_RF = 4;

BW = 0.5e9;   % bandwidth [Hz]
Ts = 1 / BW;
d = 0.5*lambda;   % element spacing
switch IRS_type
    case 'D'
        dI = 0.5*lambda;
        gI = round(A/dI);
    case 'U'
        dI = 0.25*lambda;
        gI = round(A/dI);
end

N_cp = 64;    % sub-carrier i.e. FFT-size
F_Ncp = dftmtx(N_cp)/sqrt(N_cp);

Lc = 1;
Lp = 1;
AS = 7.5; % angle spread [deg]

K_f_dB = 30;
K_f = 10^(K_f_dB/10);  % Rician factor

%% dictionary
K = N_redu_Dd * N_cp;
B_r = zeros(N_cp,K);
for k = 1:K
    B_r(:,k) = raised_cosine((0:N_cp-1).'*Ts - (k-1)*(N_cp-1)*Ts/K, 0.8, Ts);
end

%% Environment Para Setting
N_UE = 4;
R = 20;    
height = 10;

Loc_BS = [0;0;height];
Loc_IRS = [R;0;height];

%% DSC allocation scheme
switch DSC_type
    case 'U'
        K_set_total = zeros(1,N_cp);
        for u = 1:N_UE
            K_set_total((u-1)*N_cp/N_UE+1:u*N_cp/N_UE) = u:N_UE:(N_cp-N_UE+u);
        end
    case 'B'
        K_set_total = 1:N_cp;
end

%% Antenna gain at the UE
load D_UE_dB;

% D_UE_dB = zeros(n,n);
% for i = 1:n
%     for j = 1:n
%         k_opt = -1+2*(i-1)/n;
%         l_opt = -1+2*(j-1)/n;
%         F = @(x,y)NBS_DPA(x,y,d,n,lambda,k_opt,l_opt).^2.*cos(x);
%         D_UE_dB(i,j) = 10*log10(4*pi/integral2(F,-pi/2,pi/2,-pi/2,pi/2));
%     end
% end
% save D_UE_dB D_UE_dB

%% Antenna gain at the IRS
if(IRS_type == 'C')
    load D_IRS_C_dB;
elseif (IRS_type == 'D')
    load D_IRS_D_dB;
else
    load D_IRS_U_dB;
end

% D_IRS_dB = zeros(g,g);
% for i = 1:g
%     for j = 1:g
%         if(IRS_type == 'C')
%             F = @(x,y)NBS_CMS(x,y,A,lambda,-1+2*(i-1)/g,-1+2*(j-1)/g).^2.*cos(x);
%         else
%             F = @(x,y)NBS_DPA(x,y,dI,gI,lambda,-1+2*(i-1)/g,-1+2*(j-1)/g).^2.*cos(x);
%         end
%         D_IRS_dB(i,j) = 10*log10(4*pi/integral2(F,-pi/2,pi/2,-pi/2,pi/2,'RelTol',1e-2));
%     end
% end
% save(['D_IRS_' IRS_type '_dB'],'D_IRS_dB')

%% path loss between IRS-BS
PL_IRS2BS_dB = FS_PL_dB(fc,R);
F = @(x,y)NBS_DPA(x,y,d,m,lambda,0,0).^2.*cos(x);
D_BS_dB = 10*log10(4*pi/integral2(F,-pi/2,pi/2,-pi/2,pi/2));

%% No. of groups
b_IRS = 20;
g_IRS = g / b_IRS;

%% noise power
NSD = -174; % dBm/Hz
sigma2 = 10^(NSD / 10) * BW / 1000;
sigma = sqrt(sigma2);
noise_dBm = 10*log10(sigma2*1000);

%% uplink CE
iterMax = 100;
Ptx_dBm = 23;
Np_IRS_set = 10:20:110;

NMSE_IRS = zeros(1,length(Np_IRS_set));
NMSE_LS = zeros(1,length(Np_IRS_set));

for iter = 1:iterMax
    %% Environment  
    Loc_UE = [zeros(2,N_UE);1.5*ones(1,N_UE)];

    % UE position is fixed��near IRS��
    load Loc_UE_sample_nearIRS
    Loc_UE(1:2,1) = Loc_UExy(16,:).';

    if(plot_en)
       figure
       stem3(Loc_BS(1),Loc_BS(2),Loc_BS(3),'filled')
       hold on
       stem3(Loc_UE(1,1),Loc_UE(2,1),Loc_UE(3,1),'filled')
       hold on
       stem3(Loc_IRS(1,1),Loc_IRS(2,1),Loc_IRS(3,1),'filled')
       hold on
       legend('BS','UE','IRS')
       xlabel('x');ylabel('y');zlabel('z');
       ylim([-R,R])
    end

    %% molecular absorption given by ITU;
    PL_IRS_Abs_dB = zeros(1,N_UE);
    for i = 1:N_UE
        d_IRS = norm(Loc_UE(:,i)-Loc_IRS);
        PL_IRS_Abs_dB(i) = MoleAbs(fc)*(R+d_IRS)/1000;
    end

    %% LoS angles, consider the spatial frequency as "sin(ang_azi)+sin(ang_ele)cos(ang_azi)"
    ang_ele_IRSul_LoS = zeros(N_UE,1);
    for i = 1:N_UE
        ang_ele_IRSul_LoS(i) = atan((Loc_IRS(3)-Loc_UE(3,i)) / (Loc_IRS(1) - Loc_UE(1,i)));
    end

    ang_azi_IRSul_LoS = zeros(N_UE,1);
    for i = 1:N_UE
        d_3D = norm(Loc_IRS - Loc_UE(:,i));
        ang_azi_IRSul_LoS(i) = asin(Loc_UE(2,i) / d_3D);
    end

    ang_azi_UEIRSul_LoS = -pi/2+pi*rand(N_UE,1);
    ang_ele_UEIRSul_LoS = -pi/2+pi*rand(N_UE,1);

    %% NLoS
    ang_azi_BSul_NLoS = zeros(N_UE, Lc, Lp);
    ang_ele_BSul_NLoS = zeros(N_UE, Lc, Lp);
    for iu = 1:N_UE
        for ic = 1:Lc
            ang_azi_temp = -90 + AS + (180 - 2 * AS) * rand;
            ang_ele_temp = 60 + AS + (30 - 2 * AS) * rand;
            ang_azi_BSul_NLoS(iu,ic,:) = deg2rad([ang_azi_temp, ang_azi_temp - AS + 2 * AS * rand(1,Lp-1)]);
            ang_ele_BSul_NLoS(iu,ic,:) = deg2rad([ang_ele_temp, ang_ele_temp - AS + 2 * AS * rand(1,Lp-1)]);     
        end
    end

    ang_azi_UEBSul_NLoS = zeros(N_UE, Lc, Lp);
    ang_ele_UEBSul_NLoS = zeros(N_UE, Lc, Lp);
    for iu = 1:N_UE
        for ic = 1:Lc
            ang_azi_temp = -90 + AS + (180 - 2 * AS) * rand;
            ang_ele_temp = -90 + AS + (180 - 2 * AS) * rand;
            ang_azi_UEBSul_NLoS(iu,ic,:) = deg2rad([ang_azi_temp, ang_azi_temp - AS + 2 * AS * rand(1,Lp-1)]);
            ang_ele_UEBSul_NLoS(iu,ic,:) = deg2rad([ang_ele_temp, ang_ele_temp - AS + 2 * AS * rand(1,Lp-1)]);     
        end
    end

    ang_azi_IRSul_NLoS = zeros(N_UE, Lc, Lp);
    ang_ele_IRSul_NLoS = zeros(N_UE, Lc, Lp);
    for iu = 1:N_UE
        for ic = 1:Lc
            ang_azi_temp = -90 + AS + (180 - 2 * AS) * rand;
            ang_ele_temp = 45 + AS + (45 - 2 * AS) * rand;
            ang_azi_IRSul_NLoS(iu,ic,:) = deg2rad([ang_azi_temp, ang_azi_temp - AS + 2 * AS * rand(1,Lp-1)]);
            ang_ele_IRSul_NLoS(iu,ic,:) = deg2rad([ang_ele_temp, ang_ele_temp - AS + 2 * AS * rand(1,Lp-1)]);     
        end
    end    

    ang_azi_UEIRSul_NLoS = zeros(N_UE, Lc, Lp);
    ang_ele_UEIRSul_NLoS = zeros(N_UE, Lc, Lp);
    for iu = 1:N_UE
        for ic = 1:Lc
            ang_azi_temp = -90 + AS + (180 - 2 * AS) * rand;
            ang_ele_temp = -90 + AS + (180 - 2 * AS) * rand;
            ang_azi_UEIRSul_NLoS(iu,ic,:) = deg2rad([ang_azi_temp, ang_azi_temp - AS + 2 * AS * rand(1,Lp-1)]);
            ang_ele_UEIRSul_NLoS(iu,ic,:) = deg2rad([ang_ele_temp, ang_ele_temp - AS + 2 * AS * rand(1,Lp-1)]);     
        end
    end
       
    %% determine which groups the LoS angles belong to
    idx_UEIRS_opt = zeros(2,N_UE);
    idx_IRS_opt = zeros(2,N_UE);
    for iu = 1
        vec_temp = abs((-1+2*(0:n)/n)-sin(ang_azi_UEIRSul_LoS(iu)));
        idx_UEIRS_opt(1,iu) = find(vec_temp == min(vec_temp));
        if(idx_UEIRS_opt(1,iu) == 9)
           idx_UEIRS_opt(1,iu) = 1; 
        end
        vec_temp = abs((-1+2*(0:n)/n)-sin(ang_ele_UEIRSul_LoS(iu))*cos(ang_azi_UEIRSul_LoS(iu)));
        idx_UEIRS_opt(2,iu) = find(vec_temp == min(vec_temp));
        if(idx_UEIRS_opt(2,iu) == 9)
           idx_UEIRS_opt(2,iu) = 1; 
        end
        
        vec_temp = (-1+2*(((1:g_IRS+1)-1)*b_IRS)/g);
        for gg = 1:g_IRS
           if(sin(ang_azi_IRSul_LoS(iu)) >= vec_temp(gg) && sin(ang_azi_IRSul_LoS(iu)) < vec_temp(gg+1))
              idx_IRS_opt(1,iu) = gg;
              break
           end
        end
        for gg = 1:g_IRS
           if(sin(ang_ele_IRSul_LoS(iu))*cos(ang_azi_IRSul_LoS(iu,1)) >= vec_temp(gg) &&...
                   sin(ang_ele_IRSul_LoS(iu))*cos(ang_azi_IRSul_LoS(iu,1)) < vec_temp(gg+1))
              idx_IRS_opt(2,iu) = gg;
              break
           end
        end
    end   

    %% phase shift introduced by the channel
    alpha_hr = exp(1i*2*pi*rand(N_UE, 1+Lc*Lp));
    alpha_hr(:,1) =alpha_hr(:,1) * sqrt(K_f/(K_f+1));
    alpha_hr(:,2:end) = alpha_hr(:,2:end) * sqrt(1/K_f+1) / sqrt(Lc*Lp);

    %% Angel sets
    ang_azi_IRSul = zeros(N_UE,1+Lc*Lp);
    ang_ele_IRSul = zeros(N_UE,1+Lc*Lp);
    ang_azi_UEIRSul = zeros(N_UE,1+Lc*Lp);
    ang_ele_UEIRSul = zeros(N_UE,1+Lc*Lp);
    for iu = 1:N_UE
        ang_azi_IRSul(iu,:) = [ang_azi_IRSul_LoS(iu),reshape(ang_azi_IRSul_NLoS(iu,:,:),[1,Lc*Lp])];
        ang_ele_IRSul(iu,:) = [ang_ele_IRSul_LoS(iu),reshape(ang_ele_IRSul_NLoS(iu,:,:),[1,Lc*Lp])];
        ang_azi_UEIRSul(iu,:) = [ang_azi_UEIRSul_LoS(iu),reshape(ang_azi_UEIRSul_NLoS(iu,:,:),[1,Lc*Lp])];
        ang_ele_UEIRSul(iu,:) = [ang_ele_UEIRSul_LoS(iu),reshape(ang_ele_UEIRSul_NLoS(iu,:,:),[1,Lc*Lp])];
    end

    %% Optimal downlink grouping is assumed
    idx_UEIRS = idx_UEIRS_opt;
    idx_IRS = idx_IRS_opt;
    
    iu = 1; % choose one UE w.l.o.g
    %% angular domain
    CBx = -1+2*(idx_IRS(1,iu)-1)*b_IRS/g+(2*(0:b_IRS-1)/g);
    CBy = -1+2*(idx_IRS(2,iu)-1)*b_IRS/g+(2*(0:b_IRS-1)/g); 

    %% DSC
    if (DSC_type == 'R')
        K_set_total = randperm(N_cp);
    end
    K_set = K_set_total((iu-1)*N_cp/N_UE+1 : iu*N_cp/N_UE);
    F_pDFT = F_Ncp(:,K_set);

    % beam idx of UE
    k_opt = -1+2*(idx_UEIRS(1,iu)-1)/n;
    l_opt = -1+2*(idx_UEIRS(2,iu)-1)/n;

    % effective baseband channel at the UE
    hr_eff_UE = 10^(D_UE_dB(idx_UEIRS(1,iu),idx_UEIRS(2,iu))/20)*NBS_DPA(ang_azi_UEIRSul(iu,:),ang_ele_UEIRSul(iu,:),...
        d,n,lambda,k_opt,l_opt);

    % effective receive area
    d_UE2IRS = norm(Loc_UE(:,iu)-Loc_IRS);
    S_elm = (200e-6)^2;
    if(IRS_type == 'C')            
        S_eff_dB = 10*log10(A^2/(4*pi*d_UE2IRS^2));
%         S_eff_dB = 10*log10(A^2);
    else
        S_eff_dB = 10*log10(gI^2*S_elm/(4*pi*d_UE2IRS^2));
%         S_eff_dB = 10*log10(gI^2*S_elm);
    end

    % delay domain
    tau = sort(63*Ts*rand(1,1+Lc*Lp));
    B_hr = zeros(N_cp,1+Lc*Lp);
    B_hr(:,1) = raised_cosine((0:N_cp-1).'*Ts-tau(1), 0.8, Ts);
    for icp = 1:N_cp
        B_hr(icp,2:end) = raised_cosine((icp-1)*Ts-tau(2:end), 0.8, Ts);
    end

    % angular domain channel
    D = D_IRS_dB((idx_IRS(1,iu)-1)*b_IRS+1:idx_IRS(1,iu)*b_IRS,(idx_IRS(2,iu)-1)*b_IRS+1:idx_IRS(2,iu)*b_IRS).';
    D_vec = D(:);
    hr_Ad = zeros(b_IRS^2,1+Lc*Lp);
    hr_Ad(:,1) = 10.^(D_vec/20) .* a_eff(An, An, lambda, CBx, CBy, ...
            ang_azi_IRSul_LoS(iu),ang_ele_IRSul_LoS(iu), 0, 0, IRS_type);
    for p = 1:Lc*Lp
        hr_Ad(:,p+1) = 10.^(D_vec/20) .* a_eff(An, An, lambda, CBx, CBy,...
                  ang_azi_IRSul_NLoS(iu),ang_ele_IRSul_NLoS(iu), 0, 0, IRS_type);
    end

    % AdDd channel
    Hr_AdDd = hr_Ad * diag(alpha_hr(iu,:)) * diag(hr_eff_UE) * B_hr.';
    
    tic
    for i_Np = 1:length(Np_IRS_set)
        Np_IRS = Np_IRS_set(i_Np);
        W_IRS = exp(1i*2*pi*rand(Np_IRS,b_IRS^2))/b_IRS;
        noise = noise_en*sigma*(randn(Np_IRS,N_cp/N_UE)+1i*randn(Np_IRS,N_cp/N_UE))/sqrt(2);
        Prx_dBm = Ptx_dBm-10*log10(N_cp/N_UE)...
                 -(PL_IRS_Abs_dB(iu)+PL_IRS2BS_dB)...
                 +(S_eff_dB+D_BS_dB);
        Amp = sqrt(10^(Prx_dBm/10)/1000); 
        
        Y_n = W_IRS*(Amp*Hr_AdDd)*F_pDFT + noise;
        y_n = Y_n(:);
        Phi = kron(F_pDFT.',W_IRS);
        [h_hat,~,Atom] = S_OMP_Algorithm(y_n, Phi, -1, 20);
        Hr_AdDd_hat = reshape(h_hat,[b_IRS^2,N_cp]);           
        temp_IRS = norm(Hr_AdDd_hat-Amp*Hr_AdDd,'fro')^2/norm(Amp*Hr_AdDd,'fro')^2;
        NMSE_IRS(i_Np) =  NMSE_IRS(i_Np) + temp_IRS;

        disp(['  {' IRS_type,',' DSC_type '}, Np_IRS = ', num2str(Np_IRS) ', iter = ' num2str(iter) '/'  num2str(iterMax)...
              ', NMSE_avg = ' num2str(10*log10(NMSE_IRS(i_Np)/iter)) 'dB.']);
    end
    toc
end

%% Save data
if(save_en)
    NMSE_dB = 10*log10(NMSE_IRS/iterMax);
    save(['NMSE_dB_' IRS_type '_' DSC_type '_Np10to110'],'NMSE_dB')
end