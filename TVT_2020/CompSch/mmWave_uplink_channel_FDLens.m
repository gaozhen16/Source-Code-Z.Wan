function [H_uplink, A_BS, A_MS, alpha, ele_MS, azi_MS, ele_BS, azi_BS] = mmWave_uplink_channel_FDLens(D_MS, D_BS, ele_max_MS, azi_max_MS, ele_max_BS, azi_max_BS, Lp, flag_pure, G)
if(flag_pure)
    %% test AoA/AoD完全取自字典的格点
    idx_temp = randperm(G);
    ele_MS = asin(-sin(ele_max_MS/2)+1/G+(2*sin(ele_max_MS/2)-2/G)*(idx_temp(1:Lp)-1)/(G-1));
    idx_temp = randperm(G);
    azi_MS = asin(-sin(azi_max_MS/2)+1/G+(2*sin(azi_max_MS/2)-2/G)*(idx_temp(1:Lp)-1)/(G-1));
    idx_temp = randperm(G);
    ele_BS = asin(-sin(ele_max_BS/2)+1/G+(2*sin(ele_max_BS/2)-2/G)*(idx_temp(1:Lp)-1)/(G-1));
    idx_temp = randperm(G);
    azi_BS = asin(-sin(azi_max_BS/2)+1/G+(2*sin(azi_max_BS/2)-2/G)*(idx_temp(1:Lp)-1)/(G-1));
else
    %% random AoA/AoD,实际角度均匀分布
    ele_MS = -ele_max_MS/2+ele_max_MS*rand(1,Lp);
    azi_MS = -azi_max_MS/2+azi_max_MS*rand(1,Lp);

    ele_BS = -ele_max_BS/2+ele_max_BS*rand(1,Lp);
    azi_BS = -azi_max_BS/2+azi_max_BS*rand(1,Lp);

    %% random AoA/AoD,虚拟角度均匀分布
    % ele_MS = asin(-sin(ele_max_MS/2)+2*sin(ele_max_MS/2)*rand(1,Lp));
    % azi_MS = asin(-sin(azi_max_MS/2)+2*sin(azi_max_MS/2)*rand(1,Lp));
    % 
    % ele_BS = asin(-sin(ele_max_BS/2)+2*sin(ele_max_BS/2)*rand(1,Lp));
    % azi_BS = asin(-sin(azi_max_BS/2)+2*sin(azi_max_BS/2)*rand(1,Lp));
end

%% dimension
N_MS = MY_number_of_antennas_FDLens(D_MS, ele_max_MS, azi_max_MS);
if D_BS == 0
    N_BS = 1;
else
    N_BS = MY_number_of_antennas_FDLens(D_BS, ele_max_BS, azi_max_BS);
end

%% steering vect of MS
A_MS = zeros(N_MS,Lp);
for L = 1:Lp
    idx = 1;
    for me =  -floor(D_MS*sin(ele_max_MS/2))+1/2:floor(D_MS*sin(ele_max_MS/2))-1/2
        for ma = -floor(D_MS*sin(azi_max_MS/2)*cos(asin(me/D_MS))):floor(D_MS*sin(azi_max_MS/2)*cos(asin(me/D_MS)))
           A_MS(idx,L) = sinc(me - D_MS*sin(ele_MS(L)))*sinc(ma - D_MS*cos(ele_MS(L))*sin(azi_MS(L)));
           idx = idx + 1;
        end
    end 
    A_MS(:,L) = A_MS(:,L)/norm(A_MS(:,L),2);
end

%% steering vect of BS
A_BS = zeros(N_BS,Lp);
for L = 1:Lp
    idx = 1;
    for me =  -floor(D_BS*sin(ele_max_BS/2))+1/2:floor(D_BS*sin(ele_max_BS/2))-1/2
        for ma = -floor(D_BS*sin(azi_max_BS/2)*cos(asin(me/D_BS))):floor(D_BS*sin(azi_max_BS/2)*cos(asin(me/D_BS)))
           A_BS(idx,L) = sinc(me - D_BS*sin(ele_BS(L)))*sinc(ma - D_BS*cos(ele_BS(L))*sin(azi_BS(L)));
           idx = idx + 1;
        end
    end 
    A_BS(:,L) = A_BS(:,L)/norm(A_BS(:,L),2);
end

%% complex gain & results
alpha = sort(sqrt(N_MS*N_BS/Lp)*(randn(Lp,1)+1i*randn(Lp,1))/sqrt(2),'descend');
H_uplink = ones(1,Lp)*diag(alpha)*A_MS';
end