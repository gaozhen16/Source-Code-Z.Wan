function a = a_eff(Axn, Ayn, lambda, CB_set_x, CB_set_y, AoD_azi, AoD_ele, AoA_azi, AoA_ele,type)

    Gx = length(CB_set_x);
    Gy = length(CB_set_y);
    G = Gx*Gy;
    
    if(type == 'D')
        Nx = round(2*Axn);
        Ny = round(2*Ayn);
        d = lambda / 2;
    elseif(type == 'U')
        Nx = round(4*Axn);
        Ny = round(4*Ayn);
        d = lambda / 4;
    end
    
    a = zeros(G,1);
    idx = 1;
    for gx = 1:Gx
       for gy = 1:Gy
           if(type == 'C')
                a(idx) = sinc(Axn*(CB_set_x(gx)-sin(AoD_azi)+sin(AoA_azi))) .*...
                         sinc(Ayn*(CB_set_y(gy)-sin(AoD_ele).*cos(AoD_azi)+sin(AoA_ele).*cos(AoA_azi)));
           else
                a(idx) = diric(d*2*pi/lambda*(CB_set_x(gx)-sin(AoD_azi)+sin(AoA_azi)),Nx) .*...
                         diric(d*2*pi/lambda*(CB_set_y(gy)-sin(AoD_ele).*cos(AoD_azi)+sin(AoA_ele).*cos(AoA_azi)),Ny);
           end
           idx = idx + 1;
       end
    end
end