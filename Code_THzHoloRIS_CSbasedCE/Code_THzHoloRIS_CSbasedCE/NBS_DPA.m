function g = NBS_DPA(x,y,d,N,lambda,k_opt,l_opt)
    gx = diric(d*2*pi/lambda*(k_opt-sin(x)),N);    
    gy = diric(d*2*pi/lambda*(l_opt-cos(x).*sin(y)),N);
    g = gx.*gy;
end