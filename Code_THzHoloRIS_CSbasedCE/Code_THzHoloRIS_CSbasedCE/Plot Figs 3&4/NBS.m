function g = NBS(x,y,d,N,lambda,k_opt,l_opt)
    % x,y分别为方位、俯仰角，两者维度必须相同
    % 保证输入的向量均为列向量
    gx = diric(d*2*pi/lambda*(k_opt-x),N);
    gy = diric(d*2*pi/lambda*(l_opt-y),N);
    
    g = gx.*gy;
end