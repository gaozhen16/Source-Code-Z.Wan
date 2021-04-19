function g = SBF(x,y,d,N,lambda,k_min,k_max,l_min,l_max)
    % x,y分别为方位、俯仰角，两者维度必须相同
    % 保证输入的向量均为列向量
    F = @(k)diric(d*2*pi/lambda*(k-x),N);
    gx = quadv(@(k)F(k),k_min,k_max);
    
    F = @(l)diric(d*2*pi/lambda*(l-y),N);
    gy = quadv(@(l)F(l),l_min,l_max);
    
    g = gx.*gy;
end