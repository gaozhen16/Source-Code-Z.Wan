function g = NBS_CMS(x,y,A,lambda,k_opt,l_opt)
    gx = sinc(A/lambda*(k_opt-sin(x)));    
    gy = sinc(A/lambda*(l_opt-cos(x).*sin(y)));
    g = gx.*gy;
end