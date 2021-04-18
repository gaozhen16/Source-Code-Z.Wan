function y = raised_cosine(tau, alpha, Ts)
    x = tau/Ts;
    Num = cos(alpha*pi*x);
    Den = (1-(2*alpha*x).^2);
    DenZero = find(abs(Den)<eps);
    Result = Num ./ Den;
    Result(DenZero) = pi/4;
    y = sinc(x) .* Result;
end