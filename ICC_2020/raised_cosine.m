function y = rasied_cosine(tau, alpha, Ts)
    x = tau/Ts;
    y = sinc(x) * cos(pi*alpha*x)/(1-4*(alpha*x)^2);
    if isinf(y)
       y = 0; 
    end
end