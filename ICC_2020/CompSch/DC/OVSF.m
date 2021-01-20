function W = OVSF(n)
    W = 1;
    for i=1:log2(n)
       W = [W,W;W,-W]; 
    end
    W = W/sqrt(n);
end