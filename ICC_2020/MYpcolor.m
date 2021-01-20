function MYpcolor(A)
    [R C] = size(A);
    A = [A zeros(R,1)];
    A = [A;zeros(1,C+1)];
    pcolor(A);
    axis ij;
end