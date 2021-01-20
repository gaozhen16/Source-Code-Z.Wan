function idx = idx_DB_2dDFT(Nx, Ny)
% 找到2dDFT矩阵中对应sin俯仰角大于0的idx
% Nx Ny为2dDFT两个维度，俯仰角向下对应第二个维度中的角度大于0
% 为方便，保证输入为偶数

    idx = [];
    for i = 1:Nx
        for j = 1:Ny/2+1
            idx = union(idx,sub2ind([Nx,Ny],j,i));
        end
    end
end
% [1:Ny/2+2,Ny]
% Ny/2+3:Ny-1