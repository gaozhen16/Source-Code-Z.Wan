% »æÖÆ absAbhh 

function  PlotH( absAbhh )
maxNum = 64;
N1 = 32;
N2 = 32;

figure();
[~,maxPos] = sort(absAbhh,'descend');
maxPos2 = maxPos(1:maxNum);
polarAbhh = zeros(N1*N2,1);
polarAbhh(maxPos2) = 1;
polarAbhh2 = reshape(polarAbhh,N2,N1);
stem3(abs(polarAbhh2),'filled',':g');
hold on;
y = ceil(maxPos2/N2);
x = maxPos2 - N2 * (y - 1);
for i = 1:maxNum
    stem3(y(i),x(i),abs(polarAbhh2(x(i),y(i))),'filled',':k');
    text(y(i)+0.2,x(i),num2str(i));
end
stem3(y(1),x(1),abs(polarAbhh2(x(1),y(1)))+0.001,'filled',':r');
end

