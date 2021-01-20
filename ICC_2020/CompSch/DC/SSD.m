function [h_out,cout] = SSD(x,Phi,L,N1,N2)
n = N1*N2;
% V1_i = 8;
% V2_i = 8;
% V = V1_i*V2_i;
V1 = 8;
V2 = 8;
% V = 2;
% R = 2*V^2+2*V+1;
x_temp = x;
support_tot = [];
for l = 1:L
    r = Phi'*x_temp;
%     r = abs(h);
    z = reshape(r,N1,N2);
    [n1 n2] = find(z==max(max(z)));

    %%% adaptive support detection
%     [V1,V2,temp,support(l,:),cout] = AS(x,Phi,n1,n2,V1,V2);
    [~,~,temp,support,cout] = AS(x_temp,Phi,n1,n2,V1,V2);
    support_tot = unique([support_tot support]);
 
%     for i = -V2/2:1:V2/2-1
%         i1 = i + V2/2 + 1;
%         for j = -V1/2:1:V1/2-1
%             j1 = j + V1/2 +1;
%             select((i1-1)*V1 + j1,1) = n1 + i;
%             select((i1-1)*V1 + j1,2) = n2 + j;
%         end
%     end
%     
%     
% %     for i = -V1/2:1:V1/2
% %         i1 = i + V1/2 + 1;
% %         for j = -V2/2:1:V2/2
% %             j1 = j + V2/2 +1;
% %             select((i1-1)*(V1+1) + j1,1) = n1 + i;
% %             select((i1-1)*(V1+1) + j1,2) = n2 + j;
% %         end
% %     end
% %     select =
% %     [n1-2,n2;n1-1,n2-1;n1-1,n2;n1-1,n2+1;n1,n2-2;n1,n2-1;n1,n2;n1,n2+1;n1,n2+2;n1+1,n2-1;n1+1,n2;n1+1,n2+1;n1+2,n2];
%     for i = 1:length(select)
%         if select(i,1) > N1
%            select(i,1) = select(i,1) - N1;
%         elseif select(i,1)<1
%            select(i,1) = select(i,1) + N1;
%         end
%     end
%     for i = 1:length(select)
%         if select(i,2) > N2
%            select(i,2) = select(i,2) - N2;
%         elseif select(i,2)<1
%            select(i,2) = select(i,2) + N2;
%         end
%     end
%     for i = 1:length(select)
%         support(l,i) = 32*(select(i,2)-1) + select(i,1);
%     end
%     Phi2 = Phi(:,unique(support(l,:)));
%     h_hat2 = inv(Phi2'*Phi2)*Phi2'*x_temp;
%     temp = zeros(n,1);
%     temp(unique(support(l,:))) = h_hat2;
    if l>=1
       x_temp = x_temp - Phi*temp;
    end
end
% h_out = temp;
% support_tot = unique(reshape(support,L*(V1)*(V2),1));
if length(support_tot)>L*(V1)*(V2) support_tot = support_tot(1:L*(V1)*(V2)); end
if length(support_tot)>size(Phi,1) support_tot = support_tot(1:size(Phi,1)); end
Phi_final = Phi(:,support_tot);
est =  inv(Phi_final'*Phi_final)*Phi_final'*x;
% est =
% inv(Phi_final'*Phi_final+sigma2*eye(length(select_final)))*Phi_final'*x;
h_out = zeros(n,1);
h_out(support_tot) = est;

