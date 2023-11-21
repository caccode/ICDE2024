function [F,obj,runtime] = function(F,B,c,gamma)
% FFGC1 Fast Fuzzy Graph cluertering _ Fuzzy cmeans optimization
%   此处显示详细说明
NITR=100;
[n,m]=size(B);
% F = rand([n, c]);
% for i = 1:n
%     F(i, :) = F(i, :)./sum(F(i, :));
% end
tic
M=zeros(c,m);
for iter=1:NITR
    % update M
%     for j=1:c
%         M(j,:) = F(:,j)'*B/(F(:,j)'*ones(n,1));
%     end
    M= diag(sum(F).^(-1))*F'*B;
    % update F
    distBM=pdist2(B,M);
    for i=1:n
        hnew=-distBM(i,:)/(2*gamma);
        F(i,:)=EProjSimplex_new(hnew);
    end

    obj(iter)=trace(F'*distBM)+gamma*norm(F,'fro')^2;
    if iter>2 && abs((obj(iter)-obj(iter-1))/obj(iter))<1e-10
        break;
    end
    if iter>30 && sum(abs(obj(iter-9:iter-5)-obj(iter-5+1:iter)))<1e-10
        break;
    end
end
runtime=toc;
end

