
% function res = fdot_2(x,y)
% 
% if (min(x) * max(x) < 0) || (min(y) * max(y) < 0)
%     error('product not defined for the input values')
% else
%     M = [];
%     Extraction of the sign of the product
%     if sum(abs(x)) * sum(abs(y)) == 0        
%         sign_prod = sign(sum(x) + sum(y));
%     else
%         sign_prod = sign(sum(x) * sum(y));
%     end
%     for i = 1:length(x)
%         for j = 1:length(y)
%             M = [M, fdot_sing(abs(x(i)),abs(y(j)))];
%         end
%     end
% end
% a = sign_prod * [min(M):max(M)];
% res = unique(a);


function res = fdot_k(x,y,k)

res = fdot(x,y) + k;


function res = fdot(x,y)

if x*y ~= 0
    res = [x + y, x + y + 1];
else
    res = [0 : max(x,y) + 1];
end