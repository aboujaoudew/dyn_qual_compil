%%% Commutativity of fdot operation
%%% The state domain on which this operation applies is the set of intervals
%%% which are either positive or negative

function res = fdot_2(x,y)

if (min(x) * max(x) < 0) || (min(y) * max(y) < 0)
    error('product not defined for the input values')
else
    M = [];
    % Extraction of the sign of the product
    if sum(abs(x)) * sum(abs(y)) == 0        
        sign_prod = sign(sum(x) + sum(y));
    else
        sign_prod = sign(sum(x) * sum(y));
    end
    for i = 1:length(x)
        for j = 1:length(y)
            M = [M, fdot_sing(abs(x(i)),abs(y(j)))];
        end
    end
end
a = sign_prod * [min(M):max(M)];
res = unique(a);

function res = fdot_sing(x,y)

if x*y ~= 0
    res = sign(x * y) * min(abs(x)*abs(y),3);
elseif x*y == 0 && abs(x) <= 1 && abs(y) <= 1
    res = 0;
elseif x*y == 0 && (abs(x) > 1 || abs(y) > 1)
    res = sign(x + y)*[0:abs(x + y)];
end