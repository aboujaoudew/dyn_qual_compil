%%% Commutativity and associativity of fplus operation
function res = fplus_2(x,y)

M = [];
for i = 1:length(x)
    for j = 1:length(y)
        M = [M, fplus_2_sing(x(i),y(j))];
    end
end

res = [min(M):max(M)]; 

% equivalent to res = [fplus_sing(min(x),min(y)), fplus_sing(max(x),max(y))]; ?

res = unique(res);
   
function res = fplus_2_sing(x,y)

if x + y ~= 0
    res = sign(x + y) * max(abs(x),abs(y));
elseif (x + y == 0) && (x * y ~= 0)
    res = [-abs(x):abs(x)];
elseif (x + y == 0) && (x * y == 0)
    res = 0;
end