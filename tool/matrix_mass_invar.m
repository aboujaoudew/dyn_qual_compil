%Computation of the moiety conservation matrix according to the algorithm
%described in Schuster and Hofer, J Chem Soc Faraday Trans, 1991
function res = matrix_mass_invar(S)

size_S = size(S);
n = size_S(1);
r = size_S(2);
T = [S eye(size_S(1))]; % We concatenate T with the Identity matrix

for l = 1:r % at the end of each step of the loop, a column of T is set to 0.
            %The loop stops when the first r columns of T are set to 0. 
    a = size(T);
    m = a(1);
    L = [];
    for i = 1:m
        b = find(T(i,r+1:end));
        c = [1:n];
        c(b)= [];
        L = [L {c}];
    end
    
    d = T(:,l); %Column l of T
    dd = d*d';
    
    %Determination of the indeces i for which T(i,l) = 0
    I_0 = find(d==0);
    
    %Computation of the pairs of indeces (i,k) fullfilling the 1st condition
    %T(i,l)*T(k,l)<0
    [e,f] = find(dd<0);
    g = [e,f];
    size_g = size(g);
    I_neg = [];
    for i = 1:size_g(1)
        if g(i,1) - g(i,2) > 0
            I_neg = [I_neg; g(i,:)];
        end
    end
    
    %Computation of the pairs of indeces fullfilling the 1st and 2nd
    %condition.
    %2nd condition: for all (i,k) determined above, intersect(L(i),L(k))
    %dos not belong to L(j) for all j not equal to (i,k,I_0)
    M = [];
    size_I_neg = size(I_neg);
    for i = 1:size_I_neg(1)
        a = I_neg(i,1); b = I_neg(i,2);
        c = intersect(cell2mat(L(a)),cell2mat(L(b)));
        for j = 1:m
            if j ~= a && j ~= b && ismember(j,I_0) == 0
                d = intersect(c, cell2mat(L(j)));
                if length(d) == length(c)
                    e = 0; 
                    break
                end
            else
            end
        end
        if e ~= 0
            M = [M; abs(T(b,l))*T(a,:) + abs(T(a,l))*T(b,:)];
        end
    end
    T = [M; T(I_0,:)];
end

res =  T(:,r+1:end);
