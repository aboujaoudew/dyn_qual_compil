%%% The state domain, D, on which this operation applies is the set of intervals
%%% which are either positive or negative
%%% test

%%% Associative in D - {0}
%%% Not associative in D. Ex: (2 x -1) x 0 not equal to 2 x (-1 x 0)

%%% Commutativity for fdot_3: comes from the commutativity for fdot_2
%%% INPUT: x : list of vectors of length > 2
function res = fdot(x)

if length(x) < 2
    error('Less than 2 vectors in the input of the dot product')
else
    %Extraction of the non-zero terms of the product
    P = [];
    for i = 1:length(x)
        if sum(abs(x{i})) ~= 0
            P = [P {x{i}}];
        else    
        end
    end
    if length(P) == 0
        res = 0;
    elseif length(P) == 1
        res = fdot_2(P{1},0);
    else
        %Recurrence to compute the product if more than 2 terms
        a = P{1};
        for i = 2:length(P)
            b = fdot_2(a,P{i});
            a = b;
        end
        if length(P) < length(x) % there is a zero in the product
            res = fdot_2(b,0);
        else
            res = b;
        end
    end
end


