function res = fplus(x)

if length(x) < 2
    error('Less than 2 vectors in the input of the dot product')
else
    %Recurrence to compute the addition if more than 2 terms
    a = x{1};
    for i = 2:length(x)
        b = fplus_2(a,x{i});
        a = b;
    end
    res = b;
end


