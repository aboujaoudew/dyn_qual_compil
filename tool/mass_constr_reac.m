%% Computation of the mass constraints induced by a reaction v at a state x
%%INPUTS: x: state
%         v: a reaction
function P = mass_constr_reac(x,v,nb_val)

%Determination of the reactants and product indeces
reac_ind = []; prod_ind = [];
for i = 1: length(v)
    if v(i) < 0
        reac_ind = [reac_ind abs(v(i))];
    elseif v(i) > 0
        prod_ind = [prod_ind v(i)];        
    end
end

%Computation of the mass constraints induced by v at x
MCR_reac = [];
MCR_prod = [];
for i = 1:length(prod_ind)
    a = min(nb_val,fplus([{x{prod_ind(i)}}, {min([x{reac_ind}])}]));
    MCR_prod = [MCR_prod; prod_ind(i) a];
end

for i = 1:length(reac_ind)
    a = max(0,min(fplus([{x{reac_ind(i)}}, {-min([x{reac_ind}])}])));
    MCR_reac = [MCR_reac; reac_ind(i) a];
end

%Elimination of the redondancies (which appear for example in a
%dimerization reaction: 2 A -> B)
MCR_prod = unique(MCR_prod, 'rows');
MCR_reac = unique(MCR_reac, 'rows');

P.prod = MCR_prod;
P.reac = MCR_reac;