%%% In this version, we implemented the updating scheme by scheduling the
%%% reactions according to their time scale.
%%% We take into account the constraints from the mass invariants and from 
%%% the borns on the derivative (i.e. limiting reactant constraint)
%%% The updating policy is asynchronous.

%%%INPUTS
%%% v: vector of reactions
%%% ci: vector of initial conditions
%%% MMI_sign: matrix of mass invariants

function res = stg_discr_v1_interp(ci,MMI_sign,v)

%%Computation of the mass invariant
MI = mass_invar(ci,MMI_sign);

a = ci;
%Computation of the discrete propensity for each reaction
b = next_states(a,v); % 1st column: next states - 2nd column: associated propensities
size_b = size(b);
ns = []; % matrix of next states
if length(b) == 0
    FP = [FP; a]; %a has no successor so it is a fixed point
else
    c = cell2mat(b(:,1)); %extraction of the state from b
    c = num2cell(c);
    for j = 1:size_b(1)
        cc = invar_check(MI, MMI_sign, c(j,:));
        if cc == 0 %b is not in the invariant
        else
            ns = [ns ; b(j,:)];
        end
    end
end

res.ns = ns;
res.a = a;

%save_report(ci,res);
%dot_file(ci,S,P);

function N = invar_check(MI, MMI_sign, x) 
% N = 0: x is not in the invariant; N = 1: x is in the invariant
N = 1;
a = mass_invar(x,MMI_sign);
for i = 1:length(MI)
    b = intersect(a{i},MI{i});
    if length(b) ~= 0 % at least one element of a{i} is contained in MI{i}
    else
        N = 0;
        break        
    end
end

function N = next_states(x,v)
    
%%Computation of the propensities of the reactions
M = prop_comp(x,v); % 1st column: reaction index - 2nd column: associated propensity

%Computation of the next states of state y: we consider an asynchronous
%policy between the variables of a reaction. The transitions are labeled with the
%propensity of the reaction
N = [];
for i = 1: length(v) % loop on the reactions
    %Computation of the mass constraints induced by the reaction v{i}
    mcr = mass_constr_reac(x,v{i});
    MCR_prod = mcr.prod; % MCR_prod(:,1): index of the product variables - MCR_prod(:,2): associated born
    MCR_reac = mcr.reac; % MCR_reac(:,1): index of the reactant variables - MCR_reac(:,2): associated born

    %Reactants and product indeces
    reac_ind = MCR_reac(:,1);
    prod_ind = MCR_prod(:,1);

    y = cell2mat(x);
    %Asynchronous updating between the variables of a reaction
    for j = 1:length(reac_ind)
        c = y;
        if y(reac_ind(j)) ~= 0 % constraint on the minimum value of the variables (=0)
            c(reac_ind(j)) = y(reac_ind(j)) - 1;
            if c(reac_ind(j)) < MCR_reac(j,2) % constraint on the maximum consumption of the reactants 
            else
                N = [N; {c} M{i,2}];
            end
        end
    end

    for j = 1:length(prod_ind)
        c = y;
        if y(prod_ind(j)) ~= 2 % constraint on the maximum value of the variables (=2)
            c(prod_ind(j)) = y(prod_ind(j)) + 1;
            if c(prod_ind(j)) > MCR_prod(j,2) % constraint on the maximum production of the products
            else
                N = [N; {c} M{i,2}];
            end
        end
    end    
end

function M = prop_comp(x,v)

% Computation of the propensities of the reactions v(i) at state x
M = [];
for i = 1:length(v)
ind = [];
a = v{i}; % reaction i
    for j = 1: length(a)
        if a(j) < 0
            ind = [ind abs(a(j))]; %list of the reactant indeces of reaction i
        end
    end

    if length(ind) == 1
        M = [M; {i} {x{ind}}];
    else
        a = fdot([{x{ind}}]);
        M = [M; {i} {a}]; % 1st column: reaction index - 2nd column: associated propensity
    end    
end

