%%% In this version, we implemented the updating scheme by scheduling the
%%% reactions according to their time scale (priority classes). 
%%% We take into account the constraints from the mass invariants and from 
%%% the borns on the derivative (i.e. limiting reactant constraint)
%%% The updating policy is an extended asynchronous policy.
%%% We introduce the number of values as a parameter (nb_val), and kinetic
%%% constants for reactions (kin)
%%% WE SUPPRESS THE BORN TO THE DERIVATIVE CONSTRAINT AND IMPLEMENTED
%%% INSTEAD THE INTERVAL CROSSING CONSTRAINT (TO UPDATE!)

%%%INPUTS
%%% v: vector of reactions
%%% ci: vector of initial conditions
%%% MMI_sign: matrix of mass invariants

function res = main_v1(ci,MMI_sign,v,kin,nb_val)

%%Computation of the mass invariant
MI = mass_invar(ci,MMI_sign);

%Generation of the state space
nb_var = length(ci);
n = nb_var;
A = repmat(0:nb_val,1,nb_var);
B = unique(nchoosek(A,nb_var),'rows');
B = num2cell(B);
b = size(B);
nb_state = b(1);

%Computation of the next states
FP = []; P = []; ns = 0; S = [];
for i = 1:nb_state
    a = B(i,:); %state
    b = invar_check (MI, MMI_sign, a);

    if b == 0 % b is not in the invariant
    %elseif (cell2mat(a(6)) == 1 || cell2mat(a(6)) == 2)  %state filtering prozone case ci242000
    elseif (cell2mat(a(3)) == 7 || (cell2mat(a(1)) >= 6 && cell2mat(a(3)) >= 1) )  %state filtering rescaling case ci700
    elseif b == 1
        S = [S; a];
        ns = ns + 1;
        P = [P ; {cell2mat(a)} {-1}]; %list of states which are in the mass invariant
        %Computation of the discrete propensity for each reaction
        b = next_states(a,v,nb_val,kin,ci); % 1st column: next states - 2nd column: associated propensities
        size_b = size(b);
        d = []; % matrix of next states
        if length(b) == 0
            FP = [FP; {a} {length(S(:,1))}]; %a has no successor so it is a fixed point
        else
            c = cell2mat(b(:,1)); %extraction of the state from b
            c = num2cell(c);
            for j = 1:size_b(1)
                cc = invar_check(MI, MMI_sign, c(j,:));
                if cc == 0 %b is not in the invariant
                else
                    d = [d ; b(j,:)];
                end
            end
        end
        P = [P ; d; {-1} {-1}];
        size_d = size(d);
    end
end

S = cell2mat(S);

res.S = S; % list of the admissible states (respecting the invariant mass constraints)
res.P = P; % list of the admissible states and their successors
res.FP = FP;
res.ns = ns;


ci = cell2mat(ci);

dot_file(ci,S,P,FP);

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

function N = next_states(x,v,nb_val,kin,ci)
    
%%Computation of the propensities of the reactions
M = prop_comp(x,v,kin); % 1st column: reaction index - 2nd column: associated propensity

%Computation of the next states of state y: we consider an asynchronous
%policy between the variables of a reaction. The transitions are labeled with the
%propensity of the reaction
N = [];

etat = cell2mat(x)

%--------------------------------------------------------------------------
%Computation of the variables which are a product of a rule and a reactant
%of another (or the same) rule. This will be used in the computation of the
%interval crossing constraint
reac_ind = []; prod_ind = [];
for i = 1: length(v)
    mcr = mass_constr_reac(x,v{i},nb_val);
    MCR_prod = mcr.prod; % MCR_prod(:,1): index of the product variables - MCR_prod(:,2): associated born
    MCR_reac = mcr.reac; % MCR_reac(:,1): index of the reactant variables - MCR_reac(:,2): associated born
    
    %Reactants and product indeces of reaction v(i)
    a = MCR_reac(:,1);
    b = MCR_prod(:,1);
    
    reac_ind = [reac_ind a];
    prod_ind = [prod_ind b];
end
reac_ind = unique(reac_ind);
prod_ind = unique(prod_ind);
intersec_reac_prod = intersect(reac_ind, prod_ind);
%--------------------------------------------------------------------------

for i = 1: length(v) % loop on the reactions
    P = [];
    
    reaction = v{i}
    y = cell2mat(x);

    %%%====================================================================
    %Full asynchronous updating
    a = v{i};
    v_uniq = unique(a','rows');
    v_uniq = v_uniq';
    for j = 1:length(v_uniq)
        vv = nchoosek(v_uniq,j);
        for k = 1: length(vv(:,1))
            k;
            flag = 0;
            vvv = vv(k,:);
            vvv_reac_ind = abs(vvv(vvv<0));
            vvv_prod_ind = abs(vvv(vvv>0));
            if length(vvv_reac_ind)~=0
                for l = 1:length(vvv_reac_ind)
                    if y(vvv_reac_ind(l)) == 0 % constraint on the minimum value of the variables (=0)
                        flag = 1;
                        break;
                    end
                end
            end            
            if length(vvv_prod_ind)~=0
                for l = 1:length(vvv_prod_ind)
                    if y(vvv_prod_ind(l)) == nb_val % constraint on the maximum value of the variables (=0)
                        flag = 1;
                        break;
                    end
                end
            end
            
            if flag == 1
            else
               c = y;
               c(vvv_reac_ind) = y(vvv_reac_ind) - 1;
               c(vvv_prod_ind) = y(vvv_prod_ind) + 1;
               P = [P; {c} M{i,2}];
               cell2mat(P);
            end    
        end   
    end
%-------------------------------------------------------------------------
%Suppression of the transitions which are forbidden by the interval
%crossing constraint
    
    avant_supp_inter_cross = cell2mat(P)
    
    if length(P) > 0
        list_a = [];
        PP = cell2mat(P);
        PP = PP(:,1:end-2); %suppression of the last 2 columns of NN which correspond to the min and max of the propensities
        size_P = size(P);
        cii = cell2mat(ci);
        y = cell2mat(x);
        
        a = y - cii; 
        %Determination of the list of the indeces of the variables of the state y which have been positively updated compared to the initial condition
        for i = 1:length(a) 
            if a(i) > 0
                list_a = [list_a i]; 
            end              
        end
        list_a;
        b = abs(v_uniq(v_uniq < 0)); %list of the reactant indeces of reaction v_uniq

        list_c = [];
        for i = 1:size_P(1)
            list_b = [];
            c = PP(i,:) - y;
            %Determination of the list of the indeces of the variables of the
            %transitions which have been positively updated compared to the state y
            for j = 1:length(PP(i,:))
               if c(j) > 0
                    list_b = [list_b j]; %list of the indeces of the variables of the state x which have been positively updated compared to the initial condition
               end              
            end
            list_b = unique(list_b);
            
            %Intersection of list_a and list_b
            list_ab = intersect(list_a,list_b);
            
            %Determination of the transitions which violate the constraint
            %on the interval crossing
            for j = 1:length(list_ab)
                reac_min = [];
                y_b = y(b);
                for k = 1:length(y_b) % determination of the reactants indeces of min value
                    if y_b(k) == min(y(b))
                        reac_min = [reac_min, b(k)];
                    end
                end
                if y(list_ab(j)) > min(y(b)) % 1st condition to violate the constraint
                   if  length(intersect(reac_min, intersec_reac_prod)) == 0 % 2nd condition to violate the constraint
                       list_c = [list_c i];
                       break
                   end
                end
            end
        end
        list_c;
        P;
        P(list_c,:) = [] %suppression of the transitions which violate the constraint
                      %on the crossing intervals
    end
    
    apres_supp_inter_cross = cell2mat(P)
    
    N = [N; P]; % on stocke les transitions permises dans N
end

% if cell2mat(x) == [1 0 2]
%     error('stop')
% end
    
avant_supp_prior = cell2mat(N)

%-------------------------------------------------------------------------
%Suppression of the successor states of lower priorities
list = [];
if ~isempty(N)
    NN = cell2mat(N);
    size_N = size(N);
    for i = 1:size_N(1)
        for j = 1:size_N(1)
            if max(NN(i,end-1:end)) < min(NN(j,end-1:end))
                list = [list, i];
                break; % lower priorities successor states are suppressed
            else
            end
        end
    end
end
N(list,:) = [];

apres_supp_prior = cell2mat(N)

%-------------------------------------------------------------------------
if length(N)~=0
    for i = 1:length(N(:,1))
        N{i,end} = 1;
    end
end

% if y == [2 1 2 0 0 0]
%     error('stop')
% end

function M = prop_comp(x,v,kin)

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
    if length(ind) == 1 %monomolecular reaction
        if x{ind} ~= 0
            M = [M; {i} {[x{ind} + kin(i), x{ind} + kin(i)]}]; %for consistency when concatenating the propensities
        else
            aa = [0 : kin(i)];
            a = [min(aa), max(aa)];          
            M = [M; {i} {a}]; % 1st column: reaction index - 2nd column: associated propensity
        end
    elseif (length(ind) == 2 && length(unique(ind)) == 2) % bimolecular reaction of the form A + B ->*
        %ind(1)
        x_ind = [x{ind}];
        %x{ind}
        %sum(x_ind)
        %error('stop')
        if prod(x_ind) == 0
            aa = [0 : sum(x_ind) + 1 + kin(i)];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        else
            aa = [sum(x_ind) + kin(i), sum(x_ind) + 1 + kin(i)];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        end
    elseif (length(ind) == 2 && length(unique(ind)) == 1) % bimolecular reaction of the form 2A ->*
        x_ind = [x{ind}];
        if prod(x_ind) ~= 0
            aa = [2*sum(x{ind})-1+kin(i), 2*sum(x{ind})+kin(i), 2*sum(x{ind})+kin(i)+1];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        else
            aa = [0 : kin(i) + 1];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        end 
    end
    
% % Computation of the propensities of the reactions v(i) at state x
% M = [];
% for i = 1:length(v)
% ind = [];
% a = v{i}; % reaction i
%     for j = 1: length(a)
%         if a(j) < 0
%             ind = [ind abs(a(j))]; %list of the reactant indeces of reaction i
%         end
%     end
% 
%     if length(ind) == 1
%         M = [M; {i} {[x{ind} + kin(i), x{ind} + kin(i)]}]; %for consistency when concatenating the propensities
%     else
%         aa = fdot_k(x{ind},kin(i));
%         a = [min(aa),max(aa)];
%         M = [M; {i} {a}] % 1st column: reaction index - 2nd column: associated propensity
%     end    
end

 