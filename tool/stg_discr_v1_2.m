%%% In this version, we implemented the updating scheme by scheduling the
%%% reactions according to their time scale (priority classes). 
%%% We take into account the constraints from the mass invariants
%%% The updating policy is an extended asynchronous policy.
%%% We introduce the number of values as a parameter (nb_val), and kinetic
%%% constants for reactions (kin)
%%% We use the additive law defined in fplus_R 

%%%INPUTS
%%% v: vector of reactions
%%% ci: vector of initial conditions
%%% MMI_sign: matrix of mass invariants

function res = stg_discr_v1_2(ci,MMI_sign,v,kin)

%%Computation of the mass invariant
%MI = mass_invar(ci,MMI_sign)

% DEFINITION OF THE MASS INVARIANT
% MI = [{1}; {1}; {1}]; % for the prozone
MI = [{5}]; % for the rescaling

%Generation of the state space
nb_val = 5;
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
    a = B(i,:);
    % a = [{2} {2} {0}];
    % a = [{2} {1} {2} {0} {0} {0}];
    % a = [{1} {1} {0} {0} {1}];

    b = invar_check (MI, MMI_sign, a);
    if b == 0 % b is not in the invariant
    elseif b == 1
        S = [S; a];
        ns = ns + 1;
        P = [P ; {cell2mat(a)} {-1}]; %list of states which are in the mass invariant
        %Computation of the discrete propensity for each reaction
        b = next_states(a,v,nb_val,kin); % 1st column: next states - 2nd column: associated propensities
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

%save_report(ci,res);
dot_file(ci,S,P,FP);

function N = invar_check(MI, MMI_sign, x) 
% N = 0: x is not in the invariant; N = 1: x is in the invariant
N = 1;
% a = mass_invar(x,MMI_sign);
a = mass_invar_R(x, MMI_sign);
for i = 1:length(MI)
    b = intersect(a{i},MI{i});
    if length(b) ~= 0 % at least one element of a{i} is contained in MI{i}
    else
        N = 0;
        break        
    end
end

function N = next_states(x,v,nb_val,kin)
    
%%Computation of the propensities of the reactions
M = prop_comp(x,v,kin); % 1st column: reaction index - 2nd column: associated propensity

%Computation of the next states of state y: we consider an asynchronous
%policy between the variables of a reaction. The transitions are labeled with the
%propensity of the reaction
N = [];
etat = cell2mat(x)

for i = 1: length(v) % loop on the reactions
    %Computation of the mass constraints induced by the reaction v{i}
    mcr = mass_constr_reac(x,v{i},nb_val);    
    reaction = v{i};
    MCR_prod = mcr.prod; % MCR_prod(:,1): index of the product variables - MCR_prod(:,2): associated born
    MCR_reac = mcr.reac; % MCR_reac(:,1): index of the reactant variables - MCR_reac(:,2): associated born

    %Reactants and product indeces
    reac_ind = MCR_reac(:,1);
    prod_ind = MCR_prod(:,1);
    
    y = cell2mat(x);

    %%%====================================================================
    %Full asynchronoues updating
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
                flag = 0;
                c = y;           
                c(vvv_reac_ind) = y(vvv_reac_ind) - 1;
                c(vvv_prod_ind) = y(vvv_prod_ind) + 1;
                c;
                if length(vvv_reac_ind) == 0
                else
                    for l = 1:length(vvv_reac_ind)
                        a = MCR_reac(:,1);
                        b = find(a==vvv_reac_ind(l));
                        % if c(vvv_reac_ind(l)) < MCR_reac(b,2) % constraint on the maximum consumption of the reactants 
                        %     flag = 1;
                        %     break;
                        % else
                        % end
                    end
                end
                if length(vvv_prod_ind) == 0
                else            
                    for l = 1:length(vvv_prod_ind)
                        a = MCR_prod(:,1);
                        b = find(a==vvv_prod_ind(l));
                        % if c(vvv_prod_ind(l)) > MCR_prod(b,2) % constraint on the maximum production of the products
                        %     flag = 1;
                        %     break;
                        % else    
                        % end
                    end
                end
                if flag == 1
                else
                    N = [N; {c} M{i,2}];
                    cell2mat(N);
                end    
            end
        end
    end
    %%%====================================================================
%     %Extended asynchronous updating between the variables of a reaction
%     %%% Asynchronous updating between the reactants of a reaction
%     for j = 1:length(reac_ind)
%         c = y;
%         if y(reac_ind(j)) ~= 0 % constraint on the minimum value of the variables (=0)
%             c(reac_ind(j)) = y(reac_ind(j)) - 1;
%             if c(reac_ind(j)) < MCR_reac(j,2) % constraint on the maximum consumption of the reactants 
%             else
%                 N = [N; {c} M{i,2}];
%             end
%         end
%     end
%     
%     %%% Synchronous updating of the reactants of a reaction when the
%     %%% values of the reactants equal the minimum value of all the
%     %%% reactants
%     b = [];
%     a = min(y(reac_ind));
%     if a ~= 0
%         % Determination of the reactants which equal the min of the reactant
%         % values
%         for j = 1:length(reac_ind)
%             if y(reac_ind(j)) == a
%                 b = [b reac_ind(j)];
%             end
%         end
%         % Determination of the next states
%         if length(b) >= 2
%             for j = 2:length(b)
%                 d = nchoosek(b,j);
%                 c = y;
%                 c(d) = y(d) - 1;
%                 for k = 1:length(c(:,1))
%                     N = [N; {c(k,:)} M{i,2}];
%                 end
%             end         
%         end
%     end
%     
%     %%% Asynchronous updating between the products of a reaction
%     for j = 1:length(prod_ind)
%         c = y;
%         if y(prod_ind(j)) ~= nb_val % constraint on the maximum value of the variables (=2)
%             c(prod_ind(j)) = y(prod_ind(j)) + 1;
%             if c(prod_ind(j)) > MCR_prod(j,2) % constraint on the maximum production of the products
%             else
%                 N = [N; {c} M{i,2}];
%             end
%         end
%     end
    %%%====================================================================

end

avant_supp = cell2mat(N)

%Suppression of the successor states of lower priorities
list = [];
if length(N) > 0
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
apres_supp = cell2mat(N)

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
        if prod(x{ind}) ~= 0
            aa = [0 : sum(x{ind}) + 1 + kin(i)];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        else
            aa = [sum(x{ind}) + kin(i), sum(x{ind}) + 1 + kin(i)];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        end
    elseif (length(ind) == 2 && length(unique(ind)) == 1) % bimolecular reaction of the form 2A ->*
        x_ind = x{ind};
        if prod(x_ind) ~= 0
            aa = [2*sum(x{ind})-1+kin(i), 2*sum(x{ind})+kin(i), 2*sum(x{ind})+kin(i)+1];
            a = [min(aa),max(aa)];
            M = [M; {i} {a}];
        else
            aa = [0 : kin(i) + 1]
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

 