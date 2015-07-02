%%% In this version, we take into account the time scale separation
%%% constraint without lattence.
%%% We take into account the constraints of the mass invariants.
%%% The updating policy is the full asynchronous policy (i.e. we take into account all types of updatings).
%%% We implemented an approximation of the upwards crossing interval
%%% constraint

%%% INPUTS
%%% v: vector of reactions
%%% ci: vector of initial conditions
%%% MMI_sign: matrix of mass invariants
%%% kin: vector of reaction rates
%%% nb_val: number of sampling intervals

function res = main_v3(ci,MMI_sign,v,kin,nb_val,model)
clc
%%Computation of the mass invariant
MI = mass_invar(ci,MMI_sign);

%Generation of the state space
nb_var = length(ci);
n = nb_var;
A = repmat(0:nb_val-1,1,nb_var);
B = unique(nchoosek(A,nb_var),'rows');
B = num2cell(B);
b = size(B);
nb_state = b(1);

%--------------------------------------------------------------------------
%Computation of the variables which are a product of a rule and a reactant
%of another (or the same) rule. This will be used in the computation of the
%interval crossing constraint

%Determination of the reactants and product indeces
reac_ind = []; prod_ind = [];
for i = 1: length(v)
    v_i = v{i};
    for j = 1:length(v_i)
        if v_i(j) < 0
            reac_ind = [reac_ind abs(v_i(j))];
        elseif v_i(j) > 0
            prod_ind = [prod_ind v_i(j)];        
        end   
    end
end
reac_ind = unique(reac_ind);
prod_ind = unique(prod_ind);
%intersec_reac_prod = intersect(reac_ind, prod_ind);

%--------------------------------------------------------------------------
%GENERATION OF THE ANNOTATION SPACE (will be used in the upwards crossing interval constraint)
An_sp = [];

%%Indeces of the variables which are not products (these variables cannot
%%be (-) annotated)
a = 1:nb_var;
a(prod_ind) = [];
ind_not_prod = a;

%%Generation of the annotation space
A = repmat(-1:0,1,nb_var);
B_annot = nchoosek(A,nb_var);
B_annot(:,ind_not_prod) = 0;
B_annot = unique(B_annot,'rows');

%--------------------------------------------------------------------------
% For each product, determination of the reactions {r_i} which produce the
% product (will be used in the upwards crossing interval constraint)

react_for_a_prod = [];
for i = 1:length(prod_ind)  %loop on the indeces of the products
    a = [];
    for j = 1: length(v)    %loop on the indeces of the reactions
        b = v{j};
        c = unique(b(b>0));
        for k = 1: length(c)%loop on the indeces of the products of the reaction of index j
            if prod_ind(i) - c(k) == 0
                a = [a, j];  %list of the indeces of the reactions which produce the product of index prod_ind(i)
                break
            end
        end    
    end
    d = [{prod_ind(i)} {a}];
    react_for_a_prod = [react_for_a_prod; d]; % 1st column of react_for_a_prod: indeces of the products
                                              % 2nd column of react_for_a_prod: indeces of the reactions which
                                              % produce the corresponding product
end
B;
%error('stop')
%--------------------------------------------------------------------------
FP = []; P = []; ns = 0; S = []; S_annot = [];
for i = 1:nb_state
    %----------------------------------------------------------------------
    %Computation of the states and its successors, which respect the mass
    %invariant constraint
    a = B(i,:); %state
    %a =[{0} {4} {0} {2} {2} {1}];
    %a =[{0} {1} {2}];
    b = invar_check (MI, MMI_sign, a);
    if b == 0 % b is not in the invariant
    %elseif (cell2mat(a(6)) == 1 || cell2mat(a(6)) == 2)  %state filtering prozone case ci242000
    %elseif (cell2mat(a(3)) == 7 || (cell2mat(a(1)) >= 6 && cell2mat(a(3)) >= 1) )  %state filtering rescaling case ci700
    elseif b == 1        
        S = concat_row_even_if_empty(S,a);
        ns = ns + 1;
        %list of states which are in the mass invariant
        %Computation of the discrete propensity for each reaction
        b = next_states(a,v,nb_val,kin,ci); % 1st column: next states - 2nd column: associated propensities
        size_b = size(b);
        d = []; % matrix of next states
        if length(b) == 0
            if length(FP)==0
                FP = [{a} {length(S(:,1))}];
            else
                FP = [FP; {a} {length(S(:,1))}]; %a has no successor so it is a fixed point
            end; 
        else
            c = cell2mat(b(:,1)); %extraction of the state from b
            c = num2cell(c);
            for j = 1:size_b(1)
                cc = invar_check(MI, MMI_sign, c(j,:));
                if cc == 0 %b is not in the invariant
                else
                    d = concat_row_even_if_empty(d,b(j,:)); % successors of the state a allowed by the mass invariant
                end
            end
        end
        
        %------------------------------------------------------------------
        %Computation of the successors allowed by the upwards crossing interval
        %constraint of all the combinations of annotated states        
        d_uci = [];
        for j = 1: length(B_annot(:,1)) %loop on the annotation space  
            list = [];
            list = [{[cell2mat(a) B_annot(j,:)]} {-1} {-1}];  %state a and its annotation
            S_annot = concat_row_even_if_empty(S_annot,{[cell2mat(a) B_annot(j,:)]});
            if length(d) == 0 %there is no successor of state a                
            else
                succ_a = cell2mat(d(:,1)); %extraction of the successors of state a
                for k = 1:length(succ_a(:,1))  %loop on the successors of state a
                    res = up_cross_int_const(a, B_annot(j,:), succ_a(k,:), react_for_a_prod,v);
                    if res.a == 1   %the transition is allowed by the upwards crossing interval constraint                    
                        b = [succ_a(k,:) res.b]; %annotation of the authorized successor
                        list = [list; {b} d(k,2) d(k,3)]; 
                    end    
                end
            end
            cell2mat(list(:,1));
            list = [list; {-1} {-1} {-1}];
            d_uci = [d_uci; list];
            cell2mat(list(1:end-1,1));
        end
        %------------------------------------------------------------------
        
         P = concat_row_even_if_empty(P,d_uci);
%        %P = [P ; d; {-1} {-1}];
         size_d = size(d);
    end   
end

%--------------------------------------------------------------------------


P(:,end) = []; %suppression of the last column corresponding to reaction associated to a transition
P;

S = cell2mat(S);
S_annot = cell2mat(S_annot);

%res.S = S; % list of the admissible states (respecting the invariant mass constraints)
res.S_annot = S_annot; % list of the admissible states (respecting the invariant mass constraints)
res.P = P; % list of the admissible states and their successors
res.FP = FP;
res.ns = ns;

ci = cell2mat(ci);

%backend(ci,S,P,FP,model);
backend(ci,S_annot,P,model);

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

for i = 1: length(v) % loop on the reactions
    P = [];
    
    reaction = v{i};
    y = cell2mat(x);
    
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
               
               %P, [{c} {M{i,2}}]
               %P = concat_row_even_if_empty(P,[{c} M{i,2}]);
               P = concat_row_even_if_empty(P,[{c} {M{i,2}} {v{i}}]);               

               cell2mat(P);               
            end
        end   
    end
    
    avant_supp_inter_cross = cell2mat(P);

    if length(P) == 0 
    else    
        N = concat_row_even_if_empty(N,P); % on stocke les transitions permises dans N
    end

end

%
% % if cell2mat(x) == [1 0 2]
% %     error('stop')
% % end
%     
% %if length(N) == 0 
% %		    avant_supp_prior = [];
% %  else 		   
% %avant_supp_prior = cell2mat(N);
% %end

%-------------------------------------------------------------------------
%Suppression of the successor states of lower priorities
list = [];
if ~isempty(N)
    NN = N(:,2); %extraction of the column which stores the min and max of the propensities
    NN = cell2mat(NN);
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

%if length(N)=0 
%    apres_supp_prior = [];
%  else
apres_supp_prior = cell2mat(N(:,1:2));
%end
%-------------------------------------------------------------------------
%Setting of the propensities to 1
if length(N)~=0
    for i = 1:length(N(:,1))
        N{i,2} = 1;
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
	    M = concat_row_even_if_empty(M,[{i} {a}]);
        else
            aa = [sum(x_ind) + kin(i), sum(x_ind) + 1 + kin(i)];
            a = [min(aa),max(aa)];
	    M = concat_row_even_if_empty(M,[{i} {a}]);
        end
    elseif (length(ind) == 2 && length(unique(ind)) == 1) % bimolecular reaction of the form 2A ->*
        x_ind = [x{ind}];
        if prod(x_ind) ~= 0
            aa = [2*sum(x{ind})-1+kin(i), 2*sum(x{ind})+kin(i), 2*sum(x{ind})+kin(i)+1];
            a = [min(aa),max(aa)];
	    M = concat_row_even_if_empty(M,[{i} {a}]);
        else
            aa = [0 : kin(i) + 1];
            a = [min(aa),max(aa)];
	    M = concat_row_even_if_empty(M,[{i} {a}]);
        end 
    end

end

function res = max_resource_prod(state, react_for_a_prod, v)
%Function which computes the resource available for the production of each
%product at a given state, and the indeces of the reactants which equal the
%min of the reactants of the reactions which produce each product
res = [];
for i = 1:length(react_for_a_prod(:,1))
    aa = cell2mat(react_for_a_prod(i,2)); %indeces of the reactions which produce the product of index i
    b = []; c = [];
    for j = 1:length(aa)
        vv = v{aa(j)};
        vvv = unique(abs(vv(vv<0)));         %indeces of the reactants of the reaction vv{aa(j)}
        b = [b, min(cell2mat(state(vvv)))];  %vector of the min of the reactants of the reactions producing the product of index i
        
        for k = 1: length(vvv)
            d = cell2mat(state(vvv));
            if d(k) == min(d)
                c = [c, vvv(k)]; %indeces of the reactants which equal the min of the reactants of the reactions producing the product of index i
            end
        end
        
    end
    c = unique(c);
    e = [react_for_a_prod(i,1), {max(b)}, {c}];
    res = [res; e]; % 1st column: indeces of the products; 2nd column: resources available for the corresponding product; 
                    % 3rd column: indeces of the reactants which equal
                    % the min of the reactants of the reactions producing
                    % the corresponding product
end

function res = up_cross_int_const(state, annot, succ_a, react_for_a_prod,v)

% Elimination of the transitions forbidden by the upwards crossing interval
% constraint and annotation of the allowed transitions
%%% Two conditions have to be fullfilled to suppress a transition:
%%% 1) none of the reactants  which value is min(reactants) for each r_i is a
%%% product of a reaction
%%% 2) there exists an annotated product updated by the transition for which its value
%%% is > max(min(reactants))
%%% A non annotated variable is annotated after a transition if it is positively updated after the transition
%%% An annotated variable keeps annotated after a transition if it is not
%%% negatively updated after the transition

annot_succ = zeros(1,length(state));
bool = 1;

a = find(annot); %extraction of the indeces of the variables which are annotated
products = cell2mat(react_for_a_prod(:,1)); %indeces of the products
[prod_annot,i_prod,i_a] = intersect(products,a); %indeces of the products which are annotated

delta = succ_a - cell2mat(state);
b1 = []; b2 = [];
for i = 1:length(delta)
    if delta(i) > 0
        b1 = [b1, i]; %extraction of the indeces of the positively updated variables in the transition state -> succ_a
    end
end

if length(a) == 0
else
    for i = 1: length(a)
        if delta(a(i)) == 0
            b2 = [b2, a(i)]; %extraction of the indeces of the annotated variables which are not updated in the transition state -> succ_a
        end
    end
end

annot_succ(b1) = -1; annot_succ(b2) = -1;  %annotation of the successor succ_a
[c, i_pa, i_b] = intersect(prod_annot, b1); %indeces of the annotated products which are updated in the transition state -> succ_a

d = max_resource_prod(state, react_for_a_prod, v);
e = d(i_prod(i_pa),:); %extraction of the rows corresponding to the annotated products which are updated in the transition state -> succ_a

for j = 1:length(c)
    if cell2mat(state(cell2mat(e(j,1)))) > cell2mat(e(j,2))  %1st condition to suppress a transition            
        if length(intersect(cell2mat(e(j,3)),products)) == 0 %2nd condition to suppress a transition
            bool = 0;
            break
        end
    end
end
res.a = bool;
res.b = annot_succ;

