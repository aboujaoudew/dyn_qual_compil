%%% In this version, we take into account the time scale separation
%%% constraint without lattence.
%%% We take into account the constraints of the mass invariants.
%%% The updating policy is the full asynchronous policy (i.e. we take into account all types of updatings).
%%% We updated the approximation of the upwards crossing interval
%%% constraint (which is still not correct)
%%% We made a function for the constraint on time scale separation

%%% INPUTS
%%% v: vector of reactions
%%% ci: vector of initial conditions
%%% MMI_sign: matrix of mass invariants
%%% kin: vector of reaction rates
%%% nb_val: number of sampling intervals

function res = main_v4(ci,MMI_sign,vect_reac,kin,nb_val,model)
clc

%Generation of the state space
nb_var = length(ci);
n = nb_var;
A = repmat(0:nb_val-1,1,nb_var);
B = unique(nchoosek(A,nb_var),'rows');
state_space = num2cell(B);
b = size(state_space);
nb_state = b(1);

%--------------------------------------------------------------------------
%Computation of the variables which are a product of a rule and a reactant
%of another (or the same) rule. This will be used in the computation of the
%interval crossing constraint

%Determination of the reactants and product indeces
reac_ind = []; prod_ind = [];
for i = 1: length(vect_reac)
    v_i = vect_reac{i};
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

FP = []; P = []; ns = 0; S = []; S_annot = [];
%--------------------------------------------------------------------------
%Computation of the set of transitions
trans_set = trans_set(state_space, vect_reac, nb_state, nb_val);

%--------------------------------------------------------------------------
%Computation of the state space and the set of transitions refined with mass invariant
res = trans_set_mi(state_space, trans_set, MMI_sign, ci);
state_space = res.ss;
trans_set = res.ts;

%--------------------------------------------------------------------------
%Computation of the set of transitions refined with scale separation
trans_set = trans_set_scalesep(trans_set, vect_reac, kin);

%%--------------------------------------------------------------------------
%Computation of the set of transitions refined with the upwards crossing
%interval constraint
res = trans_set_upcrossint(state_space, trans_set, prod_ind, vect_reac);
trans_set = res.ts;
state_space = res.ss;

%%%Concatenation of the states of a transition and its annotation
a = [cell2mat(trans_set(:,1)) cell2mat(trans_set(:,4))];
b = mat2cell(a, [zeros(1,length(a(:,1))) + 1], [length(a(1,:))]);

a = [cell2mat(trans_set(:,2)) cell2mat(trans_set(:,5))];
c = mat2cell(a, [zeros(1,length(a(:,1))) + 1], [length(a(1,:))]);

trans_set = [b, c, trans_set(:,3)];

%--------------------------------------------------------------------------

ci = cell2mat(ci);

%backend(ci,S,P,FP,model);
%backend(ci,S_annot,P,model);
backend_v2(ci, state_space, trans_set, model);

function res = trans_set(state_space, v, nb_state, nb_val)
%Computation of the set of transitions
res = [];
for i = 1:nb_state
    %----------------------------------------------------------------------
    %Computation of the states and its successors, which respect the mass
    %invariant constraint

    state = state_space(i,:); %state
    trans = next_states(state,v,nb_val); % 1st column: state; 2nd column: successor; 3rd column: reaction    
    res = [res; trans];
end 

function res = next_states(state,v,nb_val)
%Computation of the next states of state y: we consider a full asynchronous
%updating policy between the variables of a reaction. 

res = [];
etat = cell2mat(state);
%--------------------------------------------------------------------------

for i = 1: length(v) % loop on the reactions
    P = [];
    
    reaction = v{i};
    y = cell2mat(state);
    
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
                    if y(vvv_prod_ind(l)) == nb_val-1 % constraint on the maximum value of the variables (=0)
                        flag = 1;
                        break;
                    end
                end
            end
            
            if flag == 1
            else
               succ = y;
               succ(vvv_reac_ind) = y(vvv_reac_ind) - 1;
               succ(vvv_prod_ind) = y(vvv_prod_ind) + 1;
               
               %P, [{c} {M{i,2}}]
               %P = concat_row_even_if_empty(P,[{c} M{i,2}]);
               res = concat_row_even_if_empty(res,[{y} {succ} {i}]);   %1st column: states; 2nd column: successors; 3rd column: reaction indeces
               cell2mat(P);            
            end
        end   
    end
    
    avant_supp_inter_cross = cell2mat(P);

%     if length(P) == 0 
%     else    
%         res = concat_row_even_if_empty(N,P); % on stocke les transitions permises dans N
%     end

end

function res = trans_set_mi(state_space, trans_set, MMI_sign, ci)
%Computation of the set of transitions refined with the mass invariant
%(MMI_sign, ci)

%%Computation of the mass invariant constant
MI = mass_invar(ci,MMI_sign);

list1 = []; list2 = [];
for i = 1:length(state_space(:,1))    
    a = invar_check (MI, MMI_sign, state_space(i,:));
    if a == 1           % state_space(i,1) is in the invariant
        list1 = [list1, i];
    end
end

for i = 1:length(trans_set(:,1))
    a = invar_check (MI, MMI_sign, num2cell(cell2mat(trans_set(i,1))));
    b = invar_check (MI, MMI_sign, num2cell(cell2mat(trans_set(i,2)))); 
    if a == 1 && b == 1 % trans_set(i,1) and trans_set(i,2) are in the invariant
        list2 = [list2, i];
    end        
end

res.ss = state_space(list1,:);
res.ts = trans_set(list2,:);

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

function res = trans_set_scalesep(trans_set, v, kin)
%Computation of the set of transitions refined with scale separation

%%Computation of the propensities of the reactions
concat_prop = [];
for i = 1:length(trans_set(:,1))
    prop = prop_comp(trans_set(i,1), v{cell2mat(trans_set(i,3))}, kin(cell2mat(trans_set(i,3))));    
    concat_prop = [concat_prop; prop];
end
trans_set_prop = [trans_set, concat_prop];

%%Suppression of the transitions of lower priorities

%Extraction of the transitions of common initial state
[a, ia, ib] = unique(cell2mat(trans_set(:,1)),'rows');

res = []; concat_stock = [];
for i = 1:max(ib)

    trans_set_state = [];
    b = find(ib==i);
    trans_set_state = trans_set_prop(b,:); % extraction of a block of trans_set of common initial state
        
    stock = [];
    prop_state = cell2mat(trans_set_state(:,4));
    size_ps = size(prop_state);

    for j = 1:size_ps(1)
        for k = 1:size_ps(1)
            if max(prop_state(j,:)) < min(prop_state(k,:))
                stock = [stock, b(j)]; %we store the transition of lower priorities
                break; 
            else
            end
        end
    end

    concat_stock = [concat_stock, stock]; 
    
end

res = trans_set;
res(concat_stock,:) = []; %suppression of the transition of lower priorities

function res = prop_comp(x, v, kin)
% Computation of the propensities of the reaction v of kinetic constant kin at state x

ind = [];
a = v;
for j = 1: length(a)
    if a(j) < 0
        ind = [ind abs(a(j))]; %list of the reactant indeces of reaction i
    end
end

if length(ind) == 1 %monomolecular reaction
    y = cell2mat(x);
    if y(ind) ~= 0
        res = {[y(ind) + kin, y(ind) + kin]};
    else
        aa = [0 : kin];
        res = {[min(aa), max(aa)]};
    end
 
elseif (length(ind) == 2 && length(unique(ind)) == 2) % bimolecular reaction of the form A + B ->*
    y = cell2mat(x);
    if prod(y(ind)) == 0
        aa = [0 : sum(y(ind)) + 1 + kin];
        res = {[min(aa),max(aa)]};
    else
        aa = [sum(y(ind)) + kin, sum(y(ind)) + 1 + kin];
        res = {[min(aa),max(aa)]};
    end
    %cell2mat(res)
elseif (length(ind) == 2 && length(unique(ind)) == 1) % bimolecular reaction of the form 2A ->*
    y = cell2mat(x);
    if prod(y(ind)) ~= 0
        aa = [2*y(ind) - 1 + kin, 2*y(ind) + kin, 2*y(ind) + kin + 1];
        res = {[min(aa),max(aa)]};
    else
        aa = [0 : kin + 1];
        res = {[min(aa),max(aa)]};
    end 
end

function res = trans_set_upcrossint(state_space, trans_set, prod_ind, vect_reac)
%--------------------------------------------------------------------------
%Computation of the successors allowed by the upwards crossing interval
%constraint of all the combinations of annotated states        

%--------------------------------------------------------------------------
%GENERATION OF THE ANNOTATION SPACE

%%Indeces of the variables which are not products (these variables cannot
%%be (-) annotated)
a = 1:length(state_space(1,:));
a(prod_ind) = [];
ind_not_prod = a;

%%Generation of the annotation space
A = repmat(-1:0,1,length(state_space(1,:)));
AA = nchoosek(A,length(state_space(1,:)));
AA(:,ind_not_prod) = 0;
AA = unique(AA,'rows');
annot_space = mat2cell(AA, [zeros(1,length(AA(:,1))) + 1], [length(AA(1,:))]);
zero_annot_space = mat2cell(zeros(length(AA(:,1)),length(AA(1,:))), [zeros(1,length(AA(:,1))) + 1], [length(AA(1,:))]);
%--------------------------------------------------------------------------
%GENERATION OF THE ANNOTATED STATE SPACE
state_space_annot = [];
for i = 1:length(state_space(:,1))
    a = repmat(state_space(i,:),[length(annot_space(:,1)),1]);
    state_space_annot = [state_space_annot; a num2cell(cell2mat(annot_space))];
end

%--------------------------------------------------------------------------
%Determination of the reactions which produce a product

react_for_a_prod = [];
for i = 1:length(prod_ind)  %loop on the indeces of the products
    a = [];
    for j = 1: length(vect_reac)    %loop on the indeces of the reactions
        b = vect_reac{j};
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
%--------------------------------------------------------------------------

trans_set_preannot = [];
for i = 1:length(trans_set(:,1))
    a = repmat(trans_set(i,:),length(annot_space(:,1)),1);
    b = [a, annot_space, zero_annot_space];
    trans_set_preannot = [trans_set_preannot; b]; %4th column: annotation of the state; last column: annotation of the successor set to 0
end

cell2mat(trans_set_preannot);

list = [];
for i = 1: length(trans_set_preannot(:,1))    %loop on the transition of trans_set_preannot
   
    res = up_cross_int_const(trans_set_preannot(i,1), trans_set_preannot(i,4), trans_set_preannot(i,2), react_for_a_prod, vect_reac);
    
    if res.a == 1   %the transition is allowed by the upwards crossing interval constraint                    
        trans_set_preannot(i,5) = {res.b};           %annotation of the successor of the authorized transitions
        list = [list, i];                            %storage of the indices of the authorized transitions
    end
end    

res.ts = trans_set_preannot(list,:);
res.ss = state_space_annot;

function res = up_cross_int_const(state, state_annot, succ_state, react_for_a_prod, vect_reac)

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

annot_succ = zeros(1,length(cell2mat(state)));
bool = 1;

a = find(cell2mat(state_annot)); %extraction of the indeces of the variables which are annotated
products = cell2mat(react_for_a_prod(:,1)); %indeces of the products
[prod_annot,i_prod,i_a] = intersect(products,a); %indeces of the products which are annotated

delta = cell2mat(succ_state) - cell2mat(state);
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

[c, i_pa, i_b] = intersect(prod_annot, b1); %indeces of the annotated products which are updated in the transition: state -> succ_a

d = max_resource_prod(state, react_for_a_prod, vect_reac);  %computation of the max of resources available for the production of each product
e = d(i_prod(i_pa),:);  %extraction of the rows corresponding to the annotated products which are updated in the transition: state -> succ_a

%--------------------------------------------------------------------------
%Annotation of the successor succ_state
annot_succ(b1) = -1; annot_succ(b2) = -1;  

%--------------------------------------------------------------------------
%Suppression of the forbidden transitions
state_num = cell2mat(state);
for j = 1:length(c)
    if state_num(cell2mat(e(j,1))) > cell2mat(e(j,2))  %1st condition to suppress a transition            
        if length(intersect(cell2mat(e(j,3)),products)) == 0 %2nd condition to suppress a transition
            bool = 0;
            break
        end
    end
end

res.a = bool;
res.b = annot_succ;

function res = max_resource_prod(state, react_for_a_prod, vect_reac)
%Function which computes the resource available for the production of each
%product at a given state, and the indeces of the reactants which equal the
%min of the reactants of the reactions which produce each product
cell2mat(state)
res = [];
for i = 1:length(react_for_a_prod(:,1))
    react_for_a_prod(i,1)
    aa = cell2mat(react_for_a_prod(i,2)); %indeces of the reactions which produce the product of index: react_for_a_prod(i,1)
    b = []; c = [];
    for j = 1:length(aa)
        vv = vect_reac{aa(j)};
        vvv = unique(abs(vv(vv<0)));         %indeces of the reactants of the reaction vv{aa(j)}        
        num_state = cell2mat(state);
        b = [b, min(num_state(vvv))];  %vector of the min of the reactants of the reactions producing the product of index: react_for_a_prod(i,1)
        
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



