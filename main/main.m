%%% In this version, we take into account the time scale separation
%%% constraint without lattence.
%%% We take into account the constraints of the mass invariants.
%%% The updating policy is the full asynchronous policy (i.e. we take into account all types of updatings).
%%% We implemented a correct approximation of the upwards crossing interval
%%% constraint
%%% We made a function for the constraint on time scale separation

%%% Compared to the previous version, we made a separated module for the 
%%% quotienting of the states

%%% INPUTS
%%% v: vector of reactions
%%% ci: vector of initial conditions
%%% MMI_sign: matrix of mass invariants
%%% kin: vector of reaction rates
%%% nb_val: number of sampling intervals

function res = main(ci,MMI_sign,vect_reac,kin,nb_val,model)

%Generation of the state space
nb_var = length(ci);
n = nb_var;
A = repmat(0:nb_val-1,1,nb_var);
B = unique(nchoosek(A,nb_var),'rows');
state_space = num2cell(B);
b = size(state_space);
nb_state = b(1);
disp('generation of the state space done')

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
%Computation of the set of abstract transitions
trans_set = trans_set(state_space, vect_reac, nb_state, nb_val);
disp('computation of the abstract transitions done')

%--------------------------------------------------------------------------
%Computation of the state space and the set of transitions refined with mass invariant
res = trans_set_mi(state_space, trans_set, MMI_sign, ci);
state_space = res.ss;
trans_set = res.ts;
disp('refinement with mass invariant done')

%--------------------------------------------------------------------------
%Computation of the set of transitions refined with scale separation
trans_set = trans_set_scalesep(trans_set, vect_reac, kin);
disp('refinement with scale separation done')

%%--------------------------------------------------------------------------
%Computation of the set of transitions refined with the upwards crossing
%interval constraint
res = trans_set_upcrossint(state_space, trans_set, prod_ind, vect_reac);
trans_set = res.ts;
state_space = res.ss;
disp('refinement with the upwards crossing interval constraint and quotienting of states done')

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
backend(ci, state_space, trans_set, model);

%==========================================================================
%COMPUTATION OF THE SET OF ABSTRACT TRANSITIONS
%==========================================================================

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
               %cell2mat(P);            
            end
        end   
    end
    
    %avant_supp_inter_cross = cell2mat(P);

%     if length(P) == 0 
%     else    
%         res = concat_row_even_if_empty(N,P); % on stocke les transitions permises dans N
%     end

end

%==========================================================================
%COMPUTATION OF THE SET OF ABSTRACT TRANSITIONS REFINED WITH MASS
%INVARIANTS
%==========================================================================

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

%==========================================================================
%COMPUTATION OF THE SET OF ABSTRACT TRANSITIONS REFINED WITH TIME SCALE
%SEPARATION
%==========================================================================

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

%==========================================================================
%COMPUTATION OF THE SET OF ABSTRACT TRANSITIONS REFINED WITH THE UPWARDS
%CROSSING INTERVAL CONSTRAINT
%==========================================================================

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
for i = 1:length(prod_ind)   %loop on the indeces of the products
    a = [];
    for j = 1: length(vect_reac)    %loop on the indeces of the reactions
        b = vect_reac{j};
        c = unique(b(b>0));
        for k = 1: length(c) %loop on the indeces of the products of the reaction of index j
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
%==========================================================================
%COMPUTATION OF THE REGULAR TRANSITIONS
%cell2mat(trans_set)

trans_set_preannot = [];
for i = 1:length(trans_set(:,1))
    a = repmat(trans_set(i,:),length(annot_space(:,1)),1);
    b = [a, annot_space, zero_annot_space];
    trans_set_preannot = [trans_set_preannot; b]; %4th column: annotation of the state; last column: annotation of the successor set to 0
end

cell2mat(trans_set_preannot(:,[1 4 2]));

%--------------------------------------------------------------------------
%Keeping only the regular transitions for which none of the updated
%variables are annotated and updating of the annotation of the successors
%of the allowed regular transitions

trans_set_reg = trans_set_preannot;
list = [];
for i = 1 : length(trans_set_preannot(:,1))    
    res = reg_trans_filt(trans_set_preannot(i,1), trans_set_preannot(i,4), trans_set_preannot(i,2));
    bool = res.a;
    annot_succ = res.b;
    trans_set_reg(i,5) = {res.b};  %updating of the annotation of the successor of the authorized regular transitions
    if res.a == 1         %the regular transition is allowed                     
        list = [list, i]; %storage of the indices of the authorized regular transitions
    end
end
trans_set_reg = trans_set_reg(list,:);

cell2mat(trans_set_reg(:, [1 4 2 5]));

%==========================================================================
%COMPUTATION OF THE SINGULAR TRANSITIONS
%Computation of the singular transitions restricted to the states with no 
%annotation which can be updated by a regular transition (states
%corresponding to the 1st column of trans_set_preannot)

states_trans_set_preannot = trans_set_preannot(:,[1 4]); %extraction of the states and their annotations from trans_set_preannot
a = cell2mat(states_trans_set_preannot);
[b, ia, ib] = unique (a,'rows');
states_trans_set_preannot = states_trans_set_preannot(ia,:); %suppression of the redondant annotated states

trans_set_sing_prefilt = []; 
for i = 1:length(states_trans_set_preannot(:,1))
    annot = cell2mat(states_trans_set_preannot(i,2));
    for j = 1:length(annot)
        if annot(j) == -1
            b = annot;
            b(j) = 0;  %loss of the annotation of the variable j
            trans_sing = [states_trans_set_preannot(i,1) states_trans_set_preannot(i,1) {0} states_trans_set_preannot(i,2) {b}];  %reaction index set to 0 by default
            trans_set_sing_prefilt = [trans_set_sing_prefilt; trans_sing];
        end
    end
end

cell2mat(trans_set_sing_prefilt(:,[1 4 5]));

%--------------------------------------------------------------------------
%Keeping only the singular transitions which are allowed by the upwards
%crossing interval constraint
list = [];
for i = 1: length(trans_set_sing_prefilt(:,1))    %loop on the transition of trans_set_preannot
    res = sing_trans_filt(trans_set_sing_prefilt(i,1), trans_set_sing_prefilt(i,4), trans_set_sing_prefilt(i,5), react_for_a_prod, vect_reac);
    if res == 1
        list = [list i];
    end    
end
trans_set_sing = trans_set_sing_prefilt(list,:);

cell2mat(trans_set_sing(:,[1 4 5]));

res = state_quotienting(prod_ind, state_space, state_space_annot, trans_set_sing_prefilt, trans_set_sing, trans_set_reg);

%==========================================================================
%QUOTIENTING OF ANNOTATED STATES
%==========================================================================

function res = state_quotienting(prod_ind, state_space, state_space_annot, trans_set_sing_prefilt, trans_set_sing, trans_set_reg)
%Quotienting of the annotated states for which no singular transitions involving
%a variable which can be positively updated by a regular transition have been
%suppressed by the upwards crossing interval constraint

%--------------------------------------------------------------------------
%Number of singular transitions without taking into account the upwards
%crossing interval constraint
%numb_var = length(state_space(1,:));
nb_annot_var = length(prod_ind);

%--------------------------------------------------------------------------
%DETERMINATION OF THE STATES TO QUOTIENT
states_set_sing = cell2mat(trans_set_sing_prefilt(:,1));
states_set_sing_uniq = unique(states_set_sing,'rows');   %extraction of the initial states of the transitions stored in trans_set_sing

%Extraction of the set of the annotations of the singular transitions
%before filtering
b = (cell2mat(trans_set_sing_prefilt(:,1)) - repmat(states_set_sing_uniq(1,:),[length(cell2mat(trans_set_sing_prefilt(:,1))), 1]));
c = find(sum(abs(b')) == 0);    
d_prefilt = trans_set_sing_prefilt(c,[4 5]);

list1 = [];
for i = 1:length(states_set_sing_uniq(:,1))

    state = states_set_sing_uniq(i,:);
    
    %Determination of the variables of 'state' which are not positively updated by a
    %regular transition
    b = cell2mat(trans_set_reg(:,1)) - repmat(state,[length(trans_set_reg(:,1)), 1]);
    c = find(sum(abs(b')) == 0);
    cell2mat(trans_set_reg(c,1:2));

    list = [];
    for j = 1:length(c)
        d = cell2mat(trans_set_reg(c(j),2)) - cell2mat(trans_set_reg(c(j),1));
        e = find(d <= 0);
        list = [list e];       
        list = intersect(list, e);  %list of the variables which are not positively updated by a regular transition
    end
    list = unique(list);
    var_annot_not_posit_updat = intersect(list, prod_ind);  %list of the annotated variables which are not updated by a regular transition

    %Computation of the number of singular transitions which involve
    %variables which can be positively updated by a regular transition before filtering with the
    %uci constraint
    iter1 = 0;
    for j = 1:length(d_prefilt(:,1))
        e = find([cell2mat(d_prefilt(j,2)) - cell2mat(d_prefilt(j,1))] ~= 0);
        f = intersect(e, var_annot_not_posit_updat);
        if length(f) == 0
            iter1 = iter1 + 1;
        end        
    end
    
    %Extraction of the set of the annotations of the singular transitions
    %of state 'state' after filtering
    b = (cell2mat(trans_set_sing(:,1)) - repmat(state,[length(trans_set_sing(:,1)), 1]));
    c = find(sum(abs(b')) == 0);
    d_postfilt = trans_set_sing(c,[4 5]);

    %Computation of the number of singular transitions which involve
    %variables which can be updated by a regular transition after filtering with the
    %uci constraint
    iter2 = 0;
    for j = 1:length(d_postfilt(:,1))
        e = find([cell2mat(d_postfilt(j,1)) - cell2mat(d_postfilt(j,2))] ~= 0);
        f = intersect(e,var_annot_not_posit_updat);
        if length(f) == 0
            iter2 = iter2 + 1;
        end
    end

    if iter1 - iter2 == 0
        list1 = [list1; state];
    end   
end

states_set_reg = cell2mat(trans_set_reg(:,1));
states_set_reg_uniq = unique(states_set_reg,'rows');   %extraction of the initial states of the transitions stored in trans_set_sing
num_state_space = cell2mat(state_space);

[d,ia,ib] = intersect(states_set_reg_uniq, num_state_space,'rows');
dd = num_state_space;
dd(ib,:) = [];
states_no_outgo_trans = dd;

% Add to list1 the states which do not have outgoing transitions (states which
% belong to state_space but not to a_uniq)
list1 = [list1; states_no_outgo_trans];

%--------------------------------------------------------------------------
%Suppression of the singular transitions in trans_sing for which its state
%is in the array list1
a = cell2mat(trans_set_sing(:,1));
list2 = [];
for i = 1:length(list1(:,1))
    b = a - repmat(list1(i,:),[length(a(:,1)), 1]);
    c = find(sum(abs(b')) == 0);
    list2 = [list2 c];
end

trans_set_sing(list2,:) = [];

cell2mat(trans_set_sing(:,1));

cell2mat(state_space_annot);
%--------------------------------------------------------------------------
%Quotienting of rows in trans_reg and of the state space
cell2mat(trans_set_reg);
for i = 1:length(list1(:,1))
    
    %Quotienting of the state space
    numb_var = length(cell2mat(trans_set_reg(1,1)));
    a = cell2mat(state_space_annot(:,1:numb_var));
    a = a(:,1:numb_var);
    b = a - repmat(list1(i,:),[length(a(:,1)), 1]);
    c = find(sum(abs(b')) == 0);
    for j = 1:length(c)
        d = cell2mat(state_space_annot(c(j),:));
        d(numb_var+1:end) = [zeros(1,numb_var) + 1];
        d = num2cell(d);
        state_space_annot(c(j),:) = d;
    end

    %Quotienting of rows in trans_reg 
    a = cell2mat(trans_set_reg(:,1));
    b = a - repmat(list1(i,:),[length(a(:,1)), 1]);
    c = find(sum(abs(b')) == 0);
    for j = 1:length(c)
        trans_set_reg(c(j),4) = {[zeros(1,length(cell2mat(trans_set_reg(1,1)))) + 1]};
    end
    
    a = cell2mat(trans_set_reg(:,2));
    b = a - repmat(list1(i,:),[length(a(:,1)), 1]);
    c = find(sum(abs(b')) == 0);
    for j = 1:length(c)
        trans_set_reg(c(j),5) = {[zeros(1,length(cell2mat(trans_set_reg(1,1)))) + 1]};
    end

end

a = cell2mat(trans_set_reg);
[b, ia,ib] = unique(a, 'rows');
trans_set_reg = trans_set_reg(ia,:); %Elimination of the redondant transitions

cell2mat(trans_set_reg);

a = cell2mat(state_space_annot);
[b, ia,ib] = unique(a, 'rows');
state_space_annot = state_space_annot(ia,:); %Elimination of the redondant transitions

cell2mat(state_space_annot);

cell2mat(trans_set_sing);

res.ts = [trans_set_reg; trans_set_sing];
res.ss = state_space_annot;

function res = reg_trans_filt(state, state_annot, succ_state)
%Filtering of a regular transition + annotation of the allowed
%successor of a regular transition

%%% A non annotated variable is annotated after a regular transition if it 
%%% is positively updated after the transition.
%%% An annotated variable keeps annotated after a regular transition if it 
%%% is not negatively updated after the transition.

bool = 1;
annot_succ = zeros(1,length(cell2mat(state)));

%--------------------------------------------------------------------------
%Filtering of the regular transition
annot_var = find(cell2mat(state_annot));   %extraction of the indeces of the variables which are annotated
delta = cell2mat(succ_state) - cell2mat(state);
for i = 1:length(delta)
    if delta(i) > 0
        if length(intersect(i, annot_var)) ~= 0  % check whether the updated variable of index i is annotated
            bool = 0; break;
        end
    end
end

%--------------------------------------------------------------------------
%Annotation of the successor if the regular transition is allowed
if bool == 1
    b1 = []; b2 = [];
    for i = 1:length(delta)
        if delta(i) > 0
            b1 = [b1, i]; %extraction of the indeces of the positively updated variables in the transition state -> succ_a
        end
    end

    if length(annot_var) == 0
    else
        for i = 1: length(annot_var)
            if delta(annot_var(i)) == 0
                b2 = [b2, annot_var(i)]; %extraction of the indeces of the annotated variables which are not updated in the transition state -> succ_a
            end
        end
    end
    annot_succ(b1) = -1; annot_succ(b2) = -1;
end

res.a = bool;
res.b = annot_succ;

function res = sing_trans_filt(state, annot_state, annot_succ, react_for_a_prod, vect_reac)
%Filtering of the singular transition according to the upwards crossing
%interval constraint

bool = 1;
num_state = cell2mat(state);
num_annot_state = cell2mat(annot_state);
num_annot_succ = cell2mat(annot_succ);

mrp = max_resource_prod(state, react_for_a_prod, vect_reac);  %computation of the max of resources available for the production of each product (assuming none of the limiting reactants are product of a reaction) 

index_updated_state = find(num_annot_state - num_annot_succ);
a = find([cell2mat(mrp(:,1)) - index_updated_state] == 0);
mrp_us = mrp(a,:);  %max ressource available for the production of the updated annotated product in the singular transition

%--------------------------------------------------------------------------
%Filtering of the singular transition (this encodes an approximation of the
%parameter: escape#)
products = cell2mat(react_for_a_prod(:,1)); %indeces of the products of the reaction network
if length(intersect(cell2mat(mrp_us(3)), products)) == 0      %1st condition to suppress a singular transition: none of the reactants of min reactants value producing the updated annotated species are a product of a reaction
    if num_state(cell2mat(mrp_us(1))) > cell2mat(mrp_us(2))   %2nd condition to suppress a singular transition: there is not enough resources for the updated annotated species to cross his interval upwards (i.e. to lose its annotation)
        bool = 0;
    end
end

res = bool;

function res = max_resource_prod(state, react_for_a_prod, vect_reac)
%Function which computes the resource available for the production of each
%product at a given state, and the indeces of the reactants which equal the
%min of the reactants of the reactions which produce each product
res = [];
for i = 1:length(react_for_a_prod(:,1))
    %react_for_a_prod(i,1)
    aa = cell2mat(react_for_a_prod(i,2)); %indeces of the reactions which produce the product of index: react_for_a_prod(i,1)
    b = []; c = [];
    for j = 1:length(aa)
        vv = vect_reac{aa(j)};
        vvv = unique(abs(vv(vv<0)));   %indeces of the reactants of the reaction vv{aa(j)}        
        num_state = cell2mat(state);
        b = [b, min(num_state(vvv))];  %vector of the min of the reactants of the reactions producing the product of index: react_for_a_prod(i,1)
        for k = 1: length(vvv)
            d = num_state(vvv);
            if d(k) == min(d)
                c = [c, vvv(k)];       %indeces of the reactants which equal the min of the reactants of the reactions producing the product of index i
            end
        end       
    end
    c = unique(c);
    e = [react_for_a_prod(i,1), {max(b)}, {c}];
    res = [res; e]; % 1st column: indeces of the products; 2nd column: maximum of the minimum of the reactant values over the reactions which produce the product of 1st column index ; 
                    % 3rd column: indeces of the reactants which equal
                    % the min of the reactants of the reactions producing
                    % the product of 1st column index
end