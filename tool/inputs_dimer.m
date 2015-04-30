%%%Inputs for the "rescaling" model

%% Stoechiometric matrix

S = [{-2} {1}
     {2} {-1}]

%%% Matrix of sign of the stoechiometric matrix
S_sign = [];
a = size(S);
for i = 1:a(1)
    for j = 1:a(2)
        if S{i,j} > 0
            S_sign{i,j} = 1;
        elseif S{i,j} < 0
            S_sign{i,j} = -1;
        elseif S{i,j} == 0
            S_sign{i,j} = 0;
        end
    end
end
%% Initial conditions
 
ci = [{2} {0}];

%% Matrix of mass invariant

MMI = [{1} {2}]

%%%% Matrix of sign of the mass invariant matrix
size_MMI = size(MMI);
MMI_sign = [];
for i = 1:size_MMI(1)
    for j = 1:size_MMI(2)
        if MMI{i,j} > 0
            MMI_sign{i,j} = 1;
        elseif MMI{i,j} == 0
            MMI_sign{i,j} = 0;
        elseif MMI{i,j} < 0
            MMI_sign{i,j} = -1;
        end
    end
end

%% Vector of reactions
 
v = [{[ 1 1 ]} {[ 2 ]}]
S_sign{1,1}
v{1}