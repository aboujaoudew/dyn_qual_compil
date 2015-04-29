%%%Inputs for the "rescaling" model

%% Stoechiometric matrix

S = [{-1} {0} {0} {-1}
     {-1} {-1} {0} {0}
     {0} {-1} {-1} {0}
     {1} {0} {-1} {0}
     {0} {1} {0} {-1}
     {0} {0} {1} {1}];

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

ci = zeros(1,16);
ci(7) = 1;

ci = num2cell(ci);

%% Matrix of mass invariant
MMI = zeros(4,16);
MMI(1,7:13) = 1; % Total amount of EGF moiety
MMI(2,[1:6 8:13]) = 1; % Total amount of EGFR moiety
MMI(3, [2 4:6 9:10 12:14 16]) = 1; % Total amount of SOS moiety
MMI(4,[3:6 10:13 15:16]) = 1; % Total amount of Shc moiety

MMI = num2cell(MMI);

%%%% Matrix of sign of the mass invariant matrix
size_MMI = size(MMI);
MMI_sign = [];
for i = 1:size_MMI(1)
    for j = 1:size_MMI(2)
        if MMI{i,j} > 0
            MMI_sign{i,j} = 1;
        elseif MMI{i,j} == 0
            MMI_sign{i,j} = 0;
        end
    end
end

%% Vector of reactions
 
v = [{[-8 1 7]},
     {[-9 2 7]},
     {[-10 4 7]}, 
     {[-11 3 7]}, 
     {[-12 6 7]},
     {[-13 5 7]},
     {[8 -1 -7]},
     {[9 -2 -7]},
     {[10 -4 -7]}, 
     {[11 -3 -7]}, 
     {[12 -6 -7]},
     {[13 -5 -7]},
     {[-8 -14 9]},
     {[-10 -14 13]},
     {[-11 -14 10]},
     {[-2 ]},
     {[]},
     ];
 
% a = [2, 5, 6, 9, 12, 13];
% b = [1, 3, 4, 8, 10, 11];
% 
% for i = 1:length(a)
%     for j = 1:length(b)
%         v = [v, {[-]}]
% 
%     end
% end
 
S_sign{1,1}
v{1}