function frontend(model)
%function res = frontend(model,init_cond,mass_invar)

%file = strcat('.\input\',model,'\');
file = strcat('./input/',model,'/');

%%%========================================================================
%%% Parsing of the model file to extract the reactions
%%%========================================================================
fid = fopen(strcat(file,strcat(model),'_reactions.txt'));
v = [];
%var = [];
tline = fgetl(fid);
while ischar(tline)
    a = [];
    [d, b] = strtok(tline);
    % Parsing of the reactants 
    while d ~= '->'
        if d ~= '+' 
            a = [a -str2num(d)];
        end
        [d, b] = strtok(b);
    end
   
    % Parsing of the products
    [d, b] = strtok(b); % suppress the character '->'

    while length(d) ~= 0
        if d ~= '+'
            a = [a str2num(d)];
        end
        [d, b] = strtok(b);
    end
    v = [v , {a}];   % List of the reactions (parsed)
    tline = fgetl(fid);
end
fclose(fid);

%%%========================================================================
%%% Parsing of the variablels
%%%========================================================================

fid = fopen(strcat(file,strcat(model),'_reactions_alphabet.txt'));
var = [];
tline = fgetl(fid);
while ischar(tline)
    [d, b] = strtok(tline);
    while d ~= '='
        [d, b] = strtok(b);
    end
    b = strtrim(b);
    var = [var {b}];
    tline = fgetl(fid);
end
fid(close);

%%%========================================================================
%%% Parsing of the initial conditions file
%%%========================================================================
fid = fopen(strcat(file,'init_cond.txt'));

i_c_pars = [];
tline = fgetl(fid);
while ischar(tline)
    [d, b] = strtok(tline);
    % Parsing of the initial conditions 
    a = [];
    while d ~= '='
        a = {d};
        [d, b] = strtok(b);
    end
   % Parsing of the products
    [d, b] = strtok(b); % suppress the character '='        
    i_c_pars = [i_c_pars; {a} {d}];
    tline = fgetl(fid);
end

%List of the initial conditions in the input format of the core program
i_c = zeros(1,length(var));
size_ic = size(i_c_pars);
for i = 1: size_ic(1)
    b = i_c_pars{i,1};
    for k = 1:length(var)
        c = [b, var(k)];
        c = unique(c);
        if length(c) ==  1
            i_c(k) = str2num(i_c_pars{i,2});
            break;
        end    
    end
end
i_c = num2cell(i_c);
fclose(fid);

%%%========================================================================
%%% Parsing of the mass invariant file
%%%========================================================================
fid = fopen(strcat(file,'mass_invar.txt'));

MI_pars = [];
tline = fgetl(fid);
while ischar(tline)
    [d, b] = strtok(tline);
    % Parsing of the mass invariants
    a = [];
    while d ~= '='
        if d ~= '+' 
            a = [a {d}];           
        end
        [d, b] = strtok(b); % b: character left after suppression of the firts character - d: first character
    end
    
    if length(b) == 0
    else
        while length(b) ~= 0
            [d, b] = strtok(b);
            if d ~= '+'
                a = [a {strcat('-',d)}];
            end
        end    
    end
    MI_pars = [MI_pars; {a}];
    tline = fgetl(fid);
end

fclose(fid);

%List of the mass invariants in the input format of the core program
size_MI = size(MI_pars);
MI = zeros(size_MI(1),length(var));
for i = 1: size_MI(1)
    b = MI_pars{i};
    for j = 1:length(b)
        for k = 1:length(var)
            c = [];
            c = [b(j), var(k)];
            c = unique(c);
            if length(c) ==  1
                MI(i,k) = 1;
                break;
            else
                c = [b(j), strcat('-',var(k))];
                c = unique(c);
                if length(c) ==  1
                    MI(i,k) = -1;
                    break;
                end
            end
        end    
    end
end
MI = num2cell(MI);

%==========================================================================
%%% Loading of the number of samplibg intervals and of the reaction rates
%==========================================================================

nb_val = load(strcat(file,'nb_sampl_interv.txt'));

kin = load(strcat(file,'reaction_rates.txt'));

%--------------------------------------------------------------------------

res.v = v;
res.var = var;
res.i_c = i_c;
res.MI = MI;

main_v1(i_c,MI,v,kin,nb_val);


% % Construction of the model from the text file
% m = sbiomodel(model);
% for i = 1: length(C)
%     m.addreaction(C{i});
% end
% 
% % Computation of the stoechiometric matrix using matlab function
% [M,sp_stoech,reactions] = getstoichmatrix(m);
% S = zeros(length(sp_stoech),length(reactions));
% S = S + M;

% % Computation of the matrix of mass invariant with positive coefficients
% [g sp_mat_inv] = sbioconsmoiety(m,'semipos');

%res.a = g;
%res.b = sp_mat_inv;
%res.S = S;
%res.sp_stoec = sp_stoech;
