function backend(ci, state_space, trans_set, model)

%--------------------------------------------------------------------------
%Creation of a folder "output_files", to stock the results, if it does not exist
newSubFolder = strcat('../case_studies/',model,'/output_files/');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end
%--------------------------------------------------------------------------

size_S = size(state_space);
ss = cell2mat(state_space);

a = [];
for i=1:length(ci)
    a = strcat(a,int2str(ci(i)));
end
fid = fopen(strcat('../case_studies/',model,'/output_files/','ci',a,'.dot'), 'wt');

fprintf(fid, '%s\n\n', 'digraph G{');
fprintf(fid, '%s', '{');

for i = 1:size_S(1)
    fprintf(fid, '%s', strcat('node_',int2str(i),' [label="'));
    a = zeros(1,size_S(2)/2) + 1;
    b = ss(i,:);
    if b(end - size_S(2)/2 + 1:end) - a == 0
        b = ss(i, 1:size_S(2)/2);
    end
    for j = 1:length(b)
        fprintf(fid, '%s', int2str(b(j)));
    end
    fprintf(fid,'%s\n','"];');
end

fprintf(fid, '%s\n', '}');

%cell2mat(trans_set(:,1:2))

for i = 1:length(trans_set(:,1))
    i;
    %cell2mat(trans_set(i,1))
    for j = 1:size_S(1)
        %ss(j,:)
        if cell2mat(trans_set(i,1)) - ss(j,:) == 0   

            a = strcat('node_',int2str(j),'->');
            
            break;
        else
        end
        
    end
    
    for j = 1:size_S(1)
        j;
        cell2mat(trans_set(i,2)) - ss(j,:);
        
        if cell2mat(trans_set(i,2)) - ss(j,:) == 0          
            b = strcat('node_', int2str(j),'[penwidth=','1',']');
            break;
        end
        
    end
   
    fprintf(fid, '%s\n', strcat(a,b));
    
end

% % Positioning the steady state at the bottom of the graph
% a = strcat(' node_',int2str(FP{1,2}));
% for i = 2:length(FP(:,1))
%     a = strcat(a,' node_',int2str(FP{i,2}));
% end
% fprintf(fid, '%s\n', strcat('{rank=same',a,'}'));
% 
fprintf(fid, '%s\n', '}');

fclose(fid);
