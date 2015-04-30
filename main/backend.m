function backend(ci,S,P,FP,model)

size_S = size(S);
a = [];
for i=1:length(ci)
    a = strcat(a,int2str(ci(i)));
end
fid = fopen(strcat('../case_studies/',model,'/output_files/','ci',a,'.dot'), 'wt');

fprintf(fid, '%s\n\n', 'digraph G{');
fprintf(fid, '%s', '{');

for i = 1:size_S(1)
    fprintf(fid, '%s', strcat('node_',int2str(i),' [label="'));
    for j = 1:size_S(2)
        fprintf(fid, '%s', int2str(S(i,j)));
    end
    fprintf(fid,'%s\n','"];');    
end

fprintf(fid, '%s\n', '}');

Q = P; 

%Replacement of the following propensities: [0 1] -> 0.5 - [0 1 2] -> 1 - [0 1 2 3] -> 1.5 

size_Q = size(Q); 
n = size_Q(1);
a = 0;
while n ~= 0
    a = a + 1;
    k = 0; b = [];
    for i = 2:100
        if Q{i,1} ~= -1
            for j = 1:size_S(1)
                if cell2mat(Q(i,1)) - S(j,:) == 0
                    if cell2mat(Q(i,2)) == 0
                        % fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                        % 'node_', int2str(j),'[penwidth=',int2str(1),'; style=dashed]')); % 0 -> dashed
                    break;
                    elseif cell2mat(Q(i,2)) == 1
                        fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                        'node_', int2str(j),'[penwidth=',int2str(cell2mat(Q(i,2))),']')); % 1 -> 1
                    break;
                    elseif cell2mat(Q(i,2)) == 2
                        fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                        'node_', int2str(j),'[penwidth=',int2str(cell2mat(Q(i,2))),']')); % 2 -> 2
                    break;
                    elseif cell2mat(Q(i,2)) == 3
                        fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                        'node_', int2str(j),'[penwidth=',int2str(cell2mat(Q(i,2))),']')); % 3 -> 3
                    break;
                    elseif length(cell2mat(Q(i,2))) == 2    
                        if cell2mat(Q(i,2)) == [0 1]
                            fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                            'node_', int2str(j),'[penwidth=',int2str(0.5),']')); % [0 1] -> 0.5
                            break;
                        end
                    elseif length(cell2mat(Q(i,2))) == 3    
                        if cell2mat(Q(i,2)) == [0 1 2]
                            fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                            'node_', int2str(j),'[penwidth=',int2str(1),']')); % [0 1 2] -> 1
                            break;
                        end
                    elseif length(cell2mat(Q(i,2))) == 4    
                        if cell2mat(Q(i,2)) == [0 1 2 3]
                            fprintf(fid, '%s\n', strcat('node_',int2str(a),'->',...
                            'node_', int2str(j),'[penwidth=',int2str(1.5),']')); % [0 1 3] -> 1.5
                            break;
                        end    
                    end
                end
            end
            k = k + 1;
        else
            Q(1:k+2,:) = []; m = size(Q); n = m(1);
            break;
        end
    end
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
