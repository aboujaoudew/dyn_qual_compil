function save_report(ci,res)

P = res.P;
FP = res.FP; size_FP = size(FP);
ns = res.ns;

fid = fopen(strcat('./result/','report_','ci',int2str(ci),'.txt'), 'wt');

fprintf(fid, '%s\n\n', 'Number of states:');
fprintf(fid, '%d\n\n', ns);

fprintf(fid, '%s\n\n', 'Fixed points:');
if sum(size_FP) == 0
    fprintf(fid, '%d\n\n', 'No fixed points');
else
    for i = 1:size_FP(1)
        fprintf(fid, '%d', FP(i,:));
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n');

fprintf(fid, '%s\n\n', 'States and their successors:');
n=size(P);
for i=1:n(1)
    if P(i,1) == -1
        fprintf(fid, '\n');
    else
        fprintf(fid, '%d ', P(i,:)'); fprintf(fid, '\n');    
    end
end
fclose(fid);