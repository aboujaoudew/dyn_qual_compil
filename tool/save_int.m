function save_int(ci,res)

P = res.P;
FP = res.FP; a = sum(size(FP));
ns = res.ns;

fid = fopen(strcat('.\result\','report_','ci',int2str(ci),'.txt'), 'wt');

fprintf(fid, '%s\n\n', 'Number of states:');
fprintf(fid, '%d\n\n', ns);

fprintf(fid, '%s\n\n', 'Number of fixed points:');
if a == 0
    fprintf(fid, '%d\n\n', 0);
else
    fprintf(fid, '%d', FP);
end

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