%% Computation of the mass invariant
function MI = mass_invar_R(ci,MMI_sign)

size_MMI_sign = size(MMI_sign);
%%% Computation of the mass invariant
b = [];
for i = 1:size_MMI_sign(1)
    for j = 1: size_MMI_sign(2)
        b{i,j} = [MMI_sign{i,j}*ci{j}];
    end
end
MI = [];
for i = 1:size_MMI_sign(1)
    % MI = [MI ; {fplus({b{i,:}})}];
    MI = [MI ; {fplus_R({b{i,:}})}];
end