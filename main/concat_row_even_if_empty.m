function res = concat_row_even_if_empty(a,b)
  if length(a) == 0
    res = [b];
  else 
    %a,b
    res = [a; b];
  end
end


 
