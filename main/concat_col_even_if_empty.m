function res = concat_col_even_if_empty(a,b)
  if length(a) == 0 
    res = [b];
  else 
    res = [a, b];
  end
end


 
