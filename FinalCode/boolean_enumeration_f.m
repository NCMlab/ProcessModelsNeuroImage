function combo_matrix = boolean_enumeration_f(nOfLayers)
%function combo_matrix = boolean_enumeration_f(nOfLayers)
%%%%% Construct a binary tree structure in a matrix form.
%%%%% First row represents one; last row represents (2^nOfLayers - 1).
%%%%% "nOfLayers" equals the number of columns.



combo_matrix = zeros(2^nOfLayers - 1, nOfLayers);

for combo_num=1:(2^nOfLayers - 1)
   temp_char_string = dec2bin(combo_num);

   for j=1:length(temp_char_string)
      if temp_char_string(j)=='1'
         combo_matrix(combo_num, nOfLayers - length(temp_char_string) + j) = 1;
      end
   end

end
