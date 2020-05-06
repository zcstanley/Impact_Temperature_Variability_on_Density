function [lg] = len_grad(T)
% input: mean temperature, 3d field
% output: centered differences squared, 3d field
% first and last rows are NaNs because differences in y derivatives cannot be calculated there.

  lg = 0.25 .* ( circshift(T, [0 1 0]) - circshift(T, [0 -1 0 ]) ).^2 + 0.25 .* ( circshift(T, [1 0 0]) - circshift(T, [-1 0 0]) ).^2 ;
  lg(1,:, :) = NaN ;
  lg(size(T,1), :, :) = NaN ;

end