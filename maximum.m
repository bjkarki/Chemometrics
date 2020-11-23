function [v, i] = maximum (X)
% maximum gives the maximum matrix element and its position.
%
% Syntax:
% [v, i] = maximum (X)
%
% Input Argument:
% X = the given matrix to be searched.
%
% Output Arguments:
% v = the minimum value in the matrix X.
% i = the location of minimum value.

v = max (max (X));      % minimum value of the given matrix.  

[r, c] = find (X == v);   % returns the index of minimum value.

i = [r, c];

end