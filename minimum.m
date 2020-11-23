function [v, i] = minimum (X)
% minimum gives the minimum value of a matrix and its location.
%
% Syntax:
% [v, i] = minimum (X)
%
% Input Argument:
% X = the given matrix to be searched.
%
% Output Arguments:
% v = the minimum value in the matrix X.
% i = the location of minimum value.

v = min (min (X));      % minimum value of the given matrix.  

[r, c] = find (X == v);   % returns the index of minimum value.

i = [r, c];

end