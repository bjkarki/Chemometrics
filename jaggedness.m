function result = jaggedness (vector)
%
% Created by:
% Bibek Karki on 08/19/2014
% 
% jaggedness calculates jaggedness of given vector
%
% Syntax:
% result = jaggedness (vector)
%
% Input arguments:
% vector = a vector whose jaggedness is to be calculated.

% variable initialization
jag = 0;

% Loop iteration to pick right vector elements
for j = 2 : length (vector)
    
    % neighbour index of j
    i = j - 1;
    
    % summing of the values
    jag = jag + (vector(j) - vector(i))^2;
end

% jaggedness
result = sqrt (jag);

end