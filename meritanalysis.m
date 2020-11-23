function  merit = meritanalysis(output, lambda, n)
%
% Created by:
% Bibek Karki on 08/06/2014
% 
% 
% meritanalysis calculates ratios of different merits for analysis
% Syntax:
% merit = meritanalysis(output, lambda, n)
%
% Input arguments:
% Xcaluc = Calibration Matrix in primary condition (un-mean centered).
% lambda = initial lambda selection for picking the model
% n = no. of eigs to be used.

% start time
start_time = clock;
disp(['start time: ', num2str(start_time (4)), ' : ', num2str(start_time (5)) ])

% Initial verification step
if nargin > 3
    error ('Too many input arguments')
elseif nargin < 3
    error ('Not enough input arguments')
end

% Picking up the right merit from output file
%--------------------------------------------------------------------------------------------------------------------------------------------------
[merit.brval2norm, i] = maximum (output.brval2norm (lambda, :, n));

merit.Rb = output.Rb (:, lambda, i(2), n);

merit.rmsev = output.rmsev (lambda, i(2), n);

merit.bhat = output.bhat (:, lambda, i(2), n);

merit.twonorm = output.twonorm (lambda, i(2), n);

merit.bxvalcal2norm = output.bxvalcal2norm (lambda, i(2), n);

merit.cosAlpha = output.cos_alpha (lambda, i(2), n);

merit.rval = output.rval (:, n);

merit.residual2norm_ratio = output.rval2norm (n) ./ output.rcalm2norm (n);
%--------------------------------------------------------------------------------------------------------------------------------------------------



% Calculating new merits for post analysis
%--------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating ||brval|| over ||rval|| to estimate the magnitude of projection.
merit.brvalOverR = merit.brval2norm ./ output.rval2norm (n);

% Calculating ||rval|| over mean ||rcal|| to estimate how big rval is.
merit.rvalOverRcal = output.rval2norm (n) ./ output.rcalm2norm (n);

% Monitor the ||Rb|| obtained from cal set using the same model vector.
[merit.RbRatioUsingMean, merit.Rbcal_mean, merit.RbRatioUsingMax, merit.Rbcal_max] = Rbloocv1 (output.Xcaluc, output.xvaluc, merit.bhat, n);

% Angle between Rbval and Rbcal_mean
merit.angleRbMean = (merit.Rbcal_mean' * merit.Rb) / (norm(merit.Rbcal_mean) * norm(merit.Rb));

% Angle between Rbval and Rbcal_max
merit.angleRbMax = (merit.Rbcal_max' * merit.Rb) / (norm(merit.Rbcal_max) * norm(merit.Rb));

% Jaggedness measurement
merit.jagg_Rbval = jaggedness (merit.Rb);
merit.jagg_Rbval_Rbcalmean = merit.jagg_Rbval ./ jaggedness (merit.Rbcal_mean);
merit.jagg_Rbval_Rbcalmax = merit.jagg_Rbval ./ jaggedness (merit.Rbcal_max);

%--------------------------------------------------------------------------------------------------------------------------------------------------

% stop time
stop_time = clock;
disp(['stop time: ', num2str(stop_time (4)), ' : ', num2str(stop_time (5)) ])

end

function [usingmean, mean_vector, usingmax, max_vector] =Rbloocv1 (Xcaluc, xvaluc, bhat, n)

% Initial verification step
if nargin > 4
    error ('Too many input arguments')
elseif nargin < 4
    error ('Not enough input arguments')
end

% Xcal and xval
[Xcal, xmean] = mncn (Xcaluc);
xval = scale (xvaluc, xmean);
clear xmean;

% Identity Matrix introduction
I = eye (size (Xcal, 2));

% val SVD
[~, ~, V] = svd (Xcal);
R = diag ((I - (V (:, 1 : n) * V (:, 1 : n)')) * xval');

% ||Rb|| of val
Rb2normval = norm (R * bhat);

% Clear repeating variables
clear R V;

for j = 1 : size (Xcaluc, 1);
    
    % Pickup a cal sample
    Xcalrequc = Xcaluc;
    xvalrequc = Xcalrequc (j, :);
    Xcalrequc (j, :) = [];
    
    % Mean centering
    [Xcalreq, xmean] = mncn (Xcalrequc);
    xvalreq = scale (xvalrequc, xmean);
    
    % cal SVD
    [~, ~, V] = svd (Xcalreq);
    R = diag ((I - (V (:, 1 : n) * V (:, 1 : n)')) * xvalreq');
   
    % Rb cal vector
    Rb (:, j) = R * bhat;
    
    % ||Rb|| of cal
    Rb2normcal (j) = norm (Rb(:, j));
    clear R V;
    
end         % ends j

% mean vector
mean_vector = mean (Rb, 2);

% max vector
[~, i] = maximum (Rb2normcal);
max_vector = Rb (:, i (2));

% mean and max ratio
usingmean = Rb2normval ./ mean (Rb2normcal);
usingmax = Rb2normval ./ max (Rb2normcal);

end

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
