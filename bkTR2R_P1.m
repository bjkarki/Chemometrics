function output = bkTR2R_P1 (Xcaluc, ycaluc, xvaluc, yvaluc, lambdas, etas, n)
% bkTR2R_P1 performs Residual Twonorm Tikhonov Regularization (TR2-R) using Process 1 on a given data set.
% In this process, Rval is used for TR2-R procedure.
%
% Syntax:
% output = bkTR2R_P2 (Xcaluc, ycaluc, xvaluc, yvaluc, lambdas, etas, n)
%
% Input arguments:
% Xcaluc = Calibration Martix in primary condition (un-mean centered).
% ycaluc = Analyte in primary condition (un-mean centered).
% xvaluc = A Validation sample in secondary condition (un-mean centered).
% yvaluc = Analyte in secondary condition (un-mean centered).
% lambdas = The first tuning meta-parameters (tunes  Identity matrix).
% etas = The second tuning meta-parameters (tunes  Residual matrix).
% n = no. of eigenvectors to be used. 

% Initial verification step
if nargin > 7
    error ('Too many input arguments')
elseif nargin < 7
    error ('Not enough input arguments')
end

% Validation sample check
if size (xvaluc, 1) ~= 1
    if size (xvaluc, 2) ~= 1
    error ('Only 1 validation sample can be used at a time')
    end
end

% Save input arguments for future purpose
output.Xcaluc = Xcaluc; 
output.ycaluc = ycaluc;
output.xvaluc = xvaluc; 
output.yvaluc = yvaluc;
output.lambda = lambdas;
output.eta = etas;
output.no_of_eigs = n;

% mean of rcal for given no. of eigs
rcal = bkloocv1 (Xcaluc);
rcalm = mean (squeeze (rcal (:, n, :)), 2);
output.rcalm = rcalm;   
output.rcalm2norm = norm (rcalm);

% Mean-center X and y
[Xcal, xmean] = mncn (Xcaluc);
xval = scale (xvaluc, xmean);
[ycal, ymean] = mncn (ycaluc);
yval = scale (yvaluc, ymean);

% Introduction of Identity matrix
I = eye (size (Xcal, 2));

% SVD
[~, ~, V] = svd (Xcal);
Vreq = V (:, 1 : n);
clear V;

% rval
rval = (I - (Vreq * Vreq')) * xval';
Rval = diag (rval);

% Save rvals
output.rval = rval;
output.rval2norm = norm (rval);
output.Rval = Rval;

% xvalcal
xvalcal = (Vreq * Vreq') * xval';
output.xvalcal = xvalcal;


% Loop iteration to provide different weights to Residual matrix
for  i = 1 : length (etas)
    
    % Pick up specific nus
    eta = etas (i);
    
    %Loop iteration stabilize the inverse
    for ii = 1 : length (lambdas)
        
        % Pick up specific lambda
        lambda = lambdas (ii);
        
        % Rank check
        rankcheck (ii, i) = rank ((Xcal' * Xcal) + (lambda ^ 2 .* I) + (eta ^ 2 .* (Rval ^ 2)));   
        
        % bhat for specific lambda's
        if rankcheck (ii, i) == size (I, 2)
            
            % Estimation of bhat
            bhat (:, ii, i) = ((Xcal' * Xcal) + (lambda ^ 2 .* I) + (eta ^ 2 .* (Rval ^ 2))) \ (Xcal' * ycal);
            twonorm (ii, i) = norm (bhat (:, ii, i));
            
            % Minimazation equation: Rval*bhat
            Rb (:, ii, i) = Rval * bhat (:, ii, i);
            Rb2norm (ii, i) = norm (Rb (:, ii, i));
            
            % yr: yhat residual
            yr (ii, i) = rval' * bhat (:, ii, i);
            
            % brval and brcal
            brval (:, ii, i) =  ((rval * rval') / (rval' * rval))* bhat (:, ii, i);
            brcal (:, ii, i) =  ((rcalm * rcalm') / (rcalm' * rcalm))* bhat (:, ii, i);
            
            brval2norm (ii, i) = norm (brval (:, ii, i));
            brcal2norm (ii, i) = norm (brcal (:, ii, i));
            
            % bval and bcal
            bval (:, ii, i) = ((bhat (:, ii, i) * bhat (:, ii, i)') / (bhat (:, ii, i)' * bhat (:, ii, i))) * rval;
            bcal (:, ii, i) = ((bhat (:, ii, i) * bhat (:, ii, i)') / (bhat (:, ii, i)' * bhat (:, ii, i))) * rcalm;
            
            bval2norm (ii, i) = norm (bval (:, ii, i));
            bcal2norm (ii, i) = norm (bcal (:, ii, i));
            
            % bxvalcal and bxcal
            bxvalcal (:, ii, i) = ((xvalcal * xvalcal') / (xvalcal' * xvalcal)) * bhat (:, ii, i);
            bxcal (:, ii, i) = (Vreq*Vreq') * bhat (:, ii, i);
            
            bxvalcal2norm (ii, i) = norm (bxvalcal (:, ii, i));
            bxcal2norm (ii, i) = norm (bxcal (:, ii, i));
            
            % Re-estimate ycal to find RMSEC
            ycalhat = Xcal * bhat (:, ii, i);
            rmsec (ii, i) = rmse (ycal, ycalhat);
            
            % Evaluate regstats
            ycalhatuc = rescale (ycalhat, ymean);
            ycalstat = regstats (ycalhatuc, ycaluc, 'linear', {'beta', 'rsquare'});
            
            r2cal (ii, i) = ycalstat.rsquare;
            slopecal (ii, i) = ycalstat.beta (2);
            interceptcal (ii, i) = ycalstat.beta (1);
            
            %Estimate yval
            yvalhat (ii, i) = xval * bhat (:, ii, i);
            rmsev (ii, i) = rmse (yval, yvalhat (ii, i));
        
        else
            
            % Estimation of bhat
            bhat (1 : (size (Xcal, 2)), ii, i) = nan;
            twonorm (ii, i) = nan;
            
            % Minimazation equation: Rval*bhat
            Rb (1 : (size (Xcal, 2)), ii, i) = nan;
            Rb2norm (ii, i) = nan;
            
            % yr
            yr (ii, i) = nan;
            
            % brval and brcal
            brval (1 : (size (Xcal, 2)), ii, i) = nan;
            brcal (1 : (size (Xcal, 2)), ii, i) = nan;
            
            brval2norm (ii, i) = nan;
            brcal2norm (ii, i) = nan;
            
            % bval and bcal
            bval (1 : (size (Xcal, 2)), ii, i) = nan;
            bcal (1 : (size (Xcal, 2)), ii, i) = nan;
            
            bval2norm (ii, i) = nan;
            bcal2norm (ii, i) = nan;
            
            % bxvalcal and bxcal
            bxvalcal (1 : (size (Xcal, 2)), ii, i) = nan;
            bxcal (1 : (size (Xcal, 2)), ii, i) = nan;
            
            bxvalcal2norm (ii, i) = nan;
            bxcal2norm (ii, i) = nan;
            
            % Re-estimate ycal to find RMSEC
            rmsec (ii, i) = nan;
            
            % Evaluate regstats
            r2cal (ii, i) = nan;
            slopecal (ii, i) = nan;
            interceptcal (ii, i) = nan;
            
            %Estimate yval
            yvalhat (ii, i) = nan;
            rmsev (ii, i) = nan;
        
        end     % ends if statement
    
    end         % ends ii loop
    
end             % ends i loop

% Angle merit
output.cos_alpha = brval2norm ./ twonorm;
output.cos_beta = brcal2norm ./ twonorm;

output.cos_theta = bxvalcal2norm ./ twonorm;
output.cos_phi = bxcal2norm ./ twonorm;

% Saving output
output.rankcheck = rankcheck;

output.bhat = bhat;
output.twonorm = twonorm;

output.Rb = Rb;
output.Rb2norm = Rb2norm;

output.yr = yr;

output.brval = brval;
output.brcal = brcal;

output.brval2norm = brval2norm;
output.brcal2norm = brcal2norm;

output.bval = bval;
output.bcal = bcal;

output.bval2norm = bval2norm;
output.bcal2norm = bcal2norm;

output.bxvalcal = bxvalcal;
output.bxcal = bxcal;

output.bxvalcal2norm = bxvalcal2norm;
output.bxcal2norm = bxcal2norm;

output.rmsec = rmsec;

output.r2cal = r2cal;
output.slopecal = slopecal;
output.interceptcal = interceptcal;

output.yvalhat = yvalhat;
output.rmsev = rmsev;

% Saving files on the output
 save ('output', 'output');

end     % ends simulatedTR2R

% rmse function
function rmse_ = rmse (y, yhat)
% rmse calculates root mean square error.
%
% Syntax
% rmse_ = rmse (y, yhat)
%
% Input arguments:
% y =  measured value.
% yhat = prediction of the same value.
% PS - They both should either be mean-centered or unmean-centered.
%
% Output arguments:
% rmse_ = root mean square error value.

d = 0;
for i = 1 : length (y)
    b = y (i); 
    c = yhat (i);
    d = d + ((b - c) ^ 2);
end
rmse_ = sqrt(d / length (y));
end

% bkloocv1 function
function rcal = bkloocv1 ( Xcaluc)
% bkloocv perfroms leave one out cross validation and finds out residual cal vector.
%
% Syntax:
% rcal = bkloocv1 (Xcaluc)
% 
% Input arguments:
% Xcaluc = Unmean-centered Xcal
%
% Output arguments:
% rcal = rcal for each cal sample and for different no. of eigenvectors (wavelength x eigenvector x sample).


% Initial verification step
if nargin > 1
    error ('Error: Too many input arguments')
elseif nargin < 1
    error ('Error: Not enough input arguments')
end

% Introduction of Identity matrix
I = eye (size (Xcaluc, 2));

% loop iteration to perfrom leave one out cross validation.
for i = 1 : size (Xcaluc, 1)
    
    % Pick out a cal sample
    Xcalrequc = Xcaluc;
    Xvalrequc = Xcalrequc (i, :);
    Xcalrequc (i, :) = [];
    
    % Mean Centering
    [Xcalreq, Xmean] = mncn (Xcalrequc);
    Xvalreq = scale (Xvalrequc, Xmean);
    
    % SVD Xcalsvd
    [~, ~, V] = svd (Xcalreq);
    
    
    % loop to get norm of rcal for different no. of eigenvectors.
    for ii = 1 : size (Xcalreq, 1)
        
        % estimation of rcal
        Vreq = V (:, 1 : ii);
        rcal (:, ii, i) = (I - (Vreq * Vreq')) * Xvalreq';
        
    end             % ends ii

end                 % ends i

end                 % ends function