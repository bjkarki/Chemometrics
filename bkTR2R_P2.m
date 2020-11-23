function output = bkTR2R_P2 (Xcaluc, ycaluc, xvaluc, yvaluc, lambdas, etas, ieig, feig)

% bkTR2R_P2 performs Residual Twonorm Tikhonov Regularization (TR2-R) using Process 2 on a given data set.
% In this process, Rval is used for TR2-R procedure.
%
% Syntax:
% output = bkTR2R_P2 (Xcaluc, ycaluc, xvaluc, yvaluc, lambdas, etas, ieig, feig)
%
% Input arguments:
% Xcaluc = Calibration Martix in primary condition (un-mean centered).
% ycaluc = Analyte in primary condition (un-mean centered).
% xvaluc = A Validation sample in secondary condition (un-mean centered).
% yvaluc = Analyte in secondary condition (un-mean centered).
% lambdas = The first tuning meta-parameters (tunes  Identity matrix).
% etas = The second tuning meta-parameters (tunes  Residual matrix).
% ieig = starting no. of eigs
% feig = ending no. of eigs

% start time
start_time = clock;
disp (['start time: ', num2str(start_time (4)), ':', num2str(start_time (5))])

% Initial verification step
if nargin > 8
    error ('Too many input arguments')
elseif nargin < 8
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
output.starting_no_of_eigs = ieig;

% Mean-center X and y
[Xcal, xmean] = mncn (Xcaluc);
xval = scale (xvaluc, xmean);
[ycal, ymean] = mncn (ycaluc);
yval = scale (yvaluc, ymean);

% Introduction of Identity matrix
I = eye (size (Xcal, 2));

% ieig and feig verification

if ieig >= rank (Xcal)
    fprintf ('ieig should be less than %d.\n', rank(Xcal))
end

if feig > rank (Xcal)
    fprintf ('feig should be less than or equal to %d.\n', rank(Xcal))
end

for iii = ieig : feig

    % mean of rcal for given no. of eigs
    rcal = bkloocv1 (Xcaluc);
    rcalm = mean (squeeze (rcal (:, iii, :)), 2);
    output.rcalm (:, iii) = rcalm;   
    output.rcalm2norm (iii) = norm (rcalm);

    % SVD
    [U, S, V] = svd (Xcal);
    Ureq = U (:, 1 : iii); 
    Sreq = S (1 : iii, 1 : iii); 
    Vreq = V (:, 1 : iii);
    clear U S V;

    % rval
    rval = (I - (Vreq * Vreq')) * xval';
    Rval = diag (rval);

    % Save rvals and rcals
    output.rval (:, iii) = rval;
    output.rval2norm (iii) = norm (rval);
    output.Rval (:, :, iii) = Rval;

    % xvalcal and Xcalreq
    xvalcal = (Vreq * Vreq') * xval';
    Xcalreq = Ureq * Sreq * Vreq';

    % Save xvalcal and Xcalreq
    output.xvalcal (:, iii) = xvalcal;
    output.Xcalreq (:, :, iii) = Xcalreq;

    % Loop iteration to provide different weights to Residual matrix
    for  i = 1 : length (etas)

        % Pick up specific nus
        eta = etas (i);

        %Loop iteration to identify bhat for different values of lambda
        for ii = 1 : length (lambdas)

            % Pick up specific lambda
            lambda = lambdas (ii);

            % Rank check
            rankcheck (ii, i, iii) = rank ((Xcalreq' * Xcalreq) + (lambda ^ 2 .* I) + (eta ^ 2 .* (Rval ^ 2)));   

            % bhat for specific lambda's
            if rankcheck (ii, i, iii) == size (I, 2)

                % Estimation of bhat
                bhat (:, ii, i, iii) = ((Xcalreq' * Xcalreq) + (lambda ^ 2 .* I) + (eta ^ 2 .* (Rval ^ 2))) \ (Xcalreq' * ycal);
                twonorm (ii, i, iii) = norm (bhat (:, ii, i, iii));

                % Minimazation equation: Rval*bhat
                Rb (:, ii, i, iii) = Rval * bhat (:, ii, i, iii);
                Rb2norm (ii, i, iii) = norm (Rb (:, ii, i, iii));

                % yr: yhat residual
                yr (ii, i, iii) = rval' * bhat (:, ii, i, iii);

                % brval and brcal
                brval (:, ii, i, iii) =  ((rval * rval') / (rval' * rval))* bhat (:, ii, i, iii);
                brcal (:, ii, i, iii) =  ((rcalm * rcalm') / (rcalm' * rcalm))* bhat (:, ii, i, iii);

                brval2norm (ii, i, iii) = norm (brval (:, ii, i, iii));
                brcal2norm (ii, i, iii) = norm (brcal (:, ii, i, iii));

                % bval and bcal
                bval (:, ii, i, iii) = ((bhat (:, ii, i, iii) * bhat (:, ii, i, iii)') / (bhat (:, ii, i, iii)' * bhat (:, ii, i, iii))) * rval;
                bcal (:, ii, i, iii) = ((bhat (:, ii, i, iii) * bhat (:, ii, i, iii)') / (bhat (:, ii, i, iii)' * bhat (:, ii, i, iii))) * rcalm;

                bval2norm (ii, i, iii) = norm (bval (:, ii, i, iii));
                bcal2norm (ii, i, iii) = norm (bcal (:, ii, i, iii));

                % bxvalcal and bxcal
                bxvalcal (:, ii, i, iii) = ((xvalcal * xvalcal') / (xvalcal' * xvalcal)) * bhat (:, ii, i, iii);
                bxcal (:, ii, i, iii) = (Vreq*Vreq') * bhat (:, ii, i, iii);

                bxvalcal2norm (ii, i, iii) = norm (bxvalcal (:, ii, i, iii));
                bxcal2norm (ii, i, iii) = norm (bxcal (:, ii, i, iii));

                % Re-estimate ycal to find RMSEC
                ycalhat = Xcal * bhat (:, ii, i, iii);
                rmsec (ii, i, iii) = rmse (ycal, ycalhat);

                % Evaluate regstats
                ycalhatuc = rescale (ycalhat, ymean);
                ycalstat = regstats (ycalhatuc, ycaluc, 'linear', {'beta', 'rsquare'});

                r2cal (ii, i, iii) = ycalstat.rsquare;
                slopecal (ii, i, iii) = ycalstat.beta (2);
                interceptcal (ii, i, iii) = ycalstat.beta (1);

                %Estimate yval
                yvalhat (ii, i, iii) = xval * bhat (:, ii, i, iii);
                rmsev (ii, i, iii) = rmse (yval, yvalhat (ii, i, iii));

            else

                % Estimation of bhat
                bhat (1 : (size (Xcal, 2)), ii, i, iii) = nan;
                twonorm (ii, i, iii) = nan;

                % Minimazation equation: Rval*bhat
                Rb (1 : (size (Xcal, 2)), ii, i, iii) = nan;
                Rb2norm (ii, i, iii) = nan;

                % yr
                yr (ii, i, iii) = nan;

                % brval and brcal
                brval (1 : (size (Xcal, 2)), ii, i, iii) = nan;
                brcal (1 : (size (Xcal, 2)), ii, i, iii) = nan;

                brval2norm (ii, i, iii) = nan;
                brcal2norm (ii, i, iii) = nan;

                % bval and bcal
                bval (1 : (size (Xcal, 2)), ii, i, iii) = nan;
                bcal (1 : (size (Xcal, 2)), ii, i, iii) = nan;

                bval2norm (ii, i, iii) = nan;
                bcal2norm (ii, i, iii) = nan;

                % bxvalcal and bxcal
                bxvalcal (1 : (size (Xcal, 2)), ii, i, iii) = nan;
                bxcal (1 : (size (Xcal, 2)), ii, i, iii) = nan;

                bxvalcal2norm (ii, i, iii) = nan;
                bxcal2norm (ii, i, iii) = nan;

                % Re-estimate ycal to find RMSEC
                rmsec (ii, i, iii) = nan;

                % Evaluate regstats
                r2cal (ii, i, iii) = nan;
                slopecal (ii, i, iii) = nan;
                interceptcal (ii, i, iii) = nan;

                %Estimate yval
                yvalhat (ii, i, iii) = nan;
                rmsev (ii, i, iii) = nan;

            end     % ends if statement

        end         % ends ii loop

    end             % ends i loop
    
    a = clock;
    fprintf ('%d of %d outer loop completed. \nDate: %d/%d/%d \tTime : %d:%d \n', iii, feig, a(1), a(2), a(3), a(4), a(5))
        
end                 % ends iii loop

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
 save ('output', 'output', '-v7.3');
 
 % stop time
 stop_time = clock;
 disp (['stop time: ', num2str(stop_time (4)), ':', num2str(stop_time (5))])

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