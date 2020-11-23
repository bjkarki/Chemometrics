function correlation_plot (Xcaluc, xvaluc)
% correlation_plot gives 2 different correlation plots that will help determine no. of eigenvectors to be used.
%
% Syntax:
% correlation_plot (Xcaluc, xvaluc)
%
% Input arguments:
% Xcaluc = Calibration Samples in primary condition (un-mean centered).
% xvaluc = a validation sample in secondary condition (un-mean centered).

% Obtaining mean of correlation and its std. deviation for all eigs
for i = 1 : rank (mncn (Xcaluc))
[mean_cor(i), std_cor(i)] = correlation_val (Xcaluc, xvaluc, i);
[mean_cor_cal(i), std_cor_cal(i)] = correlation_cal (Xcaluc, i);
end

% Correlation vs. eigs plot
figure, plot (mean_cor)
xlabel ('# of eigs')
ylabel ('mean of correlation')
title ('mean of correlation plot')
axis tight

% Correlation over std. deviation vs. eigs plot
figure, plot (mean_cor ./ std_cor)
xlabel ('# of eigs')
ylabel ('mean of correlation / std. deviation')
title ('mean of correlation / std. deviation plot')
axis tight

% Correlation of cal as a func of eig
figure, plot (mean_cor_cal)
xlabel ('# of eigs')
ylabel ('mean of correlation')
title ('mean of correlation CAL plot')
axis tight

% Correlation over std. dev of cals a s a func of eig
figure, plot (mean_cor_cal ./ std_cor_cal)
xlabel ('# of eigs')
ylabel ('mean of correlation / std. deviation')
title ('mean of correlation / std. deviation CAL plot')
axis tight

end

function [mean_cor, std_cor] = correlation_val (Xcaluc, xvaluc, eig)
% correlation_val gives correlation between the cal and val set for given no. of eigenvectors.
%
% Syntax:
% [cor_vec, xvalcal] = correlation_val (Xcaluc, xvaluc, eig)
%
% Input arguments:
% Xcaluc = Calibration Samples in primary condition (un-mean centered).
% xvaluc = a validation sample in secondary condition (un-mean centered).
% eig = no. of eigenvectors to be used.
%
% Output arguments:
% cor_vec = correlation vector
% xvalcal = projection vector of xval onto cal space (mean centered).

% Initial verification step
if nargin > 3
    error ('Too many input arguments')
elseif nargin < 3
    error ('Not enough input arguments')
end

% Validation sample check
if size (xvaluc, 1) ~= 1
    if size (xvaluc, 2) ~= 1
    error ('Only 1 validation sample can be used at a time')
    end
end

% Mean-center X
[Xcal, xmean] = mncn (Xcaluc);
xval = scale (xvaluc, xmean);

% keep track of no. of rows
[r, ~] = size (Xcal);

% SVD
[~, ~, V] = svd (Xcal);
Vreq = V (:, 1 : eig);
clear U S V;

% xvalcal
xvalcal = ((Vreq * Vreq') * xval')';

% Building Xvalcal matrix
Xvalcal = [];
for i = 1 : r
    Xvalcal = [Xvalcal; xvalcal];
end

% Transpose of each matrix because corr does pairwise correlation of columns
Xcal = Xcal'; 
Xvalcal = Xvalcal'; 


% correlation vector of samples
cor_vec = corr (Xcal, Xvalcal);
cor_vec = abs (diag (cor_vec));

% Print mean and std. deviation so we have better understanding of correlation.
mean_cor = mean (cor_vec);
std_cor = std (cor_vec);

end

function [mean_cor, std_cor] = correlation_cal (Xcaluc, n)
% correlation_cal gives correlation between the cal space formed using
% given no. of eigenvectors and the actual cal space.
%
% Syntax:
% cor_vec = correlation_cal (Xcaluc, n)
%
% Input arguments:
% Xcaluc = Calibration Samples in primary condition (un-mean centered).
% n = no. of eigenvectors to be used.
%
% Output arguments:
% cor_vec = correlation vector

% Initial verification step
if nargin > 2
    error ('Too many input arguments')
elseif nargin < 2
    error ('Not enough input arguments')
end

% Mean-center X
[Xcal, ~] = mncn (Xcaluc);

% SVD
[U, S, V] = svd (Xcal);
Ureq = U (:, 1 : n); 
Sreq = S (1 : n, 1 : n); 
Vreq = V (:, 1 : n);
clear U S V;

% xvalcal
Xcalreq = Ureq * Sreq * Vreq';

% Transpose of each matrix because corr does pairwise correlation of columns
Xcal = Xcal'; 
Xcalreq = Xcalreq'; 


% correlation vector of samples
cor_vec = corr (Xcal, Xcalreq);
cor_vec = abs (diag (cor_vec));

% Print mean and std. deviation so we have better understanding of correlation.
mean_cor = mean (cor_vec);
std_cor = std (cor_vec);
end
