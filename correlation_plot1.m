function correlation_plot1 (Xcaluc, xvaluc)
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
[mean_cor(i), std_cor(i)] = correlation_val1 (Xcaluc, xvaluc, i);
[mean_cor_cal(i), std_cor_cal(i)] = correlation_cal1 (Xcaluc, i);
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

function [mean_cor, std_cor] = correlation_val1 (Xcaluc, xvaluc, eig)

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

% keep track of no. of rows
[r, ~] = size (Xcaluc);

% SVD
[~, ~, V] = svd (Xcaluc);
Vreq = V (:, 1 : eig);
clear U S V;

% xvalcaluc
xvalcaluc = ((Vreq * Vreq') * xvaluc')';

% Building Xvalcal matrix
Xvalcaluc = [];
for i = 1 : r
    Xvalcaluc = [Xvalcaluc; xvalcaluc];
end

% Transpose of each matrix because corr does pairwise correlation of columns
Xcaluc = Xcaluc'; 
Xvalcaluc = Xvalcaluc'; 


% correlation vector of samples
cor_vec = corr (Xcaluc, Xvalcaluc);
cor_vec = abs (diag (cor_vec));

% Print mean and std. deviation so we have better understanding of correlation.
mean_cor = mean (cor_vec);
std_cor = std (cor_vec);

end

function [mean_cor, std_cor] = correlation_cal1 (Xcaluc, n)

% Initial verification step
if nargin > 2
    error ('Too many input arguments')
elseif nargin < 2
    error ('Not enough input arguments')
end

% SVD
[U, S, V] = svd (Xcaluc);
Ureq = U (:, 1 : n); 
Sreq = S (1 : n, 1 : n); 
Vreq = V (:, 1 : n);
clear U S V;

% Xcalrequc
Xcalrequc = Ureq * Sreq * Vreq';

% Transpose of each matrix because corr does pairwise correlation of columns
Xcaluc = Xcaluc'; 
Xcalrequc = Xcalrequc'; 


% correlation vector of samples
cor_vec = corr (Xcaluc, Xcalrequc);
cor_vec = abs (diag (cor_vec));

% Print mean and std. deviation so we have better understanding of correlation.
mean_cor = mean (cor_vec);
std_cor = std (cor_vec);
end