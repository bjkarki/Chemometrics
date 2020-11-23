function cumsum_plot (Xcaluc)
% cumsum_plot plots Cumulative sum plot derived from Singular Value Decomposition (SVD).
%
% Syntax:
% cumsum_plot (Xcaluc);
%
% Input argument:
% Xcaluc = Calibration Samples in primary condition (un-mean centered).
%
% Ouput plot:
% cumulative sum plot

% Mean-center Xcal
%[Xcal, ~] = mncn (Xcaluc);

% SVD
[~, S, ~] = svd (Xcaluc);

% Required variable manipulation.
S = diag (S);
Sd = cumsum (S);
Ssum = sum (S);

% cumulative sum vector
cu = Sd / Ssum;

% cumulative sum plot
figure, plot (cu);
xlabel ('# of eigs')
ylabel ('cumulative sum of Singular values')
title ('cumsum plot')
axis tight

end