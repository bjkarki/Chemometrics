function T2ratio = HotellingT2ratio (Xcaluc, Xvaluc, n)

% HotellingT2ratio performs calculation of Hotelling T squared of val over mean of Hotelling T squared of cal.
%
% Syntax:
% T2ratio = HotellingT2ratio (Xcaluc, xvaluc, n)
%
% Input arguments:
% Xcaluc = Calibration Matrix in primary condition (un-mean centered)
% xvaluc = a validation sample in secondary condition (un-mean centered)
% n = no. of eigs to be used

for i = 1 : size (Xvaluc, 1)
    
    % T2 for each of the val sample
    T2val (i) = HotellingT2val (Xcaluc, Xvaluc (i, :), n);

end

%T2 for cal sample
T2cal = HotellingT2cal (Xcaluc, n);

% T2ratio
T2cal = mean (T2cal);
T2ratio = T2val ./ T2cal;

end

function T2val = HotellingT2val (Xcaluc, xvaluc, n)

% HotellingT2val performs calculation of Hotelling T squared of a validation sample
%
% Syntax:
% T2val = HotellingT2val (Xcaluc, xvaluc, n)
%
% Input arguments:
% Xcaluc = Calibration Matrix in primary condition (un-mean centered)
% xvaluc = a validation sample in secondary condition (un-mean centered)
% n = no. of eigs to be used

% Validation sample check
if size (xvaluc, 1) ~= 1
    if size (xvaluc, 2) ~= 1
        error ('Only 1 validation sample can be used at a time');
    end
end

% Mean-centering
[Xcal, xmean] = mncn (Xcaluc);
xval = scale (xvaluc, xmean);
clear xmean

% Conversion into column vector
if size (xval, 1) == 1
    xval = xval';
end

% T2 for val sample
[~, S, V] = svd (Xcal);
tval = V (:, 1 : n)' * xval;
C = S (1 : n, 1 : n).^2;
T2val = tval' / C * tval;
clear S V C tval Xcal xval


end             % ends HotellingT2val function

function T2cal = HotellingT2cal (Xcaluc, n)

% HotellingT2cal performs calculation of Hotelling T squared of leave one out cal samples
%
% Syntax:
% T2cal = HotellingT2cal (Xcaluc, n)
%
% Input arguments:
% Xcaluc = Calibration Matrix in primary condition (un-mean centered)
% n = no. of eigs to be used

for i = 1 : size (Xcaluc, 1)
    
    % required Xcaluc
    reqXcaluc = Xcaluc;
    xvalcaluc = reqXcaluc (i, :);
    reqXcaluc (i, :) = [];
    
    % mean centering
    [reqXcal, xmean] = mncn (reqXcaluc);
    xvalcal = scale (xvalcaluc, xmean);
    clear xmean reqXcaluc xvalcaluc
    
    % Conversion into column vector
    if size (xvalcal, 1) == 1
        xvalcal = xvalcal';
    end
    
    % T2 for a valcal sample
    [~, S, V] = svd (reqXcal);
    tvalcal = V (:, 1 : n)' * xvalcal;
    C = S (1 : n, 1 : n).^2;
    T2cal (i) = tvalcal' / C * tvalcal;
    clear S V C tvalcal reqXcal xvalcal

end             % ends the for loop

end             % ends the HotellingT2cal function