function output = ucurve (variance, bias)
% ucurve converts lcurve into ucurve.
% In this process, the lcurve is modified into ucurve
%
% Syntax:
% ucurve (variance_merit, bias_merit)
%
% Input arguments:
% variance_merit = variance measure 
% bias_merit = bias measure
%
% Output arguments:
% ucurve vector plot obtained after function operation.

% Maximum values
biasmax = max (bias);
variancemax = max (variance);

% Minimum values
biasmin = min (bias);
variancemin = min (variance);

% loop to determine ucurve
for i = 1 : length (bias)
    
    % Assign variables
    biasi = bias (i);
    variancei = variance (i);
    
    output (i) = ((variancei-variancemin)/(variancemax-variancemin)) + ((biasi - biasmin)/(biasmax-biasmin));
    
end

% Position of minimum value in the ucurve
[~, position] = minimum (output)

% ucurve plot
figure, plot (output);

end
    
