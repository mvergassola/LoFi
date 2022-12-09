function [gamma] = gammaPredictor(width,h_bay,n_bays,D_legs,D_braces)
%   gammaPredictor : Predicts the gamma coefficients based on the
%   polymodels obtained via multivariate nonlinear regression analysis.
%
%   gamma : vector containing the 3 gamma coefficientes
%   width : width of the bay - [m]
%   h_bay : height of the bay - [m]
%   n_bays : number of bays - [m]
%   D_legs : diameter of the legs - [m]
%   D_braces : diameter of the braces - [m]

% Load polymodels
load('polymodel1_2ndOrder','polymodel1')
load('polymodel2_2ndOrder','polymodel2')
load('polymodel3_2ndOrder','polymodel3')

% Definition of the independent variables
indepVar1 = width;
indepVar2 = h_bay;
indepVar3 = n_bays;
indepVar4 = D_legs;
indepVar5 = D_braces;

% Check validity of the polymodels
if D_legs/D_braces < 1 || D_legs/D_braces > 2.5
    warning('The structure is outside the investigated domain. Please use these results carefully.')
elseif h_bay/width < 0.8 || h_bay/width > 2.5
    warning('The structure is outside the investigated domain. Please use these results carefully.')
elseif width > h_bay*n_bays
    warning('The structure is outside the investigated domain. Please use these results carefully.')
end

% Prediction of the dependent variables
indepvar = [indepVar1, indepVar2, indepVar3, indepVar4, indepVar5];                    
gamma1 = abs(polyvaln(polymodel1,indepvar));
gamma2 = abs(polyvaln(polymodel2,indepvar));
gamma3 = abs(polyvaln(polymodel3,indepvar));

% Definition of the gamma vector
gamma = [gamma1, gamma2, gamma3];

end
