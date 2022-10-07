function [wetbulbTemp] = wetbulb_AMS(P,Td,T)
%%wetbulb_AMS
    %Function to calculate wetbulb temperature given pressure, dewpoint, and
    %temperature. Uses the psychrometric formula from the American
    %Meteorological Society glossary. An alternate formula from the
    %notes to MEA312: Atmospheric Thermodynamics at NC State is commented
    %out. Both formulae output relatively similar results.
    %Tetens's formula is used for vapor pressure calculations.
    %
    %General form:
    % [wetbulbTemp] = wetbulb(P,Td,T)
    %Inputs:
    %P: pressure in hectopascals
    %Td: dewpoint in degrees C
    %T: temperature in degrees C
    %
    %Version Date: 10/4/2019
    %Last major revision: 10/4/2019
    %Written by: Daniel Hueholt
    %North Carolina State University
    %Undergraduate Research Assistant at Environment Analytics
    %
    %See also addWetbulb
    %

% epsilon = 0.622;
% Lv = 2.5*10^6; %J/kg
% Cp = 1005; %J/kg
% psychro = (Cp.*P)./(epsilon.*Lv); %Psychrometric constant
syms Tw %Declare the Tw symbol, required for symbolic toolbox operations
%T, P, e
eAct = 6.11*10.^((7.5.*Td)./(237.3+Td)); %Actual vapor pressure is calculated from dewpoint

% wetbulbTemp = vpasolve(eAct == 6.11*10.^((7.5.*Tw)./(237.3+Tw))-psychro.*(T-Tw),Tw,[-100 50]); %Solves the wetbulb equation numerically

if T>0
    wetbulbTemp = vpasolve(eAct == 6.11*10.^((7.5.*Tw)./(237.3+Tw))-6.60.*10^(-4).*(1+0.00115.*Tw).*P.*(T-Tw),Tw,[-100 50]); %Solves the wetbulb equation numerically
else
    wetbulbTemp = vpasolve(eAct == 6.11*10.^((7.5.*Tw)./(237.3+Tw))-5.82.*10^(-4).*(1+0.00115.*Tw).*P.*(T-Tw),Tw,[-100 50]); %Solves the wetbulb equation numerically
end

end
