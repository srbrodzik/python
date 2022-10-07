function [wetbulbTemp] = wetbulb(P,Td,T)
%%wetbulb
    %Function to calculate wetbulb temperature given pressure, dewpoint, and
    %temperature.
    %
    %General form:
    % [wetbulbTemp] = wetbulb(P,Td,T)
    %Inputs:
    %P: pressure in hectopascals
    %Td: dewpoint in degrees C
    %T: temperature in degrees C
    %
    %Version Date: 6/21/2018
    %Last major revision: 11/27/2017
    %Written by: Daniel Hueholt
    %North Carolina State University
    %Undergraduate Research Assistant at Environment Analytics
    %
    %See also addWetbulb
    %

epsilon = 0.622;
Lv = 2.5*10^6; %J/kg
Cp = 1005; %J/kg
psychro = (Cp.*P)./(epsilon.*Lv); %Psychrometric constant
syms Tw %Declare the Tw symbol, required for symbolic toolbox operations
%T, P, e
eAct = 6.11*10.^((7.5.*Td)./(237.3+Td)); %Actual vapor pressure is calculated from dewpoint

wetbulbTemp = vpasolve(eAct == 6.11*10.^((7.5.*Tw)./(237.3+Tw))-psychro.*(T-Tw),Tw,[-100 50]); %Solves the wetbulb equation numerically

end