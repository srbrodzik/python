% Name: prepare_data
% Purpose: Set appropriate values to missingVal, get rid of NaN's, etc
% By: Stacy Brodzik

missingVal = -99;
missingValInt = -9999;
missingValByte = 255;
noRainValInt = -1111;

% gridRainType - replace missing vals = -9999 and no rain vals = -1111
gridRainType(gridRainType < 0) = missingVal;

% gridRainTypeRaw
gridRainTypeRaw(isnan(gridRainTypeRaw)) = missingVal;

% gridSurfPhase
gridSurfPhase(gridSurfPhase < 0) = missingVal;
gridSurfPhase(gridSurfPhase == missingValByte) = missingVal;
gridSurfPhase(isnan(gridSurfPhase)) = missingVal;

% gridSurfPhaseRaw
gridSurfPhaseRaw(gridSurfPhaseRaw < 0) = missingVal;
%gridSurfPhaseRaw(gridSurfPhaseRaw == missingValByte) = missingVal;
gridSurfPhaseRaw(isnan(gridSurfPhaseRaw)) = missingVal;

% gridShallowRainType - replace missing vals = -9999 and no rain vals = -1111
gridShallowRainType(gridShallowRainType < 0) = missingVal;

% gridSurfRain
gridSurfRain(gridSurfRain < 0) = missingVal;

% gridWidthBB - replace missing vals = -9999. and no rain vals = -1111.1
gridWidthBB(gridWidthBB < 0) = missingVal;

% gridHeightBB - replace missing vals = -9999. and no rain vals = -1111.1
gridHeightBB(gridHeightBB < 0) = missingVal;

% gridNoSwath
gridNoSwath(gridNoSwath <= missingVal) = missingVal;

% gridFlagEcho
gridFlagEcho(gridFlagEcho < 0) = missingVal;

% gridRefl
gridRefl(gridRefl <= missingVal) = missingVal;





