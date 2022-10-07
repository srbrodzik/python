%Name: make_gpm_cdf_file_v1a
%Purpose: Output results of interpolation to compressed netcdf4 file
%By: Stacy Brodzik

% if set to 1, limit output to 122 vertical levels (limitation of cidd/jazz)
% if set to 0, output all 176 vertical levels
limitLevels = 0;

noRainValInt = -1111;
missingVal = -99;
missingValByte = 255;
missingValInt = -9999;
title = 'Interpolated GPM 2Ku Data';
source_uw = 'Created by Mesoscale Group, University of Washington';
source_nasa = 'Original data obtained from NASA Goddard Earth Sciences http://pmm.nasa.gov';
data_location = 'http://gpm.atmos.washington.edu';

% compute start time
strTimeStart_forNcFile = [date_box_start(1:4),'-',date_box_start(5:6),'-',date_box_start(7:8),'T',...
    date_box_start(10:11),':',date_box_start(12:13),':',date_box_start(14:15),'Z'];
t0 = datenum('19700101 000000','yyyymmdd HHMMSS');
timeStart = round((groundTrackDateStart-t0)*86400);

% compute stop time
strTimeEnd_forNcFile = [date_box_end(1:4),'-',date_box_end(5:6),'-',date_box_end(7:8),'T',...
    date_box_end(10:11),':',date_box_end(12:13),':',date_box_end(14:15),'Z'];
t0 = datenum('19700101 000000','yyyymmdd HHMMSS');
timeEnd = round((groundTrackDateEnd-t0)*86400);

% compute gridAlt
gridAlt = (0:(LEVELS-1)) * RES;

%Get current time
formatOut = 'mm/dd/yyyy HH:MM:SS';
currentTime = datestr(now,formatOut);

%Open new netcdf file
ncid = netcdf.create([dir_out,fileInfo,'nc'],'NETCDF4');
  
%Define dimensions
%dimid_time = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
dimid_time = netcdf.defDim(ncid,'time',1);
if limitLevels
    dimid_alt = netcdf.defDim(ncid,'alt',LEVELS-54);  % for 122 limit
else
    dimid_alt = netcdf.defDim(ncid,'alt',LEVELS);
end
dimid_lat = netcdf.defDim(ncid,'lat',gridRows);
dimid_lon = netcdf.defDim(ncid,'lon',gridColumns);
% dimid_bounds = netcdf.defDim(ncid,'bounds',2);
  
%Define variables
varid_ti = netcdf.defVar(ncid,'time','double',dimid_time);
varid_start = netcdf.defVar(ncid,'start_time','double',dimid_time);
varid_stop = netcdf.defVar(ncid,'stop_time','double',dimid_time);
% varid_tb = netcdf.defVar(ncid,'time_bounds','double',[dimid_bounds,dimid_time]);
varid_lon = netcdf.defVar(ncid,'lon','float',dimid_lon);
varid_lat = netcdf.defVar(ncid,'lat','float',dimid_lat);
varid_alt = netcdf.defVar(ncid,'alt','float',dimid_alt);
varid_gm = netcdf.defVar(ncid,'grid_mapping','int',[]);
varid_rt = netcdf.defVar(ncid,'rain_type','float',[dimid_lon,dimid_lat,dimid_time]);
varid_rtr = netcdf.defVar(ncid,'rain_type_raw','float',[dimid_lon,dimid_lat,dimid_time]);
varid_pt = netcdf.defVar(ncid,'phase_type','float',[dimid_lon,dimid_lat,dimid_time]);
varid_ptr = netcdf.defVar(ncid,'phase_type_raw','float',[dimid_lon,dimid_lat,dimid_time]);
varid_shrn = netcdf.defVar(ncid,'shallow_rain_type','float',[dimid_lon,dimid_lat,dimid_time]);
varid_nsr = netcdf.defVar(ncid,'near_surf_rain','float',[dimid_lon,dimid_lat,dimid_time]);
varid_bbw = netcdf.defVar(ncid,'width_bb','float',[dimid_lon,dimid_lat,dimid_time]);
varid_bbh = netcdf.defVar(ncid,'height_bb','float',[dimid_lon,dimid_lat,dimid_time]);
varid_sw = netcdf.defVar(ncid,'swath','float',[dimid_lon,dimid_lat,dimid_time]);
varid_fe = netcdf.defVar(ncid,'flag_echo','float',[dimid_lon,dimid_lat,dimid_alt,dimid_time]);
varid_re = netcdf.defVar(ncid,'refl','float',[dimid_lon,dimid_lat,dimid_alt,dimid_time]);

%Define fill values
netcdf.defVarFill(ncid,varid_rt,false,missingVal);
netcdf.defVarFill(ncid,varid_rtr,false,missingVal);
netcdf.defVarFill(ncid,varid_pt,false,missingVal);
netcdf.defVarFill(ncid,varid_ptr,false,missingVal);
netcdf.defVarFill(ncid,varid_shrn,false,missingVal);
netcdf.defVarFill(ncid,varid_nsr,false,missingVal);
netcdf.defVarFill(ncid,varid_bbw,false,missingVal);
netcdf.defVarFill(ncid,varid_bbh,false,missingVal);
netcdf.defVarFill(ncid,varid_sw,false,missingVal);
netcdf.defVarFill(ncid,varid_fe,false,missingVal);
netcdf.defVarFill(ncid,varid_re,false,missingVal);

%Indicate variable compression
netcdf.defVarDeflate(ncid,varid_rt,false,true,5);
netcdf.defVarDeflate(ncid,varid_rtr,false,true,5);
netcdf.defVarDeflate(ncid,varid_pt,false,true,5);
netcdf.defVarDeflate(ncid,varid_ptr,false,true,5);
netcdf.defVarDeflate(ncid,varid_shrn,false,true,5);
netcdf.defVarDeflate(ncid,varid_nsr,false,true,5);
netcdf.defVarDeflate(ncid,varid_bbw,false,true,5);
netcdf.defVarDeflate(ncid,varid_bbh,false,true,5);
netcdf.defVarDeflate(ncid,varid_sw,false,true,5);
netcdf.defVarDeflate(ncid,varid_fe,false,true,5);
netcdf.defVarDeflate(ncid,varid_re,false,true,5);

%Define attributes
netcdf.putAtt(ncid,varid_ti,'standard_name','time');
netcdf.putAtt(ncid,varid_ti,'long_name','Data time');
netcdf.putAtt(ncid,varid_ti,'units','seconds since 1970-01-01T00:00:00Z');
netcdf.putAtt(ncid,varid_ti,'axis','T');
% netcdf.putAtt(ncid,varid_ti,'bounds','time_bounds');
netcdf.putAtt(ncid,varid_ti,'comment',strTimeStart_forNcFile);

netcdf.putAtt(ncid,varid_start,'units','seconds since 1970-01-01T00:00:00Z');
netcdf.putAtt(ncid,varid_start,'comment',strTimeStart_forNcFile);

netcdf.putAtt(ncid,varid_stop,'units','seconds since 1970-01-01T00:00:00Z');
netcdf.putAtt(ncid,varid_stop,'comment',strTimeEnd_forNcFile); 

% netcdf.putAtt(ncid,varid_tb,'units','seconds since 1970-01-01T00:00:00Z');
% netcdf.putAtt(ncid,varid_tb,'comment','time_bounds also stored the start and stop times, provided the time variable value lies within the start_time to stop_time interval'); 

netcdf.putAtt(ncid,varid_lon,'standard_name','longitude');
netcdf.putAtt(ncid,varid_lon,'long_name','longitude');
netcdf.putAtt(ncid,varid_lon,'units','degrees_east');
netcdf.putAtt(ncid,varid_lon,'axis','X');

netcdf.putAtt(ncid,varid_lat,'standard_name','latitude');
netcdf.putAtt(ncid,varid_lat,'long_name','latitude');
netcdf.putAtt(ncid,varid_lat,'units','degrees_north');
netcdf.putAtt(ncid,varid_lat,'axis','Y');

netcdf.putAtt(ncid,varid_alt,'standard_name','altitude');
netcdf.putAtt(ncid,varid_alt,'long_name','altitude MSL');
netcdf.putAtt(ncid,varid_alt,'units','km');
netcdf.putAtt(ncid,varid_alt,'positive','up');
netcdf.putAtt(ncid,varid_alt,'axis','Z');

netcdf.putAtt(ncid,varid_gm,'grid_mapping_name','latitude_longitude');

% rain type (varid_rt)
netcdf.putAtt(ncid,varid_rt,'units','none');
netcdf.putAtt(ncid,varid_rt,'long_name','Rain Type');
netcdf.putAtt(ncid,varid_rt,'stratiform',1);
netcdf.putAtt(ncid,varid_rt,'convective',2);
netcdf.putAtt(ncid,varid_rt,'other',3);
%netcdf.putAtt(ncid,varid_rt,'missing_value',missingVal);

% raw rain type (varid_rtr)
netcdf.putAtt(ncid,varid_rtr,'units','none');
netcdf.putAtt(ncid,varid_rtr,'long_name','Rain Type Raw');
netcdf.putAtt(ncid,varid_rtr,'stratiform_range','[10000000,19999999]');
netcdf.putAtt(ncid,varid_rtr,'convective_range','[20000000,29999999]');
netcdf.putAtt(ncid,varid_rtr,'other_range','[30000000,39999999]');
netcdf.putAtt(ncid,varid_rtr,'no_rain_value',noRainValInt);
%netcdf.putAtt(ncid,varid_rtr,'missing_value',missingVal);

% phase type (varid_pt)
netcdf.putAtt(ncid,varid_pt,'units','none');
netcdf.putAtt(ncid,varid_pt,'long_name','Near Surface Phase Type');
netcdf.putAtt(ncid,varid_pt,'solid',0);
netcdf.putAtt(ncid,varid_pt,'mixed',1);
netcdf.putAtt(ncid,varid_pt,'liquid',2);
%netcdf.putAtt(ncid,varid_pt,'missing_value',missingVal);

% raw phase type (varid_ptr)
netcdf.putAtt(ncid,varid_ptr,'units','none');
netcdf.putAtt(ncid,varid_ptr,'long_name','Near Surface Phase Type Raw');
netcdf.putAtt(ncid,varid_ptr,'phaseType_lt_100','tempC = phaseType - 100');
netcdf.putAtt(ncid,varid_ptr,'phaseType_gt_200','tempC = phaseType - 200');
netcdf.putAtt(ncid,varid_ptr,'phaseType_eq_100','top of bright band');
netcdf.putAtt(ncid,varid_ptr,'phaseType_eq_125','between top and peak of BB');
netcdf.putAtt(ncid,varid_ptr,'phaseType_eq_175','between peak and bottom of BB');
netcdf.putAtt(ncid,varid_ptr,'phaseType_eq_200','bottom of bright band');
%netcdf.putAtt(ncid,varid_ptr,'missing_value',missingVal);
%netcdf.putAtt(ncid,varid_ptr,'missing_value',missingValByte);
%netcdf.putAtt(ncid,varid_ptr,'_FillValue',missingVal);

% shallow rain flag (varid_shrn)
netcdf.putAtt(ncid,varid_shrn,'units','none');
netcdf.putAtt(ncid,varid_shrn,'long_name','Shallow Rain Flag');
netcdf.putAtt(ncid,varid_shrn,'no_shallow_rain',0);
netcdf.putAtt(ncid,varid_shrn,'shallow_isolated_maybe',10);
netcdf.putAtt(ncid,varid_shrn,'shallow_isolated_certain',11);
netcdf.putAtt(ncid,varid_shrn,'shallow_nonisolated_maybe',20);
netcdf.putAtt(ncid,varid_shrn,'shallow_nonisolated_certain',21);
netcdf.putAtt(ncid,varid_shrn,'no_rain_value',noRainValInt);
%netcdf.putAtt(ncid,varid_shrn,'missing_value',missingVal);

% near surface rain (varid_nsr)
netcdf.putAtt(ncid,varid_nsr,'units','mm/hr');
netcdf.putAtt(ncid,varid_nsr,'long_name','Near Surface Rain');
%netcdf.putAtt(ncid,varid_nsr,'missing_value',missingVal);

% BB width (varid_bbw)
netcdf.putAtt(ncid,varid_bbw,'units','meters');
netcdf.putAtt(ncid,varid_bbw,'long_name','Width of Bright Band');
%netcdf.putAtt(ncid,varid_bbw,'missing_value',missingVal);

% BB height (varid_bbh)
netcdf.putAtt(ncid,varid_bbh,'units','meters');
netcdf.putAtt(ncid,varid_bbh,'long_name','Height of Bright Band');
%netcdf.putAtt(ncid,varid_bbh,'missing_value',missingVal);

% swath (varid_sw)
netcdf.putAtt(ncid,varid_sw,'units','none');
netcdf.putAtt(ncid,varid_sw,'long_name','GPM-Ku coverage area');
netcdf.putAtt(ncid,varid_sw,'swath',0);
netcdf.putAtt(ncid,varid_sw,'no_swath',1);
%netcdf.putAtt(ncid,varid_sw,'missing_value',missingVal);

% flag echo (varid_fe)
netcdf.putAtt(ncid,varid_fe,'units','none');
netcdf.putAtt(ncid,varid_fe,'long_name','Flag of precip and main/side lobe clutter');
netcdf.putAtt(ncid,varid_fe,'precip',5);
netcdf.putAtt(ncid,varid_fe,'main_lobe_clutter',16);
netcdf.putAtt(ncid,varid_fe,'precip_plus_main_lobe_clutter',21);
netcdf.putAtt(ncid,varid_fe,'side_lobe_clutter',64);
netcdf.putAtt(ncid,varid_fe,'precip_plus_side_lobe_clutter',69);
%netcdf.putAtt(ncid,varid_fe,'missing_value',missingVal);

% reflectivity (varid_re)
netcdf.putAtt(ncid,varid_re,'units','dBZ');
netcdf.putAtt(ncid,varid_re,'long_name','GPM-Ku Reflectivity');
%netcdf.putAtt(ncid,varid_re,'missing_value',missingVal);

%Define global variables
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.0');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title',title);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'source_uw',source_uw);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'source_nasa',source_nasa);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'reference',namein2Ku);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'data_location',data_location);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'orbit',orbitS);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'region',long_region);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lon_min',MinLongEast);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lon_max',MaxLongEast);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lat_min',MinLatNorth);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lat_max',MaxLatNorth);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history', ...
    ['File created ',currentTime] );

%End definitions
netcdf.endDef(ncid);

%Put data in file
netcdf.putVar(ncid,varid_ti,timeStart);
netcdf.putVar(ncid,varid_start,timeStart);
netcdf.putVar(ncid,varid_stop,timeEnd);
% netcdf.putVar(ncid,varid_tb,[timeStart,timeEnd]);
gridLongitude2 = reshape(gridLongitude,gridRows,gridColumns);
netcdf.putVar(ncid,varid_lon,gridLongitude2(1,:));
netcdf.putVar(ncid,varid_lat,gridLatitude(1:gridRows));
if limitLevels
    netcdf.putVar(ncid,varid_alt,gridAlt(:,1:122));  % for 122 limit
else
    netcdf.putVar(ncid,varid_alt,gridAlt);
end
netcdf.putVar(ncid,varid_rt,(reshape(gridRainType,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_rtr,(reshape(gridRainTypeRaw,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_pt,(reshape(gridSurfPhase,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_ptr,(reshape(gridSurfPhaseRaw,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_shrn,(reshape(gridShallowRainType,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_nsr,(reshape(gridSurfRain,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_bbw,(reshape(gridWidthBB,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_bbh,(reshape(gridHeightBB,gridRows,gridColumns))');
netcdf.putVar(ncid,varid_sw,(reshape(gridNoSwath,gridRows,gridColumns))');
gridFlagEcho = reshape(gridFlagEcho,gridRows,gridColumns,LEVELS);
gridFlagEcho = permute(gridFlagEcho,[2,1,3]);
if limitLevels
    netcdf.putVar(ncid,varid_fe,gridFlagEcho(:,:,1:122));  % for 122 limit
else
    netcdf.putVar(ncid,varid_fe,gridFlagEcho);
end
gridRefl = reshape(gridRefl,gridRows,gridColumns,LEVELS);
gridRefl = permute(gridRefl,[2,1,3]);
if limitLevels
    netcdf.putVar(ncid,varid_re,gridRefl(:,:,1:122));  % for 122 limit
else
    netcdf.putVar(ncid,varid_re,gridRefl);
end

%Close file
netcdf.close(ncid);

