% Name: getGPM_Ku_interp_impacts
% Purpose: Process raw GPM data (2Ku)
% Updates: 
% By: Ulrike Romatschke (Original version by Darren Wilton)
% Updated and Modified by MD.Zuluaga (2011/09/26)
% Updated and Modified by S.Brodzik for GPM (2015/03/02 and 2016/07/11)
% Updated and Modified by S.Brodzik for GPM subsetted data (2017/02/28)

% Input: 2Ku GPM products
% The input files are subsetted versions of the original HDF5 files that are
% supplied via a standing order with storm.pps.eosdis.nasa.gov/storm and
% contain only the fields of interest necessary for this work:
% NS/CSF/flagShallowRain
% NS/CSF/heightBB
% NS/CSF/typePrecip
% NS/CSF/widthBB
% NS/FLG/flagEcho
% NS/Latitude
% NS/Longitude
% NS/SLV/phaseNearSurface
% NS/SLV/precipRateNearSurface
% NS/SLV/zFactorCorrected
% NS/ScanTime

% These notes refer to old version numbering system:
% Modified v74 03/31/15 Output 'shallowRainTypes', 'widthBB', 'heightBB'
% Modified v75 07/11/16 Calls C routine to interpolate all fields using
%   either weighted or nearest neighbor methods.  Output only one netcdf file 
%   with all data
% Modified v76 02/28/17 - no more need to use coordinate files since input
%   files are in HDF5 format and are not compressed

% These notes refer too new version numbering system:
% Modified v2 05/04/2017 - added log file and two new checks: one to flag
%   cases where stop time is earlier than start time (lines 669-673) and 
%   one to bail out if indices are bad (lines 616-621); this does not happen
%   often but is probably due to a bug in the code
% Modified v3 05/18/2017 - changed log file from yearly to monthly; changed
%   outfile prefix to reflect change from v04A to v05A for 2Ku input data

% Output: 1 compressed netcdf4 file per case

% Other code/functions needed: datestamper.m dateSorterGPM.m datestring.m
%    interp_for_matlab.c interp_for_matlab.h prepare_data.m
%    make_gpm_cdf_file_v1a.m

% This code runs for 176 levels and takes 6 digit orbits as input. It
% takes whole orbits as input as well as pre cut orbits.

% This code cuts out the lat long boxes, corrects the geolocation,
% interpolates the data and outputs it all in a single netcdf file.

% This code does not work for processing whole orbits!!! You must specify
% the region of interest and this cant be the whole globe!

addpath /home/disk/bob/gpm/code/
clear variables;
clear global;
close all;

t=clock;
datestamper(t);  % this Prints the Actual Time
clear t ans

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!!!!!!!!!!!!!!!If your area extends over the 180 meridian set over180=1. Default is 0. The
% negative value is MinLongEast, the positive is MaxLongEast. The resulting
% output longitude will be positive for all values, e.g. -175 becomes +185.!!!!!!!!!!!!!
over180=0; %default is zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
minus40=0; %default is zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Number of levels in GPM-2Ku data
LEVELS = 176;
RES = 0.125;      % bin to bin resolution in km

% (Does not work for the whole globe! I.e., MinLongEast can never be equal to MaxLongEast.)
MinLongEast =  -105.0;  %  coordinates of whole area of interest: N1 (Lynns Ocean Region)
MaxLongEast =   -65.0;  %  (slightly larger area of interest)
MinLatNorth  =   25.0;  %
MaxLatNorth  =   50.0;  % 
region      = 'IMP';  % IMPACTS domain
long_region = 'IMPACTS Domain';
% Do not set your southern border between -60 and -65. You can set the border at -65.

% With changing boxsize you can exclude information that produces boxes smaller than a
% specific area from being produced. The unit of boxsize is degree^2. Files
% whos surrounding lat long box is smaller or equal to boxsize are not created. 
% If you set boxsize to zero (default), all files are created whos horizontal
% projection has two or more contiguous pixels.
boxsize=0;

% Input directories for GPM 2Ku
main_dir='/home/disk/bob/gpm/';
in_dir=[main_dir,'hdf_subsetted/2Ku/'];
mex_interp_path = '/home/disk/shear2/brodzik/matlab/mexFiles';

% Output directory and prefix for netcdf file
out_dir='/home/disk/bob/gpm/nam_imp_ku/classify/ex_data_v06/';
%out_dir='/home/disk/bob/gpm/n1_oly_lynn_ocean_ku/classify/ex_data/SMALLEST_DOMAIN/';
outFilePrefix = 'GPM2Ku6_uw3';

% Output directory for log file
log_dir = out_dir;

% Radius of influence (in km) and missing value for interpolation
RADIUS = 4.25;
NODATA = -99;

% Minimum dbz for reflectivity
MINDBZ = 0.0;

% The 2Ku data should be the subsetted HDF5 files (2A-UW.GPM.Ku.*.HDF5)
% and should be located in subfolders with years and months as substructure

% Years to process -- if specific months only, enter info at lines 160 & 180
% src_years=[2014,2015,2016,2017];
src_years=[2020];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create interp_for_matlab executable; Only need to do this once and make
% sure the result (interp_for_matlab.mexa64) is in MATLABPATH

% mex interp_for_matlab.c
% movefile('interp_for_matlab.mexa64',mex_interp_path);
path(mex_interp_path,path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MinLatNorth>-67 && MinLatNorth<-60
    error('Do not set the southern border between -60 and -67 but at -67.')
end

GRIDBOXSTEP = 0.05;
RAYS = 49;
if mod(RAYS,2) == 0
    error('RAYS must be an odd number in this code.')
end
SCAN_ANGLE = 17.0;   % +/-17 deg crosstrack
angleRANGE = (-SCAN_ANGLE:(SCAN_ANGLE/fix(RAYS/2)):SCAN_ANGLE); 
MinLongEastOrig=MinLongEast;
MaxLongEastOrig=MaxLongEast;

[dummy, n_years]=size(src_years);
for yyy=1:n_years
    y_file=src_years(1,yyy);   

    % here I select the months
    if y_file==2020
        month_ini=02;
        month_fin=02;
    else
        month_ini=02;
        month_fin=02;
    end
    
    for mmm=month_ini:month_fin
        if mmm<10
            m_month=['0',num2str(mmm)];
        else
            m_month=num2str(mmm);
        end

        % Open log file
        logID = fopen([log_dir,'log_gpm.',num2str(y_file),m_month],'w');

        % List of files in 2Ku and coordinates
        fileKu=dir([in_dir,num2str(y_file),'/',m_month,'/*.HDF5']);  
        files2Ku=char(fileKu(:).name);
        [n_files, lght]=size(files2Ku);
        
        %if y_file==2014 && mmm==03
        %    initial_file=1;
        %else
            initial_file=1;
        %end    

        for filenumber=initial_file:n_files
        %for filenumber=41:41     % for orbit 010019 12/2015
            clear correctZ* flagEcho*
            clear geolocation* orbit*
            clear rainTypes* shallowRainTypes*
            clear nearSurfRain* nearSurfPhase*
            clear widthBB* heightBB*
            clear scanTime*
            clear dayFile* month* year*
            clear firstTime* lastTime*
  
            namein2Ku = char(files2Ku(filenumber,:));
            fprintf('\n')
            disp(['n=',num2str(filenumber),'/',num2str(n_files),'  file analyze ',namein2Ku]);
  
            %get orbit number and convert to string
            orbit0 = str2double(namein2Ku(51:56));
            orbitS0 = namein2Ku(51:56);

            %get date info from file name -- assume filename of this type:
            %2A-UW.GPM.Ku.V6-20160118.20140528-S195703-E212935.001403.V04A.HDF5
            dayFile0 = str2double(namein2Ku(32:33));
            month0 = str2double(namein2Ku(30:31));
            year0 = str2double(namein2Ku(26:29));
  
            %*****************************************************************
            % Read the lats and lons first - only extract data within range
            %*****************************************************************
            latitude = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/Latitude');
            longitude = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/Longitude');
            
            % FOR DEBUGGING
            % figure(1); clf; 
            % plot(longitude,latitude)

            losque=find(latitude>=MinLatNorth   & latitude<=MaxLatNorth & ...
                longitude>=MinLongEast & longitude<=MaxLongEast);
      
            [ll1, ll2]=size(losque);   %dimensions of array
            em=0;      % this will tell us there are orbits within the area
            if ll1<2   % changed 12/01/2014 from ==0 to <2
                disp('Orbit is not in area.')
                em=1;  % this will tell us there are no orbits within the area
            end
            clear losque ll1 ll2
      
            if em==0   % there are orbits in the area of interest
                
                %*************************************************
                % Load the remainder of the variables of interest.
                %*************************************************
 
                [x1, x2]=size(latitude);
                geolocation0 = single(zeros(2,x1,x2));
                geolocation0(1,:,:) = latitude(:,:);
                geolocation0(2,:,:) = longitude(:,:);
                clear x1 latitude longitude

                % FOR DEBUGGING
                % figure(1); clf; 
                % lons=reshape(geolocation0(2,:,:),1,[]);
                % lats=reshape(geolocation0(1,:,:),1,[]);
                % plot(lons,lats)
                % hold on

                % this is substituted in place of original 'findcuts' function
                [~, x2, x3]=size(geolocation0);   %dimensions of arrayexit
                wholeLongs0=zeros(x2,x3);
                wholeLongs0(:,:)=geolocation0(2,:,:);   
                wholeLats0=zeros(x2,x3);
                wholeLats0(:,:)=geolocation0(1,:,:);
                longvector=wholeLongs0;
                longvector(longvector<MinLongEast)=nan;
                longvector(longvector>MaxLongEast)=nan;
                longvector=max(longvector); % vector with longitudes of orbits in study area
                latvector=wholeLats0;
                latvector(latvector<MinLatNorth)=nan;
                latvector(latvector>MaxLatNorth)=nan;
                latvector=max(latvector);   % vector with latitudes of orbits in  study area
                longvector(isnan(latvector))=nan;

                yr = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/ScanTime/Year');
                yr=single(yr);
                mo = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/ScanTime/Month');
                mo=single(mo);
                da = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/ScanTime/DayOfMonth');
                da=single(da);
                hour0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/ScanTime/Hour');
                hour0=single(hour0);
                minute0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/ScanTime/Minute');
                minute0 = single(minute0);
                second0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/ScanTime/Second');
                second0 = single(second0);
      
                [xx2, xx3]=size(hour0);
                firstTime0 = hour0(1)*3600.+minute0(1)*60+second0(1);
                lastTime0 = hour0(xx2)*3600+minute0(xx2)*60+second0(xx2);
	
                time0=int16(transpose([yr,fix(mo),fix(da),fix(hour0),fix(minute0),fix(second0)]));
                scanTimeAll0 = double(time0(4,:))*3600. + double(time0(5,:))*60. + double(time0(6,:));
                clear yr mo da
                
                correctZ0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/SLV/zFactorCorrected');
                correctZ_ScaleFactor = 1;
	            nearSurfRain0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/SLV/precipRateNearSurface');
                nearSurfPhase0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/SLV/phaseNearSurface');
                rainTypes0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/CSF/typePrecip');
                shallowRainTypes0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/CSF/flagShallowRain');
                widthBB0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/CSF/widthBB');
                heightBB0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/CSF/heightBB');
	            flagEcho0 = h5read([in_dir,num2str(y_file),'/',m_month,'/',namein2Ku],'/NS/FLG/flagEcho');
                
                hdfml('closeall')

                %***************************************
                % selected only data in the study region
                %***************************************
                geolocation0=geolocation0(:,:,isfinite(longvector));  %trim the geolocation
                correctZ0=correctZ0(:,:,isfinite(longvector));
                nearSurfRain0=nearSurfRain0(:,isfinite(longvector));
                nearSurfPhase0=nearSurfPhase0(:,isfinite(longvector));
                hour0=hour0(isfinite(longvector),1);
                minute0=minute0(isfinite(longvector),1);
                second0=second0(isfinite(longvector),1);
                scanTimeAll0=scanTimeAll0(:,isfinite(longvector));
                rainTypes0=rainTypes0(:,isfinite(longvector));
                shallowRainTypes0=shallowRainTypes0(:,isfinite(longvector));
                widthBB0=widthBB0(:,isfinite(longvector));
                heightBB0=heightBB0(:,isfinite(longvector));
                flagEcho0=flagEcho0(:,:,isfinite(longvector));
                
                % FOR DEBUGGING
                % lons=reshape(geolocation0(2,:,:),1,[]);
                % lats=reshape(geolocation0(1,:,:),1,[]);
                % plot(lons,lats,'r')

                clear filename
                clear longvector latvector wholeLongs0 wholeLats0
                clear x1 x2 x3 xx2 xx3
	
                %**************************************************
                % Check for holes in the orbit (will happen at very
                % southern extent of orbit if at all)
                %**************************************************
                [x1, x2, x3]=size(geolocation0);
                wholeLongs0=zeros(x2,x3);
                wholeLongs0(:,:)=geolocation0(2,:,:);
                wholeLats0=zeros(x2,x3);
                wholeLats0(:,:)=geolocation0(1,:,:);
    
                long4=wholeLongs0(1,:);
                long5=zeros(1,length(long4)+1);
                long5(:,2:length(long5))=long4;
                long5(1,1)=long4(1,1);
                long4(1,length(long4)+1)=long4(1,length(long4));
	
                diffTest=abs(long4-long5);
                maximumTest=max(diffTest);
                
                if maximumTest>0.15 && sum(sum(wholeLats0 < 0)) > 0  % i.e. break in orbit in southern hemisphere
                    
                    flag_hole=1;
                    donde = find(diffTest>=0.15);  % locates longitude range of the whole area of interest in the swath
                    orbit=num2str(orbit0);
                    orbitS=orbitS0;
                    dayFile=dayFile0;
                    month=month0;
                    year=year0;
                    geolocation=cat(3,geolocation0(:,:,donde:end),geolocation0(:,:,1:donde-1));
                    correctZ=cat(3,correctZ0(:,:,donde:end),correctZ0(:,:,1:donde-1));
                    nearSurfRain=cat(2,nearSurfRain0(:,donde:end),nearSurfRain0(:,1:donde-1));
                    nearSurfPhase=cat(2,nearSurfPhase0(:,donde:end),nearSurfPhase0(:,1:donde-1));
                    scanTimeAll=cat(2,scanTimeAll0(1,donde:end),scanTimeAll0(1,1:donde-1));
                    firstTime=firstTime0;
                    lastTime=lastTime0;
                    rainTypes=cat(2,rainTypes0(:,donde:end),rainTypes0(:,1:donde-1));
                    shallowRainTypes=cat(2,shallowRainTypes0(:,donde:end),shallowRainTypes0(:,1:donde-1));
                    widthBB=cat(2,widthBB0(:,donde:end),widthBB0(:,1:donde-1));
                    heightBB=cat(2,heightBB0(:,donde:end),heightBB0(:,1:donde-1));
                    flagEcho=cat(3,flagEcho0(:,:,donde:end),flagEcho0(:,:,1:donde-1));
                    
                    % FOR DEBUGGING
                    % lons1=reshape(geolocation(2,:,:),1,[]);
                    % lats1=reshape(geolocation(1,:,:),1,[]);
                    % plot(lons1,lats1,'r')
                    % fid55=fopen([out_dir,'test/orbits_with_holes.txt'],'a');
                    % fprintf(fid55,'Hole in orbit %8.1f \n',orbit0);
                    % fclose(fid55);

                else   % no holes in orbit
                    
                    flag_hole=0;
                    orbit=num2str(orbit0);
                    orbitS=orbitS0;
                    dayFile=dayFile0;
                    month=month0;
                    year=year0;
                    geolocation=geolocation0;
                    correctZ=correctZ0;
                    nearSurfRain=nearSurfRain0;
                    nearSurfPhase=nearSurfPhase0;
                    scanTimeAll=scanTimeAll0;
                    firstTime=firstTime0;
                    lastTime=lastTime0;
                    rainTypes=rainTypes0;
                    shallowRainTypes=shallowRainTypes0;
                    widthBB=widthBB0;
                    heightBB=heightBB0;
                    flagEcho=flagEcho0;
                    
                end
                
                clear orbit0 orbitS0 dayFile0 month0 year0
                clear geolocation0 correctZ0 nearSurfRain0 nearSurfPhase0 scanTimeAll0
                clear firstTime0 lastTime0
                clear rainTypes0 shallowRainTypes0 widthBB0 heightBB0 flagEcho0
                clear directory

                %************************************
                % Create latitude & longitude vectors
                %************************************
                [a1, a2, a3]=size(geolocation);
                SCANS = a3;
                clear a1 a2 a3
                % clear latitudes longitudes
	
                latitudes(1:RAYS,1:SCANS)= geolocation(1,1:RAYS,1:SCANS);
                longitudes(1:RAYS,1:SCANS)= geolocation(2,1:RAYS,1:SCANS);
                
                % FOR DEBUGGING
                % Plot the orbit
                % plot(reshape(longitudes,1,[]),reshape(latitudes,1,[]),'m')
	
                clear geolocation SCANS
                
                mm = longitudes >= MinLongEast & longitudes <= MaxLongEast;
                nn = find(sum(mm)>0);  % locates longitude range of the whole area of interest in the swath
                NN = length(nn);
                lat = double(latitudes(:,nn));
                long1 = double(longitudes(:,nn));
                
                clear mm
                
                % FOR DEBUGGING
                % plot(reshape(long1,1,[]),reshape(lat,1,[]),'g')
                % hold on
                
                %******************************
                % Search for data in the swath.
                %******************************
                refl=correctZ(:,:,nn);
                reflpos=zeros(LEVELS,RAYS,NN);
                reflpos(refl>0)=1;    %assign 1 where Zfact > 0 avoid missing data
	
                clear nn NN refl
	
                reflpos(:,lat<MinLatNorth)=0;
                reflpos(:,lat>MaxLatNorth)=0;
	
                emptyscans=sum(sum(reflpos));  %plot(squeeze(emptyscans))
                outin=zeros(size(emptyscans)); %plot(squeeze(outin))
	
                if emptyscans==0   % start to analyze scans if 0 jump to next orbit   
                    clear reflpos
                    clear emptyscans
                    clear long1 outin
                else
                    % Assign a variable called outin with type of data!!!
                    % outin = 0 for empty scans
                    % outin = 1 for beginning of contiguous data scans
                    % outin = 2 for middle data scans
                    % outin = 3 for end data scans
                    % outin = 4 for non contiguous data scans
	  
                    % Assign values to first and last elements in outin array
                    if emptyscans(:,:,1)>0 && emptyscans(:,:,2)>0 %first scan
                        outin(1)=1;
                    elseif emptyscans(:,:,1)>0 && emptyscans(:,:,2)==0
                        outin(1)=4;
                    end
	  
                    if emptyscans(:,:,length(emptyscans)-1)>0 && emptyscans(:,:,length(emptyscans))>0 %last scan
                        outin(length(emptyscans))=3;
                    elseif emptyscans(:,:,length(emptyscans)-1)==0 && emptyscans(:,:,length(emptyscans))>0
                        outin(length(emptyscans))=4;
                    end
	  
                    % assign values to outin between first and last scans
                    for i=2:length(emptyscans)-1
                        if emptyscans(:,:,i)>0
                            if emptyscans(:,:,i-1)>0 && emptyscans(:,:,i+1)>0
                                outin(i)=2;
                            elseif emptyscans(:,:,i-1)>0 && emptyscans(:,:,i+1)==0
                                outin(i)=3;
                            elseif emptyscans(:,:,i-1)==0 && emptyscans(:,:,i+1)>0
                                outin(i)=1;
                            elseif emptyscans(:,:,i-1)==0 && emptyscans(:,:,i+1)==0
                                outin(i)=4;
                            end
                        end
                    end
                    
                    % FOR DEBUGGING
                    % figure(3); clf; 
                    % cos1=find(outin==1); [d1 d2]=size(cos1);
                    % plot(reshape(long1(:,cos1),1,RAYS*d1),reshape(lat(:,cos1),1,RAYS*d1),'g')
                    % hold on
                    % cos1=find(outin==2); [d1 d2]=size(cos1);
                    % plot(reshape(long1(:,cos1),1,RAYS*d1),reshape(lat(:,cos1),1,RAYS*d1),'r')
                    % cos1=find(outin==3); [d1 d2]=size(cos1);
                    % plot(reshape(long1(:,cos1),1,RAYS*d1),reshape(lat(:,cos1),1,RAYS*d1),'b')
                    % cos1=find(outin==4); [d1 d2]=size(cos1);
                    % plot(reshape(long1(:,cos1),1,RAYS*d1),reshape(lat(:,cos1),1,RAYS*d1),'black')
	
                    %************************************************
                    % Calculate the lat and long boxes where there is
                    % data and cut these boxes out of the swath.
                    %************************************************
                    reflpos2D=sum(reflpos);         % accounts only for horizontal 2D       
                    contiguous=zeros(size(outin));  % var to account for single contiguous scans within swath    
	  
                    clear reflpos emptyscans
	  
                    indices=find(outin==4); %non contiguous scans with just a single full pixel are ignored
                    [e1, e2]=size(indices);
                    for i=1:e1
                        a=indices(i);
                        for j=1:RAYS-1
                            if reflpos2D(:,j,a)>0 && reflpos2D(:,j+1,a)>0
                                contiguous(:,:,a)=1;
                            end
                        end
                    end
	  
                    clear reflpos2D indices e1 e2 a
	  
                    %**************************************************************
                    % Localize the outerbounds of the regions where data is located
                    % and search for the longitude limits
                    %**************************************************************
                    longminimum=zeros(size(outin));
                    longminimum(:,:,:)=nan;
                    longminimum(1,1,outin==1)=min(long1(:,outin==1));
                    longminimum(1,1,contiguous==1)=min(long1(:,contiguous==1));
                    longminimum=longminimum(isfinite(longminimum));
	
                    longmaximum=zeros(size(outin));
                    longmaximum(:,:,:)=nan;
                    longmaximum(1,1,outin==3)=max(long1(:,outin==3));
                    longmaximum(1,1,contiguous==1)=max(long1(:,contiguous==1));
                    longmaximum=longmaximum(isfinite(longmaximum));
	  
                    longlimits2=cat(1,longminimum,longmaximum);
                    [q1, q2, q3]=size(longlimits2);
                    clear longlimits
                    longlimits=zeros(2,q3);
                    longlimits(:,:)=longlimits2;
	  
                    clear long1 outin contiguous
                    clear longminimum longmaximum longlimits2
	
                    %************************
                    %put small boxes together
                    %************************
                    if q3>1   %so What if q3 is less than 1????
                        for i=1:q3
                            if longlimits(1,i)>=longlimits(2,i)
                                [d1, d2]=size(longlimits);
                                fprintf('\n')
                                disp('There is a reverse maximum in swath')
                                while longlimits(2,i+1)<longlimits(1,i)
                                    longlimits(1,i)=longlimits(1,i+1);
                                    longlimits(:,i+1)=[];
                                    [d1, d2]=size(longlimits);
                                    if i+1>d2
                                        break
                                    end
                                end
                            else
                                [d1, d2]=size(longlimits);
                                while longlimits(1,i+1)<longlimits(2,i)
                                    longlimits(2,i)=longlimits(2,i+1);
                                    longlimits(:,i+1)=[];
                                    [d1, d2]=size(longlimits);
                                    if i+1>d2
                                        break
                                    end
                                end
                            end
                            if i+2>d2
                                break
                            end
                        end
                    end
                    clear q1 q2 q3 d1 d2
	  
                    % longlimits has now the longitude limites of the boxes.
                    
                    %***************************************************************************
                    % Go through all boxes listed in longlimits and calculate all needed output.
                    %***************************************************************************
                    clear f1 f2
                    [f1, f2]=size(longlimits);
	  
                    fprintf('\n')
                    disp(['Number of boxes in swath to be analyzed: ',num2str(f2)])
	  
                    if f2>0        %for more than 1 region!!!
                        for inoneorbit=1:f2
                        %for inoneorbit=2:f2  % for testing only
                            MINLONGEAST=longlimits(1,inoneorbit);
                            MAXLONGEAST=longlimits(2,inoneorbit);
	      
                            if MINLONGEAST<MinLongEast
                                MINLONGEAST=MinLongEast;
                            end
	      
                            if MAXLONGEAST>MaxLongEast
                                MAXLONGEAST=MaxLongEast;
                            end
                            
                            % SRB - possible fix for problem I had with orbit 332
                            % if MAXLONGEAST<(MinLongEast+0.06) || MINLONGEAST>(MaxLongEast-0.06)
	      
                            if MAXLONGEAST<(MinLongEast+0.05) || MINLONGEAST>(MaxLongEast-0.05)
                                fprintf('\n')
                                disp('box found in a border of the region...skiping.this.box...')
                                clear MINLONGEAST
                                clear MAXLONGEAST
                            else
                                m = longitudes >= MINLONGEAST & longitudes <= MAXLONGEAST;
                                n = find(sum(m)>0);  % locates appropriate longitude range
                                N = length(n);
                                N1 = n(1,1);  % index of the first good value in scans
                                Nn = n(1,N);  % index of the last good value in scans
                                
                                % check to make sure indices make sense
                                if N ~= (Nn-N1)+1
                                    fprintf(logID,'%s %s%s%s %s %s\n','problems with indices...skipping box', ...
                                        num2str(inoneorbit),'/',num2str(f2),'for',namein2Ku);
                                    clear MINLONGEAST
                                    clear MAXLONGEAST
                                    clear m n N N1 Nn
                                else
                                    % FOR DEBUGGING
                                    % figure(4); clf;
                                    % cosM=find(m==1);
                                    % plot(longitudes(cosM),latitudes(cosM))
                                    % hold on
                                    % lons2=longitudes(:,n); lats2=latitudes(:,n);
                                    % plot(reshape(lons2,1,RAYS*N),reshape(lats2,1,RAYS*N),'r')

                                    clear m

                                    centreScan = ceil((Nn-N1)/2);
                                    if Nn==N1
                                        centreScan=1;
                                    end

                                    %************************************
                                    % Select data within the analyzed box
                                    %************************************
                                    lat = double(latitudes(:,n));   
                                    long = double(longitudes(:,n));
                                    sR1 = double(nearSurfRain(:,n));     %nearSurfRain
                                    sP1 = double(nearSurfPhase(:,n));    %nearSurfPhase
                                    rt1 = double(rainTypes(:,n));        %rainTypes
                                    srt1 = double(shallowRainTypes(:,n));%shallowRainTypes
                                    wBB1 = double(widthBB(:,n));         %widthBB
                                    hBB1 = double(heightBB(:,n));        %heightBB
                                    
                                    %*******************************
                                    % Select the scantime of the box
                                    %*******************************
                                    scanTime = scanTimeAll(N1:Nn); 
                                    
                                    %******************************************
                                    % Calculate the correct time and date (UTC)
                                    %******************************************
                                    groundTrackDateStart=dateSorterGPM(1,dayFile,month,year,firstTime,lastTime,scanTime);
                                    groundTrackDateEnd=dateSorterGPM(length(scanTime),dayFile,month,year,firstTime,lastTime,scanTime);

                                    %******************************************************
                                    % Calculate the date and time in YYYYMMDD.hhmmss format
                                    %******************************************************
                                    date_box_start=datestring(groundTrackDateStart,1);
                                    date_box_end=datestring(groundTrackDateEnd,1);

                                    clear centreScan scanTime

                                    if groundTrackDateStart > groundTrackDateEnd
                                        fprintf(logID,'%s %s%s%s %s %s\n','problem with times for box', ...
                                            num2str(inoneorbit),'/',num2str(f2),'for',namein2Ku);
                                        fprintf(logID,'   start=%s and end=%s\n',date_box_start,date_box_end);
                                    end
                                        
                                    %*****************************************
                                    % Assign the file name of the analyzed box
                                    %*****************************************
                                    fileInfo = [outFilePrefix,'_',date_box_start,'_to_',date_box_end,'_',orbitS,'_',region,'.'];

                                    fprintf('\n')
                                    disp(['the analyzed box will be named ',fileInfo])

                                    %*************************************************************************
                                    % Create arrays of the radar reflectivity and flag echo at all height bins 
                                    % (Z### and FE###).  This does not create arrays; it creates LEVELS 
                                    % variables with Z and FE at each level
                                    %*************************************************************************
                                    for k=1:LEVELS  % LEVELS = 176 for REAL files
                                        eval(['Z',num2str(k),'(1:RAYS,1:N)=correctZ(' num2str(k),',1:RAYS,N1:Nn);']);
                                        eval(['Z',num2str(k),' = double(Z' num2str(k),');']);
                                        eval(['FE',num2str(k),'(1:RAYS,1:N)=flagEcho(' num2str(k),',1:RAYS,N1:Nn);']);
                                        eval(['FE',num2str(k),' = double(FE' num2str(k),');']);
                                    end
                                    clear N1 Nn

                                    %*********************************************************************
                                    % Correct geolocation and give latitude/longitude coordinates based on
                                    %  scan & height.  NOTE: Level 176 is the earth ellipsoid
                                    %*********************************************************************
                                    heightmatrix=zeros(RAYS,N,LEVELS);      % This is for Reflectivity
                                    for i=1:LEVELS
                                        heightmatrix(:,:,i)=(LEVELS-i)*RES; % height in km   from 22km to 0
                                    end

                                    angle2D=zeros(RAYS,N);
                                    for i=1:N
                                        angle2D(:,i)=angleRANGE;   % Assign the angle of each scan in 2D    
                                    end

                                    anglerangematrix=zeros(size(heightmatrix));
                                    for i=1:LEVELS
                                        anglerangematrix(:,:,i)=angle2D;
                                    end

                                    %******************************************
                                    % Matrix with shift distance for each level
                                    %******************************************
                                    kmshift=heightmatrix.*tan(deg2rad(anglerangematrix));  %%in km 
                                    range=abs(km2deg(kmshift));                            %% in deg  

                                    %*****************************************************************
                                    % Calculate the angle in which the geolocation needs to be shifted
                                    %*****************************************************************
                                    az = azimuth(lat(1,1:N),long(1,1:N),lat(RAYS,1:N),long(RAYS,1:N));
                                    directionshiftBOTTOM=az;
                                    if az < 180
                                        directionshiftTOP = az + 180;
                                    else
                                        directionshiftTOP = az - 180;
                                    end

                                    lat3D=zeros(size(heightmatrix));   %these are latitudes in 3D matrix   
                                    for i=1:LEVELS
                                        lat3D(:,:,i)=lat;
                                    end

                                    long3D=zeros(size(heightmatrix));  %these are longitudes in 3D matrix 
                                    for i=1:LEVELS
                                        long3D(:,:,i)=long;
                                    end

                                    clear heightmatrix angle2D anglerangematrix kmshift az

                                    %****************************************************
                                    % This is a matrix with 24 values of N-scans azimuths
                                    %****************************************************
                                    directionshiftBOTTOM2D=zeros(fix(RAYS/2),N);
                                    for i=1:fix(RAYS/2)
                                        directionshiftBOTTOM2D(i,:)=directionshiftBOTTOM;
                                    end

                                    directionshiftTOP2D=zeros(fix(RAYS/2+1),N);
                                    for i=1:fix(RAYS/2+1)
                                        directionshiftTOP2D(i,:)=directionshiftTOP;
                                    end

                                    directionshiftBOTTOM3D=zeros(fix(RAYS/2),N,LEVELS);
                                    for i=1:LEVELS
                                        directionshiftBOTTOM3D(:,:,i)=directionshiftBOTTOM2D;
                                    end

                                    directionshiftTOP3D=zeros(fix(RAYS/2+1),N,LEVELS);
                                    for i=1:LEVELS
                                        directionshiftTOP3D(:,:,i)=directionshiftTOP2D;
                                    end

                                    %******************************************************
                                    % Lats and lons corrected for Reflectivity and flagEcho
                                    %******************************************************
                                    [latsB,longsB] = reckon(lat3D(1:fix(RAYS/2),1:N,1:LEVELS),long3D(1:fix(RAYS/2),1:N,1:LEVELS),...
                                        range(1:fix(RAYS/2),1:N,1:LEVELS),directionshiftBOTTOM3D(:,:,1:LEVELS));

                                    [latsT,longsT] = reckon(lat3D(fix(RAYS/2+1):RAYS,1:N,1:LEVELS),long3D(fix(RAYS/2+1):RAYS,1:N,1:LEVELS),...
                                        range(fix(RAYS/2+1):RAYS,1:N,1:LEVELS),directionshiftTOP3D(:,:,1:LEVELS));

                                    clear range directionshift* lat3D long3D

                                    lats=cat(1,latsB,latsT);     %corrected latitude
                                    longs=cat(1,longsB,longsT);  %corrected longitude

                                    % FOR DEBUGGING
                                    % plot reflectivity at every level (earth ellipsoid at Z176)
                                    % name = 'gpm.20140705.205454.';
                                    % name = 'gpm.20140706.200809.';
                                    % for z=169:176
                                    %      zstr=num2str(z);
                                    %      km = (176 - z) * 0.125;
                                    %      eval(['figure(',zstr,')']); clf;
                                    %      eval(['R',zstr,' = Z',zstr,';']);
                                    %      eval(['ind = find(R',zstr,'==-9999);']);
                                    %      eval(['R',zstr,'(ind)=NaN;']);
                                    %      eval(['R',zstr,' = 10*log10(R',zstr,');']);
                                    %      eval(['contourf(longs(:,:,',zstr,'),lats(:,:,',zstr,'),R',zstr,')']);
                                    %      colorbar
                                    %      %caxis([30,38])
                                    %      xlabel('longitude')
                                    %      ylabel('latitude')
                                    %      title(['Reflectivity at ',num2str(km),' km'])
                                    %      %axis([-81,-77,29.5,34.5])
                                    %      %axis([-56,-50,-42,-35])   % 20140705
                                    %      axis([-43,-35,-48,-42])   % 20140706
                                    %      colormap jet
                                    %      %eval(['minR = nanmin(nanmin(R',zstr,'));']);
                                    %      eval(['maxR = nanmax(nanmax(R',zstr,'));']);
                                    %      %str={['min = ',num2str(minR],['max = ',num2str(maxR)]};
                                    %      str={['max = ',num2str(maxR)]};
                                    %      %annotation('textbox',[0.15,0.75,0.1,0.1],'String',str);
                                    %      annotation('textbox',[0.5,0.75,0.1,0.1],'String',str);
                                    %      if length(zstr)==1
                                    %          zstr=['00',zstr];
                                    %      elseif length(zstr)==2
                                    %          zstr=['0',zstr];
                                    %      end
                                    %      eval(['print -dpng ',name,zstr,'.png']);
                                    % end

                                    clear latsB latsT longsB longsT

                                    % FOR DEBUGGING
                                    % whos('lats','longs','lats_LH','longs_LH')	      
                                    % figure(5); clf;   %plot the corrected coordinates
                                    % plot(long(:,1),lat(:,1),'--rs')                   %these are the original coordinates
                                    % hold on
                                    % plot([long(1,1),long(RAYS,1)],[lat(1,1),lat(RAYS,1)],'-ks') %straigt line between edges   
                                    % plot(longs(:,1,1),lats(:,1,1),'--gs')             %these are the corrected 21.875km coord
                                    % plot(longs(:,1,LEVELS),lats(:,1,LEVELS),'--s')    %these are the corrected 0.00km coord
                                    % plot(long(:,1),lat(:,1),'--rs')                   %these are the original coordinates
                                    % plot([long(1,1),long(RAYS,1)],[lat(1,1),lat(RAYS,1)],'-ks') %straight line between edges   
                                    % plot(longs_LH(:,1,1),lats_LH(:,1,1),'--gs')       %these are the corrected 0.00km coord
                                    % plot(longs_LH(:,1,19),lats_LH(:,1,19),'--s')      %these are the corrected 17km coord

                                    %*****************************************************
                                    % The product of this section of code is the corrected 
                                    % latitude and longitude called lats## and longs##
                                    %*****************************************************

                                    %***********************************************************
                                    % The following code calculates the dimensions of the box in
                                    % the new cartesian grid on which the reflectivity will be
                                    % interpolated. Gridspacing is 0.05 degrees.
                                    %***********************************************************

                                    reflsmall=correctZ(:,:,n);          %select Zval within the analyzed box   
                                    reflpossmall=zeros(LEVELS,RAYS,N);  %create a matrix to find where Zval data exist          
                                    reflpossmall(reflsmall>0)=1;
                                    reflpos2Dsmall=sum(reflpossmall);   %where in 2Dvertical dimension data are     

                                    clear reflsmall reflpossmall

                                    minlatND = min(min(lat(reflpos2Dsmall>0)));  %to find the latitude limits    
                                    maxlatND = max(max(lat(reflpos2Dsmall>0)));         

                                    modMinLat = mod(minlatND, 0.05);         
                                    modMaxLat = mod(maxlatND, 0.05);         

                                    minlat = minlatND - modMinLat;         
                                    maxlat = maxlatND + (0.05 - modMaxLat);         

                                    clear minlatND maxlatND modMinLat modMaxLat

                                    if minlat < MinLatNorth
                                        minlat = MinLatNorth;
                                    end

                                    if maxlat > MaxLatNorth
                                        maxlat = MaxLatNorth;
                                    end

                                    minlongND = min(min(long(reflpos2Dsmall>0)));  %to find the longitude limits
                                    maxlongND = max(max(long(reflpos2Dsmall>0)));

                                    clear reflpos2Dsmall

                                    modMinLong = mod(minlongND, 0.05);
                                    modMaxLong = mod(maxlongND, 0.05);

                                    modMINLONGEAST = mod(MINLONGEAST, 0.05);
                                    modMAXLONGEAST = mod(MAXLONGEAST, 0.05);

                                    minlong = minlongND - modMinLong;
                                    maxlong = maxlongND + (0.05 - modMaxLong);

                                    clear minlongND maxlongND modMinLong modMaxLong

                                    if minlong<MINLONGEAST-modMINLONGEAST
                                        minlong=MINLONGEAST-modMINLONGEAST;
                                    end

                                    if maxlong>MAXLONGEAST+(0.05 - modMAXLONGEAST)
                                        maxlong=MAXLONGEAST+(0.05 - modMAXLONGEAST);
                                    end

                                    if minlong < MinLongEast
                                        minlong = MinLongEast;
                                    end

                                    if maxlong > MaxLongEast
                                        maxlong = MaxLongEast;
                                    end

                                    clear modMINLONGEAST modMAXLONGEAST MINLONGEAST MAXLONGEAST

                                    boundaries=[minlong maxlong minlat maxlat];
                                    boxsize1=(maxlong-minlong)*(maxlat-minlat);
                                    disp(['Box: ',num2str(boundaries)])
                                    clear boundaries

                                    [longitudeGridEast,latitudeGridEast] = ...
                                        meshgrid(minlong:GRIDBOXSTEP:maxlong,minlat:GRIDBOXSTEP:maxlat);

                                    %*****************************************************
                                    % Create the files with the information of pixels with
                                    % Z>0 systems found in the boxes previously analyzed      
                                    %*****************************************************
                                    if minlong<maxlong && minlat<maxlat && boxsize1>boxsize %check for boxsize at beginning
                                        clear boxsize1
                                        sizeGrid = size(longitudeGridEast);
                                        gridRows = sizeGrid(1,1);
                                        gridColumns = sizeGrid(1,2);
                                        numSwaths = N;
                                        size2D = numSwaths * RAYS;
                                        size3D = numSwaths * RAYS * LEVELS;

                                        %**************************************************************
                                        % Create a grid with the size of the new grid with 0 inside
                                        % the swath and 1 outside the swath. This is usefull because we     
                                        % dont have to interpolate outside the swath.
                                        %**************************************************************
                                        empty=zeros(size(longitudeGridEast));

                                        for i=1:sizeGrid(1,1)
                                            for j=1:sizeGrid(1,2)
                                                temp=find(lat>latitudeGridEast(i,j)-0.1 & lat<latitudeGridEast(i,j)+0.1...
                                                    & long>longitudeGridEast(i,j)-0.1 & long<longitudeGridEast(i,j)+0.1);
                                                if isempty(temp)
                                                    empty(i,j)=1;
                                                end
                                                temp=[];
                                            end
                                        end

                                        clear lat long sizeGrid temp  %delete the non corrected lat and long

                                        %******************************************
                                        % Create column versions of the regularly 
                                        % gridded latitude & longitude coordinates.
                                        %******************************************
                                        gridLongitude = reshape(longitudeGridEast, numel(longitudeGridEast), 1);
                                        gridLatitude = reshape(latitudeGridEast, numel(latitudeGridEast), 1);
                                        clear longitudeGridEast latitudeGridEast

                                        gridNoSwath = reshape(empty, numel(empty), 1);
                                        clear empty

                                        %***********************************
                                        % Create column versions of the data
                                        %***********************************
                                        sRLEVELSc = reshape(sR1, numel(sR1), 1);   %near surface rain for all scans
                                        sPLEVELSc = reshape(sP1, numel(sP1), 1);   %near surface phase for all scans
                                        rtLEVELSc = reshape(rt1, numel(rt1), 1);   %rain type for all scans         
                                        srtLEVELSc = reshape(srt1, numel(srt1), 1);%shallow rain type for all scans         
                                        wBBLEVELSc = reshape(wBB1, numel(wBB1), 1);%width bright band for all scans         
                                        hBBLEVELSc = reshape(hBB1, numel(hBB1), 1);%height bright band for all scans         
                                        clear sR1 sP1 rt1 srt1 wBB1 hBB1

                                        for level = 1:LEVELS
                                            eval(['Z',num2str(level),'c=reshape(Z',num2str(level),...
                                                ',prod(size(Z',num2str(level),')), 1);']);
                                            eval(['FE',num2str(level),'c=reshape(FE',num2str(level),...
                                                ',prod(size(FE',num2str(level),')), 1);']);
                                        end

                                        %***********************************************
                                        % Make reflectivity data ready for interpolation
                                        %***********************************************
                                        values = cat(1, Z176c, Z175c, Z174c, Z173c, Z172c, Z171c,...
                                            Z170c, Z169c, Z168c, Z167c, Z166c, Z165c, Z164c, Z163c, Z162c, Z161c,...
                                            Z160c, Z159c, Z158c, Z157c, Z156c, Z155c, Z154c, Z153c, Z152c, Z151c,...
                                            Z150c, Z149c, Z148c, Z147c, Z146c, Z145c, Z144c, Z143c, Z142c, Z141c,...
                                            Z140c, Z139c, Z138c, Z137c, Z136c, Z135c, Z134c, Z133c, Z132c, Z131c,...
                                            Z130c, Z129c, Z128c, Z127c, Z126c, Z125c, Z124c, Z123c, Z122c, Z121c,...
                                            Z120c, Z119c, Z118c, Z117c, Z116c, Z115c, Z114c, Z113c, Z112c, Z111c,...
                                            Z110c, Z109c, Z108c, Z107c, Z106c, Z105c, Z104c, Z103c, Z102c, Z101c,...
                                            Z100c, Z99c, Z98c, Z97c, Z96c, Z95c, Z94c, Z93c, Z92c, Z91c,...
                                            Z90c, Z89c, Z88c, Z87c, Z86c, Z85c, Z84c, Z83c, Z82c, Z81c,...
                                            Z80c, Z79c, Z78c, Z77c, Z76c, Z75c, Z74c, Z73c, Z72c, Z71c,...
                                            Z70c, Z69c, Z68c, Z67c, Z66c, Z65c, Z64c, Z63c, Z62c, Z61c,...
                                            Z60c, Z59c, Z58c, Z57c, Z56c, Z55c, Z54c, Z53c, Z52c, Z51c,...
                                            Z50c, Z49c, Z48c, Z47c, Z46c, Z45c, Z44c, Z43c, Z42c, Z41c,...
                                            Z40c, Z39c, Z38c, Z37c, Z36c, Z35c, Z34c, Z33c, Z32c, Z31c,...
                                            Z30c, Z29c, Z28c, Z27c, Z26c, Z25c, Z24c, Z23c, Z22c, Z21c,...
                                            Z20c, Z19c, Z18c, Z17c, Z16c, Z15c, Z14c, Z13c, Z12c, Z11c,...
                                            Z10c, Z9c,  Z8c,  Z7c,  Z6c,  Z5c,  Z4c,  Z3c,  Z2c,  Z1c);

                                        %******************************************************************
                                        % Convert Z factor to reflectivity [dBZ -> mm6/m3] 
                                        % (ScaleFactor is a remnant from v0 when using extracted hdf files)
                                        %******************************************************************
                                        reflectivity = 10.^(values./double(10*correctZ_ScaleFactor));
                                        reflectivity(values < 0) = NODATA;

                                        %********************************************
                                        % Make flag echo data ready for interpolation
                                        %********************************************
                                        valuesFE = cat(1, FE176c, FE175c, FE174c, FE173c, FE172c, FE171c,...
                                            FE170c, FE169c, FE168c, FE167c, FE166c, FE165c, FE164c, FE163c, FE162c, FE161c,...
                                            FE160c, FE159c, FE158c, FE157c, FE156c, FE155c, FE154c, FE153c, FE152c, FE151c,...
                                            FE150c, FE149c, FE148c, FE147c, FE146c, FE145c, FE144c, FE143c, FE142c, FE141c,...
                                            FE140c, FE139c, FE138c, FE137c, FE136c, FE135c, FE134c, FE133c, FE132c, FE131c,...
                                            FE130c, FE129c, FE128c, FE127c, FE126c, FE125c, FE124c, FE123c, FE122c, FE121c,...
                                            FE120c, FE119c, FE118c, FE117c, FE116c, FE115c, FE114c, FE113c, FE112c, FE111c,...
                                            FE110c, FE109c, FE108c, FE107c, FE106c, FE105c, FE104c, FE103c, FE102c, FE101c,...
                                            FE100c, FE99c, FE98c, FE97c, FE96c, FE95c, FE94c, FE93c, FE92c, FE91c,...
                                            FE90c, FE89c, FE88c, FE87c, FE86c, FE85c, FE84c, FE83c, FE82c, FE81c,...
                                            FE80c, FE79c, FE78c, FE77c, FE76c, FE75c, FE74c, FE73c, FE72c, FE71c,...
                                            FE70c, FE69c, FE68c, FE67c, FE66c, FE65c, FE64c, FE63c, FE62c, FE61c,...
                                            FE60c, FE59c, FE58c, FE57c, FE56c, FE55c, FE54c, FE53c, FE52c, FE51c,...
                                            FE50c, FE49c, FE48c, FE47c, FE46c, FE45c, FE44c, FE43c, FE42c, FE41c,...
                                            FE40c, FE39c, FE38c, FE37c, FE36c, FE35c, FE34c, FE33c, FE32c, FE31c,...
                                            FE30c, FE29c, FE28c, FE27c, FE26c, FE25c, FE24c, FE23c, FE22c, FE21c,...
                                            FE20c, FE19c, FE18c, FE17c, FE16c, FE15c, FE14c, FE13c, FE12c, FE11c,...
                                            FE10c, FE9c,  FE8c,  FE7c,  FE6c,  FE5c,  FE4c,  FE3c,  FE2c,  FE1c);
                                        flagEchoSub = valuesFE; 
                                        flagEchoSub(valuesFE < 0) = NODATA;

                                        %*******************************************************
                                        % Make latitude & longitude data ready for interpolation
                                        %*******************************************************
                                        fliplat=flip(lats,3);
                                        latitude=reshape(fliplat(:,:,1:LEVELS),[],1);

                                        fliplong=flip(longs,3);
                                        longitude=reshape(fliplong(:,:,1:LEVELS),[],1);

                                        clear values valuesFE fliplat fliplong

                                        %*********************
                                        % Call interp function
                                        %*********************
                                        disp('Starting refl interpolation');
                                        interp_type = 'weighted';
                                        gridRefl = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            LEVELS,RADIUS,interp_type,latitude',longitude',        ...
                                            reflectivity',gridLatitude',gridLongitude',gridNoSwath');
                                        gridRefl(gridRefl == NODATA) = NaN;
                                        gridRefl = 10 * log10(gridRefl);
                                        gridRefl(gridRefl < MINDBZ) = NODATA;
                                        gridRefl(isnan(gridRefl)) = NODATA;
                                        gridRefl = gridRefl';

                                        %**********CHECK BETTER WAY TO DO THIS**********
                                        % can I save nn index for each grid point to make this faster?
                                        %**********CHECK BETTER WAY TO DO THIS**********
                                        disp('Starting flagEcho interpolation');
                                        interp_type = 'nn';
                                        gridFlagEcho = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            LEVELS,RADIUS,interp_type,latitude',longitude', ...
                                            flagEchoSub',gridLatitude',gridLongitude',gridNoSwath');
                                        gridFlagEcho = gridFlagEcho';

                                        disp('Starting 2D nearest neighbor interpolations');
                                        interp_type = 'nn';
                                        gridSurfRain = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            1,RADIUS,interp_type,latitude',longitude',        ...
                                            sRLEVELSc',gridLatitude',gridLongitude',gridNoSwath');
                                        gridSurfRain = gridSurfRain';
                                        gridRainType = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            1,RADIUS,interp_type,latitude',longitude',        ...
                                            rtLEVELSc',gridLatitude',gridLongitude',gridNoSwath');
                                        gridRainType = gridRainType';
                                        gridSurfPhase = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            1,RADIUS,interp_type,latitude',longitude',        ...
                                            sPLEVELSc',gridLatitude',gridLongitude',gridNoSwath');
                                        gridSurfPhase = gridSurfPhase';
                                        gridShallowRainType = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            1,RADIUS,interp_type,latitude',longitude',        ...
                                            srtLEVELSc',gridLatitude',gridLongitude',gridNoSwath');
                                        gridShallowRainType = gridShallowRainType';
                                        gridWidthBB = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            1,RADIUS,interp_type,latitude',longitude',        ...
                                            wBBLEVELSc',gridLatitude',gridLongitude',gridNoSwath');
                                        gridWidthBB = gridWidthBB';
                                        gridHeightBB = interp_for_matlab(size2D,gridColumns,gridRows,  ...
                                            1,RADIUS,interp_type,latitude',longitude',        ...
                                            hBBLEVELSc',gridLatitude',gridLongitude',gridNoSwath');
                                        gridHeightBB = gridHeightBB';
                                        clear sRLEVELSc rtLEVELSc sPLEVELSc srtLEVELSc wBBLEVELSc hBBLEVELc;

                                        %****************************************************************************************
                                        % 1. Add a new variable called gridRainTypeRaw which retains the orig gridRainType values
                                        % 2. Assign gridRainType values of 1 (stra),2 (conv), 3 (other)
                                        %****************************************************************************************
                                        gridRainTypeRaw=gridRainType;

                                        gridRainType(gridRainType>0 & fix(gridRainType/10000000.)==1) = 1;
                                        gridRainType(gridRainType>0 & fix(gridRainType/10000000.)==2) = 2;
                                        gridRainType(gridRainType>0 & fix(gridRainType/10000000.)==3) = 3;

                                        %****************************************************************************************
                                        % 1. Add a new variable called gridSurfPhaseRaw which retains the orig gridSurfPhase vals
                                        % 2. Assign gridSurfPhase values of 0 (solid), 1 (mixed phase), 2 (liquid) to all values
                                        %    that are not negative or missingvals (255)
                                        %****************************************************************************************
                                        gridSurfPhaseRaw=gridSurfPhase;

                                        gridSurfPhase(gridSurfPhase>0 & gridSurfPhase~=255 & fix(gridSurfPhase/100)==0) = 0;
                                        gridSurfPhase(gridSurfPhase>0 & gridSurfPhase~=255 & fix(gridSurfPhase/100)==1) = 1;
                                        gridSurfPhase(gridSurfPhase>0 & gridSurfPhase~=255 & fix(gridSurfPhase/100)==2) = 2;

                                        %******************************************************
                                        % call functions to prepare data and output netcdf file
                                        %******************************************************
                                        prepare_data

                                        dir_out=[out_dir,date_box_start(1:4),'/',date_box_start(5:6),'/'];
                                        if ~exist(dir_out,'dir')
                                           mkdir(dir_out)
                                        end

                                        make_gpm_cdf_file_v1a

                                        clear groundTrackDate* date_box* dir_out

                                    else
                                        fprintf('\n')
                                        disp('Empty box or box too small.')
                                        clear boxsize1
                                    end % end of boxsize restriction

                                    clear minlat minlong
                                    clear maxlat maxlong
                                    clear fid* fileInfo
                                    clear gridRows gridColumns
                                    clear gridLongitude gridLatitude
                                    clear gridRainType gridRainTypeRaw
                                    clear gridShallowRainType gridSurfRain 
                                    clear gridSurfPhase gridSurfPhaseRaw
                                    clear gridWidthBB gridHeightBB
                                    clear flagEchoSub gridFlagEcho
                                    clear reflectivity gridRefl
                                    clear Z* FE*
                                    clear i j k
                                    clear lats latitude
                                    clear longs longitude
                                    clear n N
                                    clear noSwath
                                    clear numSwaths
                                    clear fx1 fx2
                                end % endif of indices check
                            end % endif of restriction of limits being to close to the boundaries    
                        end  %endfor of analyzed regions
                    end   %endif of regions found having contigous data restriction for more than 1 region!!!
                end  %endif for orbits without correctZfact data in area of interest
                clear correctZ rainTypes shallowRainTypes
                clear nearSurfRain nearSurfPhase
                clear widthBB heightBB flagEcho
                clear hour second minute scanTimeAll 
                clear latitude* longitude*
                clear namein2Ku*
            else
                fprintf('\n')
                disp(['The file ',namein2Ku,' is not included in the limits']);
	
                clear filename namein2Ku* latitude longitude
                clear geolocation0 longvector latvector wholeLongs0 wholeLats0
                clear x1 x2 x3 em
            end    %endif for em=1 so no orbits within the area
        end  %end of for for different orbit files
        
        % close log file for year and month
        fclose(logID);
        
    end  %loop through the different months
    
end  %loop through the different years

t=clock;
datestamper(t);
