%Name: dateSorterGPM.m
%Purpose: gives the correct date and time (UTC) for a pice of GPM orbit
%Latest update: 2007/11/29
%By: Darren Wilton, (corrected by Ulrike Romatschke
%correction Manuel Zuluaga (09/21/2011))
%correction Stacy Brodzik (03/09/2017))
%Used in: getGPM.m

% there are potential date errors with this file ... --> shouldn't be anymore ;-) (Ulli)

function [groundTrackDate] = dateSorterGPM(centreScan, dayFile, month, year, firstTime, lastTime, scanTime)

timeIndex = centreScan;

    if firstTime < lastTime
        day2 = dayFile;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month, day2, gtHour, gtMinute, gtSecond);      
            
    elseif firstTime > lastTime && scanTime(timeIndex) > lastTime
        day2 = dayFile;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==1 && dayFile==31
        month2 = 2;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==2 && dayFile==28 && mod(year,4)==0
        month2 = 2;
        day2 = 29;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==2 && dayFile==28 && mod(year,4)~=0
        month2 = 3;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==2 && dayFile==29
        month2 = 3;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==3 && dayFile==31
        month2 = 4;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==4 && dayFile==30
        month2 = 5;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==5 && dayFile==31
        month2 = 6;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==6 && dayFile==30
        month2 = 7;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==7 && dayFile==31
        month2 = 8;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==8 && dayFile==31
        month2 = 9;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==9 && dayFile==30
        month2 = 10;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==10 && dayFile==31
        month2 = 11;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==11 && dayFile==30
        month2 = 12;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) < lastTime && month==12 && dayFile==31
        year2 = year+1;
        month2 = 1;
        day2 = 1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year2, month2, day2, gtHour, gtMinute, gtSecond);
        
    elseif firstTime > lastTime && scanTime(timeIndex) <= lastTime
        day2 = dayFile+1;
        gtHour = floor(scanTime(timeIndex)/3600);
        gtMinute = floor((mod(scanTime(timeIndex),3600))/60);
        gtSecond = scanTime(timeIndex)-(gtHour*3600)-(gtMinute*60);
        groundTrackDate = datenum(year, month, day2, gtHour, gtMinute, gtSecond);
                
    end
