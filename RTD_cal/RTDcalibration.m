function coeffs= RTDcalibration(files, numRTDs, plotsYorN)
% Determines a linear least-squares fit for RTD temperature data in terms
% of steady-state calibration temperatures.
%
% Takes inputs 'files', which should be a string of CSV file names that
% contain all necssary data (file names INCLUDE '.csv' and names are 
% separated by spaces), 'numRTDs', which is the number of RTDs that
% are to be calibrated, and 'plotsYorN', which generates plots of T_RTD
% vs T_cal if equal to 'Y' and does NOT generate plots if equal to 'N'.
%
% Outputs 'coeffs', which is a cell array of length numRTDs, where each
% entry is a vector of least squares coefficients [a,b] that correspond to
% the linear fit T_RTD = a * T_cal + b.
%
% IMPORTANT: guidance for setup of CSV data files!
% 1) Put column names in the first row. 
% 2) Column 3 should be labeled 't (msec)' and contains the total number of
% milliseconds passed since the beginning of data-taking for a given
% measurement.
% 3) Column 4 should be labeled 'T_cal' and contains the calibrated
% temperature reading at the time of the given measurement.
% 4) Columns 8 through numRTDs+7 should be labeled 'RTD_[#]' and contain
% the temperature reading of each RTD at the time of the given measurement.
% NOTE: plot titles will not take into account whether numbers have been
% skipped in labeling RTDs.

dTemp= 0.15; % allowed temperature variation for steady-state
minTimeForSS= 10*60*1000; % minimum length (in seconds) of a steady-state 
                            % temperature interval

%files= input('Enter names of CSV files separated by spaces, INCLUDING ''.csv'': ', 's');
%numRTDs= input('Enter the number of RTDs to be calibrated: ', 's');
%numRTDs= str2double(numRTDs);

idxOfSpaces= [];
for i= 1:length(files)
    if files(i) == char(32)
        idxOfSpaces(end+1)= i;
    end    
end    

numFiles= length(idxOfSpaces) + 1; % number of CSV data files supplied by user

allT_cal= []; % steady-state calibrated temperature data from ALL files
allRTDs= []; % RTD temperature data recorded during steady-state from ALL files
for i= 1:numFiles
    if i==1
        fileNames= files(1:idxOfSpaces(1)-1);
    elseif i==numFiles
        fileNames= files(idxOfSpaces(end)+1:end);
    else
        fileNames= files(idxOfSpaces(i-1)+1:idxOfSpaces(i)-1);
    end    
    
    T= readtable(fileNames); % read data from i'th file into table
    
    % obtain steady-state intervals and row indices with anomalous data
    [intervals, anoms]= desiredIntervals(T, dTemp, minTimeForSS);
    
    % convert from table to cell array, remove anomalous data from array,
    % and then remove non-steady-state data from array
    C= removeAnomaliesFromArray(table2cell(T), anoms);
    steadyC= joinSSIntervals(C, intervals);
    
    T_cal= concatT_cal(steadyC, 4);
    RTDs= concatRTDs(steadyC, numRTDs);
    
    if ~isempty(T_cal)
        allT_cal= [allT_cal; T_cal];
        allRTDs= [allRTDs; RTDs];
    end
end   

coeffs= cell(numRTDs, 1);
for i= 1:numRTDs
    % determine least squares fit and plot
    [a,b]= leastSquaresApprox(allT_cal, allRTDs(:,i));
    coeffs{i}= [a b];
    if strcmp(plotsYorN, 'Y')
        figure;
        plot(allT_cal, allRTDs(:,i), '.b');
        title(['T_{RTD} vs T_{cal}: RTD ' num2str(i) '; a = ' num2str(a) ...
            ', b = ' num2str(b)])
        xlabel('T_{cal}'); ylabel('T_{RTD}');
    end    
end 

end

function [a,b]= leastSquaresApprox(T_cal, RTDs)
% Determines a least-squares fit, T_RTD = a * T_cal + b.

    X= [T_cal ones(length(T_cal),1)];
    solns= (X'*X) \ (X'*RTDs);
    a= solns(1); b= solns(2);
end

function RTDs= concatRTDs(steadyC, numRTDs)
% Concatenates RTD temperature data corresponding to steady-state
% calibrated temperature data into one matrix, RTDs.

    [len,~]= size(steadyC);
    RTDs= zeros(len, numRTDs);
    for i= 1:len
        for j= 8:numRTDs+7
            RTDs(i,j-7)= steadyC{i,j}; 
        end
    end
end

function T_cal= concatT_cal(steadyC, T_calColNum)
% Concatenates steady-state calibrated temperature data into one vector,
% T_cal.

    [len,~]= size(steadyC);
    if len==0
        T_cal= [];
        return
    end
    T_cal= zeros(len,1);    
    for i= 1:len
        T_cal(i)= steadyC{i,T_calColNum}; 
    end 
end

function steadyC= joinSSIntervals(C, intervals)
% Joins all steady-state intervals together into one data array steadyC,
% removing rows that do not correspond to steady-state.

    [~,w]= size(C);

    steadyC= {};
    idx= 1;
    for i= 1:1:length(intervals)
        for j= intervals{i}(1):intervals{i}(2)
            for k= 1:w
                steadyC{idx,k}= C{j,k};
            end
            idx= idx+1;
        end    
    end    
end

function C= removeAnomaliesFromArray(tempC,anoms)
% Removes rows from data matrix tempC if rows have anomalous calibration
% temperatures.

    [len,w]= size(tempC);
    C= cell(len-length(anoms),w);

    idx= 1;
    for i= 1:len
        if ~ismember(i,anoms)
            for j= 1:w
                C{idx,j}= tempC{i,j};
            end    
            idx= idx+1;
        end    
    end
end

function [intervals, anoms]= desiredIntervals(C, dT, minTimeForSS)
% Returns intervals (i.e. ranges of row indices) corresponding to steady-state
% calibration temperatures, as well as the indices of rows whose
% calibration temperature is anomalous.

    intervals= {};
    
    tempTimes= C{:,'t_msec_'};
    tempT_cal= C{:,'T_cal'};
    
    [times, T_cal,anoms]= removeAnomaliesFromCols(tempTimes,tempT_cal);
    
    i= 1;
    beginInterval= 1;
    while i <= length(T_cal) && times(end)-times(i) >= minTimeForSS
        for j= i:length(T_cal)
            if beginInterval
                beginInterval= 0;
                intervalStartTime= times(i);
                intervalStartTemp= T_cal(i);
                intervalStartIdx= i;
                continue
            end    
            
            if ~beginInterval && abs(T_cal(j)-intervalStartTemp) <= dT
                if j < length(T_cal)
                    continue
                else
                    intervals{end+1}= [intervalStartIdx j];
                    i= j+1;
                end    
            elseif ~beginInterval && abs(T_cal(j)-intervalStartTemp) > dT
                if j < length(T_cal) && abs(T_cal(j+1)-intervalStartTemp) <= dT
                    continue
                end
                
                if times(j)-intervalStartTime < minTimeForSS
                    i= i+1;
                else
                    intervals{end+1}= [intervalStartIdx j-1];
                    i= j+1;
                end 
                beginInterval= 1;
                break
            end    
        end
    end    
end

function [times,T_cal,anoms]= removeAnomaliesFromCols(tempTimes,tempT_cal)
% Removes rows from both T_cal and times, based on the determination of
% whether a given calibration temperature is anomalous.

    anoms= [];
    for i= 2:length(tempT_cal)
        if abs(tempT_cal(i)-tempT_cal(i-1)) > 0.25 && ~ismember(i-1,anoms)
            anoms(end+1)= i;
        elseif ismember(i-1,anoms) && tempT_cal(i)==tempT_cal(i-1)
            anoms(end+1)= i;
        end    
    end    
    
    times= zeros(length(tempTimes)-length(anoms),1);
    T_cal= zeros(length(tempT_cal)-length(anoms),1);
    c= 1;
    for i= 1:length(tempT_cal)
        if ~ismember(i,anoms)
            times(c)= tempTimes(i);
            T_cal(c)= tempT_cal(i);
            c= c+1;
        end    
    end    
end
