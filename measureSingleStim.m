function singleStim = measureSingleStim(files,group,plotFlag)
%% Measure latencies for single-pulse tactile and electrical stimulation
% This function obtains the first-spike-latencies in response to the single
% electrical and tactile pulses (not the stimulus pairs) in the collide
% stimulation protocoll (Steffen, v4). Latency is defined as the time
% difference between the stimulus onset and the peak of the spike. For the
% cold experiments the returned latencies are sorted according to the 
% temperature (low to high). All latencies are returned next to the mean 
% latency per trial (calculated based on the 2nd, 3rd, and 4th pulse).
% Additionally, 10 °C-change-rates are calculated (similar to Q10) based on
% the mean latencies (only for cold experiments). 
% If requested the latencies and the temperature information is summarized 
% visually for each file. 
% -----
% Input
% -----
%   files       1*n vector, the number of the files to be analyzed
%   group       string, stating which experimental group should be analyzed, 
%               possible options: "Cold", "RoomTemp"
%   plotFlag    1*1 logical, indicate whether to plot the results in one
%               figure per file, 1 - visual output, 0 - no visual output
% ------
% Output
% ------
%   singleStim  n*10 cell array, containing all analysis results
%                   1. column - 4*m matrix with the tactile latencies for
%                               4 stimuli and m trials
%                   2. column - 4*m matrix with the electrical latencies 
%                               for 4 stimuli and m trials
%                   3. column - 1*m vector containing the median tactile
%                               latencies for m trials
%                   4. column - 1*m vector containing the median electrical
%                               latencies for m trials 
%                   5. column - 4*m matrix with the number of spikes 
%                               for each tactile pulse
%                   6. column - 4*m matrix with the number of spikes
%                               for each electrical pulse
%                   7. column - 1*m vector with the sorted temperature
%                               values
%                   8. column - 1*1 double, 10 °C change rate for the
%                               tactile spike latency
%                   9. column - 1*1 double, 10 °C change rate for the 
%                               electrical spike latency 
%                   10. column - 1*m vector, sort index (according to
%                               temperature)
% 
%   Please note for the cold experiments: 
%                The trials are sorted by temperature (low to high)!
%                If it is necessary to restore the original recording
%                order, the sort index can be used to do so (sort the index
%                values, obtain the new sort index and apply to the data).
% ----
% Helper function required: 'singlePulseLatency'
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 12.02.2023
% -------------------------------------------------------------------------

%% Preparations
% Store the master-path where the data folders are stored
masterPath = ['C:\Users\schul\Nextcloud\Studium\2. Master-Semester\' ...
    'ResearchModule\'];

% Paste search path depending on the dataset
if isequal(group,"Cold")
    dataPath = masterPath+"Cooled\";
elseif isequal(group,"RoomTemp")
    dataPath = masterPath+"RoomTemp\";
% Throw error if non-existing group name was given
else
    % Error with appropriate message
    error("Experiment group must be one of the following: 'Cold', " + ...
        "'RoomTemp'")
end   

% List all files in the folder
fileList = dir(dataPath+"tacData*.mat");
% Load a file with the trials to select
load(dataPath+'tempSelect.mat','tempSelect')

% Initialize a cell array to collect all results
singleStim = cell(length(files),10);

% Iterate through all files in the file list
for file = files
    % Load the next file and display name in command window
    load(dataPath+fileList(file).name,'tacData','metadata')
    fprintf("Currently working on file: %s!\n",fileList(file).name)

    % Reduce the recording data according to the 'tempSelect' variable
    tacData = tacData(tempSelect(file,1):tempSelect(file,2),:);

    % extract the stimulation protocol for tactile and electrical 
    % stimulation
    elecStim = tacData{1,1}(:,1);
    tacStim = tacData{1,1}(:,2);
    
    % Find the onsets of the stimulus pulses (tactile and electrical)
    [~,elecOnsets] = findpeaks(elecStim,metadata.timeVector);
    [~,tacOnsets] = findpeaks(tacStim,metadata.timeVector);

    % Obtain the onsets where the cell was either stimulated only 
    % electrically or tactilely (found by the spacing of the stimuli 
    % 0.1 s / 100 ms cutoff distance)
    onlyElec = elecOnsets(elecOnsets-tacOnsets < -0.1);
    onlyTac = tacOnsets(tacOnsets-elecOnsets > 0.1);

    % Calculate latencies for tactile and electrical stimuli in all trials
    % (values are retured in ms)
    [elecLatencies,tacLatencies,nElec,nTac] = cellfun( ...
        @(x)singlePulseLatency(x,onlyElec,onlyTac,0.1, ...
        metadata.timeVector),tacData(:,2),'UniformOutput',false);
    
    % Convert cell arrays to matrices
    elecLatencies = cell2mat(elecLatencies');
    tacLatencies = cell2mat(tacLatencies');
    nElec = cell2mat(nElec');
    nTac = cell2mat(nTac');
    
    % Extract the temperature data
    temperature = cell2mat(tacData(:,3));

    if isequal(group,'Cold')
        % Sort the temperature, reorder the latencies accordingly
        [temperatureSorted, sortIndex] = sort(temperature);
        tacLatenciesSorted = tacLatencies(:,sortIndex);
        elecLatenciesSorted = elecLatencies(:,sortIndex);
        % Also sort the spike counts accordingly
        nTacSorted = nTac(:,sortIndex);
        nElecSorted = nElec(:,sortIndex);
        
        % Calculate mean latencies as well as minimum and maximum latencies
        meanTacLat = median(tacLatenciesSorted,1,'omitnan');
        minTacLat = min(tacLatenciesSorted);
        maxTacLat = max(tacLatenciesSorted);
        % ... as well for the latencies for the electrical stimulaiton
        meanElecLat = median(elecLatenciesSorted,1,'omitnan');
        minElecLat = min(elecLatenciesSorted);
        maxElecLat = max(elecLatenciesSorted);
   
        % Calculate the Q10 values based on the latencies
        Q10valTac = (tacLatenciesSorted(2,end)/ ...
            tacLatenciesSorted(2,1))^(10/ ...
            (temperatureSorted(end)-temperatureSorted(1)));
        Q10valElec = (elecLatenciesSorted(2,end)/ ...
            elecLatenciesSorted(2,1))^(10/ ...
            (temperatureSorted(end)-temperatureSorted(1)));

        % Collect all results in cell array
        singleStim{file,1} = tacLatenciesSorted;
        singleStim{file,2} = elecLatenciesSorted;
        singleStim{file,3} = meanTacLat;
        singleStim{file,4} = meanElecLat;
        singleStim{file,5} = nTacSorted;
        singleStim{file,6} = nElecSorted;
        singleStim{file,7} = temperatureSorted;
        singleStim{file,8} = Q10valTac;
        singleStim{file,9} = Q10valElec;
        singleStim{file,10} = sortIndex;
    else
        % Calculate mean latencies as well as minimum and maximum latencies
        meanTacLat = median(tacLatencies,1,'omitnan');
        minTacLat = min(tacLatencies);
        maxTacLat = max(tacLatencies);
        % ... as well for the latencies for the electrical stimulaiton
        meanElecLat = median(elecLatencies,1,'omitnan');
        minElecLat = min(elecLatencies);
        maxElecLat = max(elecLatencies);

        % Collect all results in cell array
        singleStim{file,1} = tacLatencies;
        singleStim{file,2} = elecLatencies;
        singleStim{file,3} = meanTacLat;
        singleStim{file,4} = meanElecLat;
        singleStim{file,5} = nTac;
        singleStim{file,6} = nElec;
        singleStim{file,7} = temperature;
        singleStim{file,8} = nan;
        singleStim{file,9} = nan;
        singleStim{file,10} = (1:length(meanElecLat))';
    end

    %% Graphical presentation (file-wise)
    if plotFlag
        % Create a new graphic window
        figure('Name',metadata.date,'WindowState','maximized');
        % Set up a tiled layout
        tiledlayout(2,3,'TileSpacing','compact')
        
        nexttile    % 1st sub-plot
        % Plot the tactile spike latencies over the sorted temperature
        plot(temperatureSorted,tacLatenciesSorted,'o')
        % Label sub-plot
        title("Tactile stimulation")
        ylabel("Latency [ms]")
        xlabel("Temperature [°C]")
        
        nexttile(4)  % 4th sub-plot
        % Fill the area between minimal and maximal values of each trial
        fill([temperatureSorted;temperatureSorted(end:-1:1)], ...
            [minTacLat,maxTacLat(end:-1:1)],[0,0.6,1], ...
            'FaceAlpha',0.2,'LineStyle','none')
        % Add the mean tactile spike latencies as a line
        hold on; plot(temperatureSorted,meanTacLat,'Color', [0,0.6,1], ...
            'LineWidth',1); hold off

        % Label sub-plot
        title("Tactile stimulation")
        legend("min-max","mean")
        xlabel("Temperature [°C]")
        ylabel("Latency [ms]")
        
        nexttile(2)     % 2nd sub-plot
        % Plot the electrical spike latencies over the sorted temperature
        plot(temperatureSorted,elecLatenciesSorted,'o')
        % Label the sub-plot
        title("Electrical stimulation")
        ylabel("Latency [ms]")
        xlabel("Temperature [°C]")
        
        nexttile(5)     % 5th sub-plot
        % Fill the area between the minimal and maximal latnecies of each 
        % trial
        fill([temperatureSorted;temperatureSorted(end:-1:1)], ...
            [minElecLat,maxElecLat(end:-1:1)],[1,0.6,0], ...
            'FaceAlpha',0.2,'LineStyle','none')
        % Add the mean electrial spike latencies as a line
        hold on; plot(temperatureSorted,meanElecLat,'Color', [1,0.6,0], ...
            'LineWidth',1); hold off

        % Label sub-plot
        title("Electrical stimulation")
        legend("min-max","mean")
        xlabel("Temperature [°C]")
        ylabel("Latency [ms]")
        
        nexttile(3,[2,1])   % 3rd sub-plot (spanning 2 tiles on the right)
        % Plot the temperature over the trial index
        plot(temperature)
        % Label the sub-plot
        title("Temperature course")
        xlabel("Trial index")
        ylabel("Temperature [°C]")

    % End of optional plot environment
    end
% End of file iteration
end

% End of function definition
end
