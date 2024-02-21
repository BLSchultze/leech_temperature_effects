function [combiStimLat,combiLatSplit,stimTimeDiff] = ...
    measureCombiStim(files,group)
%% Measure spike latencies for the tac-elec collide stimulation protocol
% This function obtains the latencies of the spikes that were elicited by
% the 19 tactile and electrical stimulus pairs in the collide stimulation
% protocol (Steffen, v4). The spike latency is defined here as the time
% difference between the tactile stimulus onset and the peak value of the
% spike. The stimulus pairs are described by the relative timing of the
% electrical stimulus to the tactile stimulus (onset elec - onset tac). If
% the electrical stimulus came first, the time difference is negative.
% Accordingly, spikes that occurred before the tactile stimulus onset, have
% negative latencies. The same analysis can be performed for the cold and
% the room temperature datasets. The trials are re-arranged by temperature 
% for the cold experiments.
% The lateencies are returend in two forms: 1. latencies for all spikes
% occuring during one paired stimulation, 2. latencies for the 1st, 2nd and
% 3rd spike ocurring during one paired stimulation.
% -----
% Input
% -----
%   files   1*n vector, the number of the files to be analyzed
%   group   string, stating which experimental group should be analyzed, 
%           possible options: "Cold", "RoomTemp"
% ------
% Output
% ------
%   combiStimLat    n*19 cell array, contains all spike latencies 
%                   [ms] for the collide stimuli (tactile and electrical 
%                   stimuli), 1st dimension: file, 2nd dimension: 19 
%                   collide stimulus pairs. Spike latencies are stored in
%                   m*1 cell arrays, seperate for all m trials of a file.
%   combiLatSplit   n*5 cell array, contains the spike latencies [ms]
%                   splittted in 1st, 2nd and 3rd spike, 1st dimenstion:
%                   file, 2nd dimension: 
%                       1. column - 1st spike latency
%                       2. column - 2nd spike latency
%                       3. column - 3rd spike latency
%                       4. column - temperature per trial, sorted
%                       5. column - sort index (according to temperature)
%   stimTimeDiff    1*19 vector, contains the time difference between the
%                   19 tactile and electrical pulses (elec - tac) in ms
%
%    Please note for the cold experiments: 
%                 The trials are sorted by temperature (low to high)!
%                 If it is necessary to restore the original recording
%                 order, the sort index can be used to do so (sort the 
%                 index values, obtain the new sort index and apply to the 
%                 data).
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 18.03.2023
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

% Initialize a cell array to collect all latencies
combiStimLat = cell(length(files),19,1);
combiLatSplit = cell(length(files),3);

% Define a upper limit for the time window to search for spikes after a
% stimulation (added to the onset of the second stimulus) [s]
searchLen = 0.15;

%% Latency calculations and organizing

% Iterate through all data files
for file = files
    % Load the next file and display name in command window
    load(dataPath+fileList(file).name,'tacData','metadata')
    fprintf("Currently working on file: %s!\n",fileList(file).name)

    % Reduce the recording data according to the 'tempSelect' variable
    tacData = tacData(tempSelect(file,1):tempSelect(file,2),:);

    % Calculate the sampling rate
    samplingRate = metadata.samplingRate;

    % extract the stimulation protocol for tactile and electrical 
    % stimulation
    elecStim = tacData{1,1}(:,1);
    tacStim = tacData{1,1}(:,2);
    
    % Find the onsets of the stimulus pulses (for electrical and tactile
    % stimuli)
    [~,elecOnsets] = findpeaks(elecStim,metadata.timeVector);
    [~,tacOnsets] = findpeaks(tacStim,metadata.timeVector);
    % Obtain the onsets where the cell was either stimulated electrically
    % or tactilely by spacing of the stimuli (cutoff 0.1 s / 100 ms)
    onlyElec = elecOnsets(elecOnsets-tacOnsets < -0.1);
    onlyTac = tacOnsets(tacOnsets-elecOnsets > 0.1);

    % Select only the stimuli that occurred as pairs (1 tactile, 1
    % electrical)
    elec = elecOnsets(~ismember(elecOnsets,onlyElec));
    tac = tacOnsets(~ismember(tacOnsets,onlyTac));

    % Calculate the difference between the two stimuli of each stimulus
    % pair [ms]
    stimTimeDiff = (elec - tac)*1000;

    % Initialize/empty the variables to collect the spike latencies
    firstSpikeLat = zeros(length(tacData),length(stimTimeDiff));
    secondSpikeLat = firstSpikeLat;
    thirdSpikeLat = firstSpikeLat;

    % Sort according to temperature for the cold experiments
    if isequal(group,'Cold')
        % Sort the temperature values 
        [temperature_srt,sortIndex] = sort(cell2mat(tacData(:,3)));
    else
        % Do not sort for the room temperature experiments
        sortIndex = (1:length(tacData(:,3)))';
        temperature_srt = tacData(:,3);
    end

    % Set a trial counter variable to sort the latencies while saving
    newIndex = 0;

    % Iterate through all trials
    for trial = sortIndex'
        % Increase trial counter
        newIndex = newIndex + 1;

        % Find all spikes in the recording trace
        [peakValues,spikeIndices] = findpeaks(tacData{trial,2}, ...
            'MinPeakProminence',30,'MinPeakDistance',0.002*samplingRate);
        
        % Initialize a spike counter/index
        spikecount = 0;
        % Create a variable to collect the spike times
        spikeTimes = [];
        
        % Iterate through all spikes
        for spike = 1:length(spikeIndices)
            % Only execute if value 1.2 ms before spike is still in the 
            % recording
            if spikeIndices(spike) - 0.0012*samplingRate > 1
                % Get the membrane potential 1.2 ms before and after the 
                % peak
                l = tacData{trial,2}(spikeIndices(spike) - ...
                    0.0012*samplingRate);
                r = tacData{trial,2}(spikeIndices(spike) + ...
                    0.0012*samplingRate);
                % Calculate the difference to the peak value
                peakDiffl = peakValues(spike) - l;
                peakDiffr = peakValues(spike) - r;
                
                % Identify spikes by a membrane potential change of 15 mV 
                % on the left and 20 mV on the right side of the peak
                if peakDiffl > 4 && peakDiffr > 4
                    % Increase spike counter/index variable
                    spikecount = spikecount + 1;
                    % Store spike time
                    spikeTimes(spikecount) = ...
                        metadata.timeVector(spikeIndices(spike));
                end
            end
        end


        % Iterate through all stimulus pairs
        for stim = 1:length(stimTimeDiff) 
            % For all stimulations where the tactile stimulus came first or
            % at the same time as the electrical ...
            if tac(stim) <= elec(stim)
                % Select the spikes following the current stimulus pair
                currentSpikes = spikeTimes(spikeTimes > tac(stim) & ...
                    spikeTimes < tac(stim)+searchLen);
            else
                % Select the spikes following the current stimulus pair
                currentSpikes = spikeTimes(spikeTimes > elec(stim) & ...
                    spikeTimes < elec(stim)+searchLen);
            end

            % Calculate the spike latencies relative to the tactile 
            % stimulus
            latencies = currentSpikes - tac(stim);
            % Store all latencies in a cell array
            combiStimLat{file,stim}{newIndex,1} = latencies;

            % Split the latencies according to the spike order
            if isempty(latencies)
                % Assign NaN if no spikes were found
                firstSpikeLat(newIndex,stim) = nan;
                secondSpikeLat(newIndex,stim) = nan;
                thirdSpikeLat(newIndex,stim) = nan;
            elseif length(latencies) == 1
                % First spike latency, if one spike was found
                firstSpikeLat(newIndex,stim) = latencies(1);
                secondSpikeLat(newIndex,stim) = nan;
                thirdSpikeLat(newIndex,stim) = nan;
            elseif length(latencies) == 2
                % First and second spike latencies, if two spikes were
                % found
                firstSpikeLat(newIndex,stim) = latencies(1);
                secondSpikeLat(newIndex,stim) = latencies(2);
                thirdSpikeLat(newIndex,stim) = nan;
            elseif length(currentSpikes) > 2
                % First, second and third spike latencies, if more than
                % two spikes were found
                firstSpikeLat(newIndex,stim) = latencies(1);
                secondSpikeLat(newIndex,stim) = latencies(2);
                thirdSpikeLat(newIndex,stim) = latencies(3);
            end
        % End of stimului-iteration loop
        end
    % End of trial-iteration loop
    end


    % Store the latencies in a cell array [ms]
    combiLatSplit{file,1} = firstSpikeLat.*1000;
    combiLatSplit{file,2} = secondSpikeLat.*1000;
    combiLatSplit{file,3} = thirdSpikeLat.*1000;
    % Add the temperature for all trials
    combiLatSplit{file,4} = temperature_srt;
    % Add the sort index
    combiLatSplit{file,5} = sortIndex;
 
% End of file-iteration loop   
end

% End of function definition
end
