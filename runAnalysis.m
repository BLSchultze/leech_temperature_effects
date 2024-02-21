%% Perform the main analysis (latency calculations)
% This script performs the latency and spike count analysis for the data
% recorded with the elec-tac collide protocol (collide_v4, Steffen). 
% The single tactile/electrical pulses and the collide pulses are analyzed
% separately by different functions. Latency and spike count differences
% are calculated within this script. The analysis can be performed for the
% cold and the room temperature experiments separately. 
% The following results are saved in a .mat file:
%   singleStimResults   - analysis results for the single
%                         electrical/tactile pulses (latencies, counts, Q10
%                         values)
%   combiLatSplit       - spike latencies for the collide part, splitted in 
%                         1st, 2nd, 3rd
%   combiStimLat        - spike latencies for all spikes in the collide
%                         part
%   relLat_Combi_1stTrial   - combiLatSplit latencies, given as latency
%                             differences to the first trial (lowest
%                             temperature)
%   relLat_elec         - single pulse latencies (elec/tac pulses), given 
%   relLat_tac            as latency differences to the trial with the 
%                         highest temperature, temperatures also given as
%                         differences to highest temperature
%   spcActual           - spike count observed during the collide stimuli,
%                         seperate for each stimulus pair
%   spcExpected         - spike count expected for each collide stimulus
%                         pair based on the sum of the spike counts for the
%                         following pair of a single tactile and electrical 
%                         puls
%   stimTimeDiff        - vector with the time difference between both
%                         stimuli in each collide stimulus pair
% 
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 02.08.2023
% -------------------------------------------------------------------------

%% Run analysis for the cooled experiments

% Clear the workspace
clear

group = "Cold";

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
files = 1:length(dir(dataPath+"tacData*.mat"));

% Calculate the latencies for the stimulus pairs
[combiStimLat,combiLatSplit,stimTimeDiff] = measureCombiStim(files,group);

% Calculate the latencies for the single stimuli (electrical and tactile)
singleStimResults = measureSingleStim(files,group,0);

%% Single pulse analysis
% Allocate variables
relLat_elec = cell(size(singleStimResults,1),2);
relLat_tac = cell(size(singleStimResults,1),2);

% Calculate relative latencies (all latencies relative to the latencies in
% the trial with the highest temperature)
for cellN = 1:size(singleStimResults,1)
    % For electrical stimulation
    relLat_elec{cellN,1} = relativize(singleStimResults{cellN,2}(2,:));
    relLat_elec{cellN,2} = relativize(singleStimResults{cellN,7})';
    % For tactile stimulation
    relLat_tac{cellN,1} = relativize(singleStimResults{cellN,1}(2,:));
    relLat_tac{cellN,2} = relativize(singleStimResults{cellN,7})';
end

%% Analysis of combined stimulations
% Pre-define a cell array to collect the results
relLat_Combi_1stTrial = cell(length(files),3);

% Calculate all latencies relative to the first trial
for file = files
    relLat_Combi_1stTrial{file,1} = combiLatSplit{file,1} ...
        - combiLatSplit{file,1}(1,:);
    relLat_Combi_1stTrial{file,2} = combiLatSplit{file,2} ...
        - combiLatSplit{file,2}(1,:);
    relLat_Combi_1stTrial{file,3} = combiLatSplit{file,3} ...
        - combiLatSplit{file,3}(1,:);
end

%% Calculate expected and actual spike counts
% Pre-allocate two cell arrays to collect the spike counts
spcExpected = cell(length(singleStimResults),1);
spcActual = spcExpected;

% Iterate through all datasets
for dataset = 1:length(singleStimResults)
    % Pre-allocate two matrices to temporarily collect the spike counts
    spcExpected_temp = zeros(length(singleStimResults{dataset,5}), ...
        length(stimTimeDiff));
    spcActual_temp = spcExpected_temp;
    % For the first 6 stimulus pairs: combined spike count for second 
    % tactile and electrical pulse as a reference
    spcExpected_temp(:,1:6) = ...
        repmat((singleStimResults{dataset,5}(2,:) + ...
        singleStimResults{dataset,6}(2,:))',1,6);
    % For stimulus pairs 7-12: combined spike count for the third tactile
    % and electrical pulse as a reference
    spcExpected_temp(:,7:12) = ...
        repmat((singleStimResults{dataset,5}(3,:) + ...
        singleStimResults{dataset,6}(3,:))',1,6);
    % For stimulus pairs 13-19: combined spike count for the fourth tactile
    % and electrical pulse as a reference
    spcExpected_temp(:,13:19) = ...
        repmat((singleStimResults{dataset,5}(4,:) + ...
        singleStimResults{dataset,6}(4,:))',1,7);

    % Iterate through all stimulus time differences 
    for stim = 1:size(combiStimLat,2)
        % Calculate the actual spike count for all trials
        spcActual_temp(:,stim) = ...
            cellfun(@length,combiStimLat{dataset,stim});
    end

    % Save the spike counts for the current dataset
    spcExpected{dataset,1} = spcExpected_temp;
    spcActual{dataset,1} = spcActual_temp;
end

%% Save the result
% Save the remaining variables
save("analysisResults"+group+".mat", ...
    "relLat_Combi_1stTrial", ...
    "relLat_tac","relLat_elec", ...
    "combiLatSplit","combiStimLat","stimTimeDiff","singleStimResults", ...
    "spcActual","spcExpected")

clear 