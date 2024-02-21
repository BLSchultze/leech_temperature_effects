%% Check the stability of the recordings
% The resting membrane potential and the input resistance are determined as
% two stability parameters. Changes over time and with temperature are
% analyzed as well as the effect on the spike times.
% 
% Table of contents
% -----------------
% 1. Calculation of resting membrane potential and input resistance
% 2. Single pulse latencies over RMP and input resistance (IR)
% 3. RMP and input resistance (IR) over the temperature
% 4. Check for pulses without a response
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 28.02.2023
% -------------------------------------------------------------------------

%% Preparations
% Give the group to analyze ('Cold' or 'RoomTemp')
group = 'RoomTemp';

% List all dataset files for the chosen group
if isequal(group,'Cold')
    listDir = dir("Cooled\tacData*.mat");
    % Load a file with the trials to select
    load('Cooled\tempSelect.mat','tempSelect')
elseif isequal(group,'RoomTemp')
    listDir = dir("RoomTemp\tacData*.mat");
    % Load a file with the trials to select
    load('RoomTemp\tempSelect.mat','tempSelect')
end

% Define a set of colors for a consistent visualization
masterColors = [0.0     0.6     0.9; ...
                0.4     0.75    0.0; ...
                0.75    0.0     0.35];
% Define position values for new figure windows
masterPosition = [167,77,1161,702];

%% 1 .Check the resting membrane potential (RMP) and the input resistance (IR) 
% Calculate RMP and IR, plot the results seperately for each file
qualityCheck = check_RMP_IR(1:length(listDir),group,0);

% Save the results in a file
save("rmp_ir_"+group,"qualityCheck")

%% 2. Single pulse latencies over RMP and input resistance (IR)

% Load the analysis results (single pulse latencies)
load("analysisResults"+group+".mat")

% Set up a new figure window with a tiled layout 
figure('Position',masterPosition)
tiledlayout(2,2,'TileSpacing','compact')
% Create all tiles, catch the handles, set hold on for each tile
tl1 = nexttile; hold(tl1,'on')
tl2 = nexttile; hold(tl2,'on')
tl3 = nexttile; hold(tl3,'on')
tl4 = nexttile; hold(tl4,'on')

% Iterate through all datasets
for dataset = 1:length(singleStimResults)
    % 1st tile: latency (single tactile pulses) over RMP
    plot(tl1,qualityCheck{dataset,1},singleStimResults{dataset,3},'.', ...
        'Color','k')
    % 3rd tile: latency (single tactile pulses) over IR
    plot(tl3,qualityCheck{dataset,3},singleStimResults{dataset,3},'.',...
        'Color','k')

    % 2nd tile: latency (single electrical pulses) over RMP
    plot(tl2,qualityCheck{dataset,1},singleStimResults{dataset,4},'.', ...
        'Color','k')
    % 4th tile: latency (single electrical pulses) over IR
    plot(tl4,qualityCheck{dataset,3},singleStimResults{dataset,4},'.',...
        'Color','k')
end

% Titles for the subplots
title(tl1,"Tactile stimulation"); title(tl2,"Electrical stimulation")
% Label the axes
xlabel([tl1,tl2],"RMP [mV]")
xlabel([tl3,tl4],"IR [M\Omega]")
ylabel([tl1,tl2,tl3,tl4],"Latency [ms]")


%% 3. RMP and input resistance (IR) over the temperature

% Set up a new figure window with a tiled layout
figure('Position',masterPosition)
tiledlayout(2,2,'TileSpacing','compact')
% Create the tiles, catch the handles, set hold on
tl1 = nexttile; hold(tl1,'on')
tl2 = nexttile; hold(tl2,'on')
tl3 = nexttile; hold(tl3,'on')
tl4 = nexttile; hold(tl4,'on')

% Iterate through all datasets
for dataset = 1:length(singleStimResults)
    % 1st tile: RMP over temperature
    plot(tl1,singleStimResults{dataset,7},qualityCheck{dataset,1},'.', ...
        'Color','k')
    % 2nd tile: IR over temperature
    plot(tl2,singleStimResults{dataset,7},qualityCheck{dataset,3},'.',...
        'Color','k')

    % Sort the original sort indices to reverse the sorting according to
    % temperature 
    [~,reSort] = sort(singleStimResults{dataset,10});
    % 3rd tile: RMP over original trial index
    plot(tl3,qualityCheck{dataset,1}(reSort),'.-','Color','k')
    % 4th tile: IR over original trial index
    plot(tl4,qualityCheck{dataset,3}(reSort),'.-','Color','k')
end

% Titles for the subplots 
title(tl1,"RMP over temperature"); title(tl2,"IR over temperature")
title(tl3,"RMP over original trial index")
title(tl4,"IR over original trial index")
% Label the axes
xlabel([tl1,tl2],"Temperature [°C]")
xlabel([tl3,tl4],"Trial index")
ylabel([tl1,tl3],"RMP [mV]")
ylabel([tl2,tl4],"IR [M\Omega]")


%% 4. Check for pulses without a response
% Load the analysis results of the given group
load("analysisResults"+group+".mat")

% --- For the electrical pulses ---
% Iterate through all datasets
for dataset = 1:length(listDir)
    % Check if the data are in original order (1) or not (0)
    sortedLog = issorted(singleStimResults{dataset,10});
    % Find all single stimuli for which no reponse was found
    [rows,cols] = find(singleStimResults{dataset,6} == 0);
    % If there were missing responses found ...
    if ~isempty(rows)
        % Display a messange in the console stating all trials and pulses
        % for which there was no response
        fprintf("Dataset %g - spikes missing (e pulses) " + ...
            "[stimulus],[trial]: %s, %s (Sorted: %g)\n", ...
            dataset,mat2str(rows'),mat2str(cols'),sortedLog);
    end
end
% Print separator
disp("-------------------")

% --- For the tactile pulses ---
% Iterate through all datasets
for dataset = 1:length(listDir)
    % Check if data are in original order (1) or not (0)
    sortedLog = issorted(singleStimResults{dataset,10});
    % Check for missing spikes
    [rows,cols] = find(singleStimResults{dataset,5} == 0);
    % If spikes were found to be missing ...
    if ~isempty(rows)
        % Print an adequate message to the console
        fprintf("Dataset %g - spikes missing (t pulses) " + ...
            "[rows],[cols]: %s, %s (Sorted: %g)\n", ...
            dataset,mat2str(rows'),mat2str(cols'),sortedLog);
    end
end
