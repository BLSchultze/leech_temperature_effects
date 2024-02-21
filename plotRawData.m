%% Visualizations based on the raw recording traces
% 
% 1. Raw trace details for all stimulus pairs over all trials for one file
% 2. Response to single electrical and tactile pulse(s) over all trials
% 3. Individual trials of single files for visual analysis
% 4. Stimulation protocol (whole protocol + different details)
% 
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 15.03.2023 
% -------------------------------------------------------------------------

%% Preparations (common for all plots)
% Give the group to analyze ('Cold' or 'RoomTemp')
group = 'Cold';

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
                0.9     0.6     0.0];
% Define position values for new figure windows
masterPosition = [167,77,1161,702];

% Second display (student office)
% masterPosition = [2203.4,14.6,1160.8,702.4];

%% 1. Raw trace details for all stimulus pairs over all trials for one file

% Give the file number
file = 16;

% Load the dataset which should be displayed in detail
load(listDir(file).folder+"\"+listDir(file).name,'tacData','metadata')

% Define a window to cut around each stimulus pair
searchMin = -0.032; searchMin_ind = searchMin*metadata.samplingRate;
searchMax = 0.15; searchMax_ind = searchMax*metadata.samplingRate;

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
% electrical), time and indices
elec = elecOnsets(~ismember(elecOnsets,onlyElec));
tac = tacOnsets(~ismember(tacOnsets,onlyTac));
elec_ind = elecOnsets(~ismember(elecOnsets,onlyElec))* ...
    metadata.samplingRate;
tac_ind = tacOnsets(~ismember(tacOnsets,onlyTac))*metadata.samplingRate;

% Calculate the time difference between both stimuli in each pair
stimTimeDiff = (elec - tac)*1000;   % in ms

% Create a color gradient with a suitable number of colors
colors = colorGradient([0.2,0.6,0.9],[1,0.1,0.3],size(tacData,1));

% Set up a new figure window with an adaptive tiled layout
figure('WindowState','maximized')
tiledlayout('flow')

% Iterate through all stimulus pairs
for stim = 1:length(elec)
    % New tile for each stimulus pair
    nexttile
    hold on

    % Iterate through all trials 
    for trial = 1:size(tacData,1)
        % Create a suitable time vector
        timeVector = tac(stim)+searchMin:1/metadata.samplingRate: ...
            tac(stim)+searchMax; 
        % Plot the detail of the trace
        plot(timeVector, tacData{trial,2}( ...
            int64(tac_ind(stim)+searchMin_ind): ...
            int64(tac_ind(stim)+searchMax_ind)),'Color',colors(trial,:))
    end

    hold off
    % Add vertical lines indicating the time of tactile and electrical
    % stimuli
    xlTac = xline(tac(stim),'Color',[0.7,0.2,0.9]);
    xlElec = xline(elec(stim),'Color',[0.2,0.6,0.9]);

    % Label the subplot
    title("Stim offset (\Deltat) "+stimTimeDiff(stim)+" ms")
    xlabel("Time [s]")
    ylabel("Membrane potential [mV]")

    % Adapt the x axis limits for a uniform scaling of all subplots
    xlim([tac(stim)+searchMin,tac(stim)+searchMax])
end

% Add one legend for the tiled plots
legend([xlTac,xlElec],"Tactile stimulus","Electrical stimulus", ...
    'Position',[0.83,0.21,0.09,0.04])

%% 2. Response to single electrical and tactile pulse(s) over all trials

% Give the file number                                  %% ACCOUNT FOR ASCENDING AND DESCENDING TEMPERATURE
file = 11;
stimuli = 2:4;
trials = 1;

% Load the dataset which should be displayed in detail
load(listDir(file).folder+"\"+listDir(file).name,'tacData','metadata')
load("Cooled/tempSelect.mat")
tacData = tacData(tempSelect(file,1):tempSelect(file,2),:);

% Reset collection variable
traceDetail = [];

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

onlyTacInd = onlyTac .* metadata.samplingRate;
onlyElecInd = onlyElec .* metadata.samplingRate;

% Create a color gradient with a suitable number of colors
colors = colorGradient([0.2,0.6,0.9],[1,0.1,0.3],size(tacData,1));

% --- For an electrical pulse: 

% Define a window to cut around each stimulus pair
searchMin = -0.05; searchMin_ind = searchMin*metadata.samplingRate;
searchMax = 0.3; searchMax_ind = searchMax*metadata.samplingRate;

% Set up a new figure window with an adaptive tiled layout
figure('WindowState','maximized', ...
    'Name',"Stacked raw traces - Dataset "+file)
tiledlayout('flow')

% Iterate through all stimulus pairs
for stim = stimuli
    % New tile for each stimulus pair
    nexttile
    hold on

    % Iterate through all trials 
    for trial = 1:trials:size(tacData,1)
        % Create a suitable time vector
        timeVector = onlyElec(stim)+searchMin:1/metadata.samplingRate: ...
            onlyElec(stim)+searchMax; 
        % Plot the detail of the trace
        traceDetail(:,trial) = tacData{trial,2}(...
            int64(onlyElecInd(stim)+searchMin_ind): ...
            int64(onlyElecInd(stim)+searchMax_ind));
        plot(timeVector, traceDetail(:,trial),'Color',colors(trial,:), ...
            'LineWidth',2)
    end

    hold off
    % Add vertical lines indicating the time of tactile and electrical
    % stimuli
    xlTac = xline(onlyElec(stim),'Color',masterColors(3,:),'LineWidth',2);

    % Label the subplot
    title("Single electrical pulse - No. "+stim)
    xlabel("Time [s]")
    ylabel("Membrane potential [mV]")

    % Set limits for the x-axis
    xlim([min(timeVector),max(timeVector)])
end

% Add a suitable color bar
colormap(colors)
clb = colorbar('Ticks',[0,1],'TickLabels', ...
    {tacData{1,3},tacData{end,3}});
clb.Label.String = "Temperature [°C]";

% Catch axes handle
ax_traceDetail(1) = gca;

% --- For a tactile pulse: 

% Reset collection variable
traceDetail = [];

% Define a window to cut around each stimulus pair
searchMin = -0.05; searchMin_ind = searchMin*metadata.samplingRate;
searchMax = 0.3; searchMax_ind = searchMax*metadata.samplingRate;

% Set up a new figure window with an adaptive tiled layout
figure('WindowState','maximized', ...
    'Name',"Stacked raw traces - Dataset "+file)
tiledlayout('flow')

% Iterate through all stimulus pairs
for stim = stimuli
    % New tile for each stimulus pair
    nexttile
    hold on

    % Iterate through all trials 
    for trial = 1:trials:size(tacData,1)
        % Create a suitable time vector
        timeVector = onlyTac(stim)+searchMin:1/metadata.samplingRate: ...
            onlyTac(stim)+searchMax; 
        % Plot the detail of the trace
        traceDetail(:,trial) = tacData{trial,2}(...
            int64(onlyTacInd(stim)+searchMin_ind): ...
            int64(onlyTacInd(stim)+searchMax_ind));
        plot(timeVector, traceDetail(:,trial),'Color',colors(trial,:), ...
            'LineWidth',2)
    end

    hold off
    % Add vertical lines indicating the time of tactile and electrical
    % stimuli
    xlTac = xline(onlyTac(stim),'Color',masterColors(3,:),'LineWidth',2);

    % Label the subplot
    title("Single tactile pulse - No. "+stim)
    xlabel("Time [s]")
    ylabel("Membrane potential [mV]")

    % Set limits for the x-axis
    xlim([min(timeVector),max(timeVector)])
end

% Add a suitable color bar
colormap(colors)
clb = colorbar('Ticks',[0,1],'TickLabels', ...
    {tacData{1,3},tacData{end,3}});
clb.Label.String = "Temperature [°C]";

% Catch axes handle
ax_traceDetail(2) = gca;

%% 3. Individual trials of single files for visual analysis

% Give the file number
file = 2;

% Load the dataset which should be displayed in detail
load(strcat(listDir(file).folder,"\",listDir(file).name))
% Reduce the recording data according to the 'tempSelect' variable
% tacData = tacData(tempSelect(file,1):tempSelect(file,2),:);

% Give the trial numbers 
trials = 11:-1:8;

% Display the comments in the console
disp(metadata.comments)

% Iterate through all trials
for trial = trials
    % Create a new figure window
    figure('Position',masterPosition) 
    
    % 1st subplot: stimulation protocols
    subplot(2,1,1)
    plot(metadata.timeVector,tacData{trial,1})
    % Label the axes
    xlabel("Time [s]"); ylabel("Current [nA]")
    % 2nd subplot: recorded voltage trace
    subplot(2,1,2)
    plot(metadata.timeVector,tacData{trial,2})
    % Label the axes
    xlabel("Time [s]"); ylabel("Membrane potential  [mV]")
    % Combined title for the plots
    sgtitle("Trial "+trial)
end


%% 4. Stimulation protocol (whole protocol + different details)
% Please note for the tactile stimulation: We deliver a voltage to the
% stimulation device, which scales this voltage with 50 mN/V. This is then
% the stimulation force requested by the stimulation protocol. However,
% since the stimulation is very short here, the device cannot produce a
% step stimulus, but rather a steep stimulus with a little overshoot in
% amplitude.

% Load a trace with the actual tactile stimulation (recorded)
load("tac_stim_example.mat")

% Give the file number
file = 3;
% Give the trial number 
trials = 1;

% Load the dataset which should be displayed
load(strcat(listDir(file).folder,"\",listDir(file).name))
% Reduce the recording data according to the 'tempSelect' variable
tacData = tacData(tempSelect(file,1):tempSelect(file,2),:);

% --- Whole stimulation protocol
% Set up a new figure window
figure('Position',[239.4,202.6,988,557.6])

% 1st y-axis: electrical stimulus over time
yyaxis('left')
plot(metadata.timeVector,tacData{1, 1}(:,1))
% Adjust the y-axis limits and label the axis
ylim([-1.1,1.1])
ylabel("Current [nA]")
% Catch axes handle
ax_stimProt(1) = gca;

% 2nd y-axis: tactile stimulus over time
yyaxis('right')
plot(metadata.timeVector,tacData{1,1}(:,2)*50)
% Adjust the y-axis limits and label the axis
ylim([-11,11])
ylabel("Force [mN]")

% Label the x-axis and create a legend
xlabel("Time [s]")
legend("Electrical stimulation","Tactile stimulation")

% Catch axes handle
ax_stimProt(2) = gca;

% --- Detail: seconds 2-5 of the stimulation protocol
% Set up a new figure window
figure('Position',[239.4,202.6,988,557.6])

% 1st y-axis: electrical stimulus over time
yyaxis('left')
plot(metadata.timeVector,tacData{1, 1}(:,1))
% Adjust the axes limits and label the y-axis
ylim([-0.1,1.5])
xlim([1.8,5.3])
ylabel("Current [nA]")
% Catch axes handle
ax_stimProt(3) = gca;

% 2nd y-axis: tactile stimulus over time
yyaxis('right')
plot(metadata.timeVector,tacData{1,1}(:,2)*50)
% Adjust the y-axis limits and label the axis
ylim([-1,15])
ylabel("Force [mN]")

% Label the x-axis and create a legend
xlabel("Time [s]")
legend("Electrical stimulation","Tactile stimulation")

% Catch axes handle
ax_stimProt(4) = gca;

% --- Detail: seconds 3.96-4.01 of the stimulation protocol
% Set up a new figure window
figure('Position',[239.4,202.6,988,557.6])

% 1st y-axis: tactile stimulus over time
yyaxis('left')
plot(metadata.timeVector,tacData{1, 1}(:,1))
% Adjust the axes limits and label the y-axis
ylim([-0.1,1.1])
xlim([3.96,4.01])
ylabel("Current [nA]")
% Catch axes handle
ax_stimProt(5) = gca;

% 2nd y-axis: electrical stimulus over time
yyaxis('right')
plot(metadata.timeVector,tacData{1,1}(:,2)*50)
% Adjust the y-axis limits and label the axis
ylim([-1,11])
ylabel("Force [mN]")

% Label the x-axis and create a legend
xlabel("Time [s]")
legend("Electrical stimulation","Tactile stimulation")

% Catch axes handle
ax_stimProt(6) = gca;

%% Plot the recorded (actual) stimulation

% Set up a new figure window
figure('Position',[239.4,202.6,988,557.6])

% 1st y-axis: tactile stimulus over time
yyaxis('left')
plot(stim_x(stim_x > 3.965 & stim_x < 3.985), ...
    elec_stim_y(stim_x > 3.965 & stim_x < 3.985))
% Adjust the axes limits and label the y-axis
ylim([-0.1,1.1])
xlim([3.96,4.04])
ylabel("Current [nA]")
% Catch axes handle
ax_stimProt_r(1) = gca;

% 2nd y-axis: electrical stimulus over time
yyaxis('right')
plot(stim_x(stim_x > 3.995 & stim_x < 4.035), ...
    tac_stim_y(stim_x > 3.995 & stim_x < 4.035))
% Adjust the y-axis limits and label the axis
ylim([-1,11])
ylabel("Force [mN]")

% Label the x-axis and create a legend
xlabel("Time [s]")
legend("Electrical stimulation","Tactile stimulation")

% Catch axes handle
ax_stimProt_r(2) = gca;
