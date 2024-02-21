%% Visualize the results for the collide stimuli
% This script creates graphical representations for the results of the
% collide stimuli analysis. The script takes either the results for the
% cold or the room temperature experiments. The plots that are created are
% listed below.
% -----------------
% Table of contents
% -----------------
% 1. Latencies as a function of the stimulus time difference
% 2. Temperature raster plots per stimulus pair
% 3. Spike count difference over stimulus time difference
% 4. Expected and actual spike count for each stimulus pair
% 5. Spike count difference for each stimulus pair
% 
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 28.02.2023
% -------------------------------------------------------------------------

%% Preparations
% Clear the workspace
clear
% Change working directory
cd(['C:\Users\schul\Nextcloud\Studium\2. Master-Semester\' ...
    'ResearchModule_Leech\'])

% Decide for which group the results should be plotted
group = "Cold";     % "RoomTemp" or "Cold"

if isequal(group,"Cold")
    % Load the file with the analysis results (cold data)
    load("analysisResultsCold.mat")
elseif isequal(group,"RoomTemp")
    % Load the file with the analysis results (room temperature data)
    load("analysisResultsRoomTemp.mat")
end

% Define a set of colors for a consistent visualization
masterColors = [0.0     0.6     0.9; ...
                0.4     0.75    0.0; ...
                0.9     0.6     0.0];
% Define position values for new figure windows
masterPosition = [167,77,1161,702];


%% 1. Latencies as a function of the stimulus time difference
% One figure for each dataset

% Iterate through all datasets
for dataset = 12%:size(combiLatSplit,1)

    colors1st = colorGradient(masterColors(1,:),[0.8,0.95,1], ...
        size(combiLatSplit{dataset,1},1));
    colors2nd = colorGradient(masterColors(2,:),[0.8,1,0.6], ...
        size(combiLatSplit{dataset,1},1));
    colors3rd = colorGradient(masterColors(3,:),[1,0.9,0.7], ...
        size(combiLatSplit{dataset,1},1));

    % Create new figure window and set hold on
    figure('Position',masterPosition,'Name',"Collide compact - " + ...
        "Dataset " + dataset)
    hold on

    % Iterate through all trials
    for trial = 1:size(combiLatSplit{dataset,1},1)
        % First spike latency as a function of the stimulus time difference
        scatter(combiLatSplit{dataset,1}(trial,:),stimTimeDiff,36, ...
            colors1st(trial,:),'filled','MarkerFaceAlpha',0.7)
        % Add the second spike latencies
        scatter(combiLatSplit{dataset,2}(trial,:),stimTimeDiff,36, ...
            colors2nd(trial,:),'filled','MarkerFaceAlpha',0.7) 
        % Add the thrid spike latencies
        scatter(combiLatSplit{dataset,3}(trial,:),stimTimeDiff,36, ...
            colors3rd(trial,:),'filled','MarkerFaceAlpha',0.7)

        % In the first iteration: catch the handles for the scatter plots
        if trial == 1
            % Catch the axes handle
            ax = gca;
            % Get the scatter plot handles
            scatterHandles = ax.Children;
        end
    end

    % Set hold off
    hold off
    % Insert horizontal lines at stimulus time differences for better
    % visual interpreation
    yline(stimTimeDiff,'Color',[0.8,0.8,0.8,0.3])
    
    % Label the plot
    ylabel("Stimuli timing (\Deltat, elec-tac) [ms]")
    xlabel("Spike timing relative to tactile stimulus [ms]")
    
    % Define reference lines for the expected electrical and tactile spikes
    % based on the single puls latencies
    elecRefLine = stimTimeDiff + singleStimResults{dataset,4}(1,1);
    tacRefLine = singleStimResults{dataset,3}(1,1);

    % Add the reference lines to the plot
    hold on;
    eRef = plot(elecRefLine,stimTimeDiff, ...
        'Color',[0,0.4,0.8],'LineStyle','--','LineWidth',1);
    tRef = xline(tacRefLine,'Color',[0.8,0.2,0.8],'LineStyle','--', ...
        'LineWidth',1);
    hold off

    % Add a suitable legend
    legend([scatterHandles([3,2,1],1); eRef; tRef], ...
        "First spike","Second spike","Third spike", ...
        "Expected spike (electrically induced)", ...
        "Expected spike (tactilely induced)", ...
        'Location','northwest')
end


%% 2. Temperature raster plots per stimulus pair

% Iterate through all datasets
for dataset = 1:size(combiLatSplit,1)
    % Create a new figure window with  a fluid tiledlayout
    figure('WindowState','maximized','Name',"Temp raster - Dataset " ...
        +dataset)
    tiledlayout('flow'); 
    %tiledlayout(2,4,'TileSpacing','compact')

    % Find the overall maximum latency (for axis scaling)
    m_trial = cellfun(@max, combiLatSplit(dataset,1:3), ...
        'UniformOutput',false);
    m_stim = cellfun(@max, m_trial);
    absMax = round(max(m_stim),0) + 2;      % round to integer and add 2

    % Initialize a stimulus counter variable
    stimInd = 0;

    % Iterate through all stimulus time differences
    for stim = 1:size(stimTimeDiff,1)
        stimInd = stimInd + 1;
        % Create next subplot and put hold on
        ax_tempRaster(stimInd) = nexttile;
        hold on
        % Iterate through all trials
        for trial = 1:size(combiLatSplit{dataset,1},1)
            % Plot the latencies (x) with the according temperature (y)
            plot(combiLatSplit{dataset,1}(trial,stim), ...
                combiLatSplit{dataset,4}(trial), ...
                '.','Color',masterColors(1,:),'MarkerSize',10)
            plot(combiLatSplit{dataset,2}(trial,stim), ...
                combiLatSplit{dataset,4}(trial), ...
                '.','Color',masterColors(2,:),'MarkerSize',10)
            plot(combiLatSplit{dataset,3}(trial,stim), ...
                combiLatSplit{dataset,4}(trial), ...
                '.','Color',masterColors(3,:),'MarkerSize',10) 

            % In the first iteration: catch the handles for the plots
            if trial == 1
                % Catch the axes handle
                ax = gca;
                % Get the plot handles
                plotHandles = ax.Children;
            end
        end
        
        % Add vertical lines to mark the stimulus onsets
        xElec = xline(stimTimeDiff(stim), ...
            'Color',[0.2,0.8,0.8], 'LineWidth',1);
        xTac = xline(0,'Color',[0.6,0.1,0.9],'LineWidth',1);

        % Set hold of for the current subplot
        hold off

        % Label the subplot
        xlabel("Latency [ms]")
        ylabel("Temperature [°C]")
        title("Stim offset (\Deltat) "+stimTimeDiff(stim)+" ms")

        % Adjust the x limits
        xlim([-32,absMax])
    end
end

% Add one suitable legend for all subplots
legend([xElec;xTac;plotHandles([3,2,1],1)], ...
    "Onset electrical stimulus","Onset tactile stimulus", ...
    "First spike", "Second spike", "Third spike", ...
    'Location',[0.78,0.15,0.11,0.1])
 


%% 3. Spike count difference over stimulus time difference
% --- all trials + median and mean 
%     Please note: mean and median are only useful for the room temperature 
%                  experiments ---

% Iterate through all datasets
for dataset = 10%:size(combiLatSplit,1)
    % Set up a new figure window and put hold on
    figure('Position',masterPosition,'Name',"Spike count diff - " + ...
        "Dataset "+dataset) 
    % Calculate difference between expected and actual spike count
    spcDiff = spcActual{dataset,1} - spcExpected{dataset,1};
    % Plot the spike count differnce over the stimulus time difference 
    plot(stimTimeDiff,spcDiff,'.-')

    % Add the median spike count difference
    hold on
    medPlt = plot(stimTimeDiff,median(spcDiff),'.-','LineWidth',2.5, ...
        'Color',masterColors(1,:));
    meanPlt = plot(stimTimeDiff,mean(spcDiff),'.-','LineWidth',2.5, ...
        'Color',masterColors(2,:));
    hold off

    % Label the plot
    xlabel("Stimuli time difference (elec-tac) [ms]")
    ylabel("Spike count difference (actual-expected)")
    title("Spike count difference for each stimulus pair and trial")
    legend([medPlt,meanPlt],"Median","Mean")
end


%% 4. Expected and actual spike count for each stimulus pair

% Iterate through all datasets
for dataset = 14%:size(combiLatSplit,1)
    % Set up a new figure window with a flexible tiled layout
    figure ('Position',masterPosition,'Name',"Spike count diff ~ " + ...
        "temp - Dataset " + dataset)
    tiledlayout('flow')

    % Iterate through all stimulus time differences
    for stim = 1:length(stimTimeDiff)
        % Activate the nex tile
        nexttile; hold on;
        % Plot the expected spike count
        pltExp = plot(combiLatSplit{dataset,4}, ...
            spcExpected{dataset,1}(:,stim), ...
            '.-','Color',masterColors(2,:));
        % Plot the actual spike count
        pltAct = plot(combiLatSplit{dataset,4}, ...
            spcActual{dataset,1}(:,stim), ...
            '.-','Color',masterColors(1,:));
        
        % Label the subplot
        xlabel("Temperature [°C]")
        ylabel("Spike count")
        title("Stim offset (\Deltat) "+stimTimeDiff(stim)+" ms")
    end
end

% Add a suitable legend
legend([pltAct(end),pltExp(end)],"Actual spike count", ...
    "Expected spike count",'Location',[0.78,0.21,0.13,0.04])


%% Spike count difference for each stimulus pair

% Iterate through all datasets
for dataset = 2%:size(combiLatSplit,1)
    % Set up a new figure window with a flexible tiled layout
    figure ('Position',masterPosition,'Name',"Spike count diff ~ " + ...
        "temp - Dataset " + dataset)
    tiledlayout('flow')

    % Calculate difference between expected and actual spike count
    spcDiff = spcActual{dataset,1} - spcExpected{dataset,1};

    % Iterate through all spike time differences
    for stim = 1:length(stimTimeDiff)
        % Activate the next tile and plot the spike count difference 
        nexttile
        plot(combiLatSplit{dataset,4},spcDiff(:,stim),'.-','Color', ...
            masterColors(1,:))

        % Label the subplot
        xlabel("Temperature [°C]")
        ylabel("Spike count diff")
        title("Stim offset (\Deltat) "+stimTimeDiff(stim)+" ms")
    end
end
