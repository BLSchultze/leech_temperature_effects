%% Script to create the figures for the written report
% Figures are saved as eps files. 
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 27.02.2024
% -------------------------------------------------------------------------

% Plot the stimulation protocol (complete and detail)
plot_stimulation()

%% Plot the responses to the four separate stimuli for comparison
singleStimTests("RoomTemp")
singleStimTests("Cold")

%% Plot example traces of a cooling experiment

% Load the data for Figure 5
load("data_figure5_paper_13_07.mat")
% Selecte how many trials to plot (e.g., 2 means every second trial)
trialStep = 3;
% Create linear color gradient between the two given colors (45 steps)
colors = colorGradient([0.2,0.6,0.9],[1,0.1,0.3],45);

% Set up figure
fig4 = figure('Position',[257,172,1300,589]);
tiledlayout(2,1, 'TileSpacing', 'compact')

% Create tile
tileA1 = nexttile;

% --- Trace detail for electrical stimulation
% Set start and stop to only plot a part of the trace detail
% start = 400; stop = 2000;
start = 400; stop = 1950;

hold on
% Iterate through all trials 
for trial = 1:trialStep:size(traceDetail_elec,2)
     
    % Plot the detail of the trace
    plot(timeVector_elecDetail(start:stop), ...
        traceDetail_elec(start:stop,trial),'Color',colors(trial,:), ...
        'LineWidth',2)
end

hold off

% Add vertical lines indicating the time of electrical stimulus
xline(10,'Color',[0.7,0.2,0.9],'LineWidth',2)

% Label the subplot
xlabel("Time [s]")
ylabel(["Membrane","potential [mV]"])

% Adjust the axes limits
xlim([timeVector_elecDetail(start),timeVector_elecDetail(stop)])
ylim([-65,40]); yticks(-40:40:40)
% Catch axes handle
ax(1) = gca; bold(ax(1))

% Add a suitable color bar
colormap(colors)
clb = colorbar('Ticks',[0,1],'TickLabels',["9.2","15.8"]);
clb.Label.String = "Temperature [°C]";
%clb.Label.VerticalAlignment = 'middle';

% Place a title
title(tileA1,"A")
tileA1.TitleHorizontalAlignment = 'left'; 
tileA1.TitleFontSizeMultiplier = 1.5; tileA1.Title.Units = 'normalized';
tileA1.Title.Position = [-0.08,1.03,0];

% --- For a tactile pulse:
% Create a tile (1x4 space, 2nd row)
tileA2 = nexttile;

% Set start and stop to only plot a part of the trace detail
start = 400; stop = 1950;

hold on 
% Iterate through all trials 
for trial = 1:trialStep:size(traceDetail_tac,2)
    
    plot(timeVector_tacDetail(start:stop), ...
        traceDetail_tac(start:stop,trial),'Color', ...
        colors(trial,:),'LineWidth',2)
end

hold off

% Add vertical lines indicating the time of tactile and electrical
% stimuli
xlTac = xline(11,'Color',[0.7,0.2,0.9],'LineWidth',2);

% Label the subplot
xlabel("Time [s]")
ylabel(["Membrane","potential [mV]"])
title(tileA2, "B")
tileA2.TitleHorizontalAlignment = 'left';
tileA2.TitleFontSizeMultiplier = 1.5; tileA2.Title.Units = 'normalized';
tileA2.Title.Position = [-0.08,1.03,0];

% Adjust the axes limits
xlim([timeVector_tacDetail(start),timeVector_tacDetail(stop)])
ylim([-65,40]); yticks(-40:40:40)
% Catch axes handle
ax(2) = gca; bold(ax(2))

% Add a suitable color bar
colormap(colors)
clb = colorbar('Ticks',[0,1],'TickLabels', ["9.2","15.8"]);
clb.Label.String = "Temperature [°C]";
%clb.Label.VerticalAlignment = 'middle';

% Save figure 
exportgraphics(fig4,'Report/figures/example_traces_cold.eps', ...
    'ContentType','vector','BackgroundColor','white')


%% Figures for the separate stimuli
close all   % close all previous figures

% Extend color order
COrder = [  0	    0.447	0.741;
            0.851	0.325	0.0980;
            0.929	0.694	0.125;
            0.494	0.184	0.557;
            0.467	0.675	0.188;
            0.302	0.745	0.933;
            0.635	0.0780	0.184;
            0.604	0.859	0.008;
            0.851	0	    0.525;
            0	    0.757	0.812;
            0.502	0.502	0.502;
            0.910	0.910	0;
            0.525	0	    0.812;
            0.871	0	    0];

% Plot separate stim results
[fig_lat_temp, fig_lat_trial, fig_rel_lat_temp, ...
    fig_lat_trial_idx, fig_q10] = plotResultsSingleStim('Cold');
% Close unsued figures 
close([3:5,7:9,11:13,16:18])


% Modify latency temperature plot with stim amplitude highlighted
axFigE = fig_lat_temp(6).Children;
bold(axFigE(2))
axFigE(2).Title.String = "";


% Modify plot appearance
bold(fig_lat_temp(1).Children); bold(fig_lat_temp(2).Children)
fig_lat_temp(1).Children.Title.String = ""; 
fig_lat_temp(2).Children.Title.String = "";
% Change color for two datasets to enhace visibility
fig_lat_temp(1).Children.ColorOrder = COrder;
fig_lat_temp(2).Children.ColorOrder = COrder;

% Add dashed line indicating missing values
hold(fig_lat_temp(2).Children, 'on')
plot([15.7,17.4],[7.39,6.8],'--','Color',[0.00,0.76,0.81],'LineWidth',2)
hold(fig_lat_temp(2).Children, 'off')


% Single pulse latencies relative to warmest trial
bold(fig_rel_lat_temp.Children(1))
bold(fig_rel_lat_temp.Children(2), 10)
fig_rel_lat_temp.Children(2).Title.String = "";
fig_rel_lat_temp.Children(1).Box = "on";


% Modify the Q10 boxplot
fig_q10.Position = [488,440,544,350];
fig_q10.Children.LineWidth = 2; 
fig_q10.Children.FontSize = 14;
fig_q10.Children.Box = 'off';
fig_q10.Children.YLim = [-0.1, 1.9];

bxplt = fig_q10.Children.Children;
set(bxplt.Children(1), 'MarkerEdgeColor', [0, 0.6, 0.9], 'LineWidth', 2)
set(bxplt.Children(5), 'Color', [0, 0.6, 0.9], 'LineWidth', 2)
set(bxplt.Children(2), 'MarkerEdgeColor', [0.4, 0.75, 0], 'LineWidth', 2)
set(bxplt.Children(6), 'Color', [0.4, 0.75, 0], 'LineWidth', 2)
set(bxplt.Children(3:4), 'Color', [1, 0.5, 0.2], 'LineWidth', 2)
set(bxplt.Children(7:14), 'LineWidth', 2, 'Color', 'k')


% Modify spike count over temperature plot
spkcount_temp = findobj( 'Type', 'Figure', 'Name', 'SpikeCountTemp');
tileE = spkcount_temp.Children.Children(2);
tileT = spkcount_temp.Children.Children(1);
bold(tileE); bold(tileT)
% Adjust layouot
tileE.YLim = [-0.1,2.1];
tileT.YLim = [-0.1,2.1];
tileE.YTick = [0,1,2];
tileT.YTick = [0,1,2];
tileE.ColorOrder = COrder;
tileT.ColorOrder = COrder;
tileE.TitleHorizontalAlignment = 'left';
tileT.TitleHorizontalAlignment = 'left';
tileE.TitleFontSizeMultiplier = 1.5; tileE.Title.Units = 'normalized';
tileT.TitleFontSizeMultiplier = 1.5; tileT.Title.Units = 'normalized';


% Save plots
exportgraphics(fig_lat_temp(1), 'Report/figures/lat_temp_t.eps', ...
    'ContentType','vector','BackgroundColor','white')
exportgraphics(fig_lat_temp(2), 'Report/figures/lat_temp_e.eps', ...
    'ContentType','vector','BackgroundColor','white')
exportgraphics(fig_q10, 'Report/figures/q10_boxplot.eps', ...
    'ContentType','vector','BackgroundColor','white')
exportgraphics(fig_rel_lat_temp, 'Report/figures/rel_lat_rel_temp.eps', ...
    'ContentType', 'vector', 'BackgroundColor', 'white')
exportgraphics(spkcount_temp, 'Report/figures/spkcout_temp.eps', ...
    'ContentType','vector','BackgroundColor','white')
exportgraphics(fig_lat_temp(6), 'Report/figures/lat_temp_stim.eps', ...
    'ContentType','vector','BackgroundColor','white')


%% Temperature raster plot (combined stimuli)
close all    % close all previous figures
% Define colors
elec_col = [143,21,75]/255; elec_col_light = [200,146,170]/255;
tac_col = [67,138,63]/255; tac_col_light = [135,200,135]/255;

% Load data
load("analysisResultsCold.mat")
colors = colorGradient([0.2,0.6,0.9],[1,0.1,0.3],17);

dataset = 8;
% Create a new figure window 
figure('Position',[0.2,49.8,930,732.8])
tiledlayout(5,4,'TileSpacing','compact', 'Padding','compact')

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
    % bold(ax_tempRaster(stimInd))
    ax_tempRaster(stimInd).XAxis.LineWidth = 1.5;
    ax_tempRaster(stimInd).YAxis.LineWidth = 1.5;

    hold on
    % Iterate through all trials
    for trial = 1:size(combiLatSplit{dataset,1},1)
        % Plot the latencies (x) with the according temperature (y)
        plot(combiLatSplit{dataset,1}(trial,stim), ...
            combiLatSplit{dataset,4}(trial), ...
            '.','Color',colors(trial,:),'MarkerSize',17)
        plot(combiLatSplit{dataset,2}(trial,stim), ...
            combiLatSplit{dataset,4}(trial), ...
            '.','Color',colors(trial,:),'MarkerSize',17)
        plot(combiLatSplit{dataset,3}(trial,stim), ...
            combiLatSplit{dataset,4}(trial), ...
            '.','Color',colors(trial,:),'MarkerSize',17) 

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
        'Color',elec_col, 'LineWidth',3);
    xTac = xline(0,'Color',tac_col, 'LineWidth',3);

    % Set hold of for the current subplot
    hold off

    % Label the subplot
    title("\Deltat "+stimTimeDiff(stim)+" ms")

    % Adjust the x limits
    xlim([-35,absMax])
    ylim([10,18])
end

% Add a suitable color bar
colormap(colors)
clb = colorbar('Ticks',[0,1],'TickLabels',["10.8","17.4"], 'Position', ...
    [0.8669    0.5234    0.0165    0.3387]);
clb.Label.String = "Temperature [°C]";
clb.Label.Position = [2.5,0.5,0];
clb.Label.FontSize = 12;
clb.FontSize = 12;
clb.Layout.Tile = 'east';

% Add axes labels
ax_tempRaster(18).XLabel.String = 'Time after tactile stimulus [ms]';
ax_tempRaster(18).XLabel.FontSize = 14;
ax_tempRaster(18).XLabel.Position = [80    8.9   -1];
ax_tempRaster(9).YLabel.String = 'Temperature [°C]';
ax_tempRaster(9).YLabel.FontSize = 14;

legend([xElec,xTac],{'Onset electrical pulse','Onset tactile pulse'}, ...
    'Position',[0.71,0.15,0.16,0.04])

% Save figure
exportgraphics(gcf, 'Report/figures/temp_raster.eps', ...
   'BackgroundColor','white','ContentType','vector')


%% Plots for the control experiments
close all     % close all previous figures

% Extend color order
COrder = [  0	    0.447	0.741;
            0.851	0.325	0.0980;
            0.929	0.694	0.125;
            0.494	0.184	0.557;
            0.467	0.675	0.188;
            0.302	0.745	0.933;
            0.635	0.0780	0.184;
            0.604	0.859	0.008;
            0.851	0	    0.525;
            0	    0.757	0.812;
            0.502	0.502	0.502;
            0.910	0.910	0;
            0.525	0	    0.812;
            0.871	0	    0];

% Plot separate stim results, discard the unwanted ones
[~,lat_trial,~,~] = plotResultsSingleStim('RoomTemp');
close([1:4,7:11,13,14])

% Adjust the appearance of the latency ~ trial plots
lat_trial(1).Children.Title.String = '';
lat_trial(2).Children.Title.String = '';
bold(lat_trial(1).Children); bold(lat_trial(2).Children)
lat_trial(1).Children.ColorOrder = COrder;
lat_trial(2).Children.ColorOrder = COrder;

% Modify spike count over trial plot (room temp)
spkcount_trial = findobj( 'Type', 'Figure', 'Name', 'SpikeCountTrial');
tileE = spkcount_trial.Children.Children(2);
tileT = spkcount_trial.Children.Children(1);
% Adjust layout
bold(tileE); bold(tileT);
tileE.YLim = [-0.1,2.1];
tileT.YLim = [-0.1,2.1];
tileE.YTick = [0,1,2];
tileT.YTick = [0,1,2];
tileE.ColorOrder = COrder;
tileT.ColorOrder = COrder;
tileE.TitleHorizontalAlignment = 'left';
tileT.TitleHorizontalAlignment = 'left';
tileE.TitleFontSizeMultiplier = 1.5; tileE.Title.Units = 'normalized';
tileT.TitleFontSizeMultiplier = 1.5; tileT.Title.Units = 'normalized';

% Save plots
exportgraphics(spkcount_trial, 'Report/figures/spkcout_trial.eps', ...
    'ContentType','vector','BackgroundColor','white')
exportgraphics(lat_trial(2), 'Report/figures/lat_trial_e.eps', ...
    'ContentType','vector','BackgroundColor','white')
exportgraphics(lat_trial(1), 'Report/figures/lat_trial_t.eps', ...
    'ContentType','vector','BackgroundColor','white')
