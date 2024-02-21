function singleStimTests(group)
%% Test and plot the first spike latencies in response to separate stimuli
% 
% -----
% Input
%   group       giving the experimental group, 'RoomTemp' or 'Cold'
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       [last modified: 01.08.2023]
% -------------------------------------------------------------------------

% Load the room temperature results
load(sprintf("analysisResults%s.mat",group), 'singleStimResults')

% Pre-allocate variables for result-saving
sy = size(singleStimResults,1);
second_spikes_t = zeros(1,sy); second_spikes_e = zeros(1,sy);
third_spikes_t = zeros(1,sy); third_spikes_e = zeros(1,sy);
fourth_spikes_t = zeros(1,sy); fourth_spikes_e = zeros(1,sy);

% Iterate over all datasets
for cell = 1:size(singleStimResults,1)
    % Relativize all latencies to the latency of the first separate
    % stimulus
    lat_relative_t = singleStimResults{cell,1}(:,:) - ...
        singleStimResults{cell,1}(1,:);
    lat_relative_e = singleStimResults{cell,2}(:,:) - ...
        singleStimResults{cell,2}(1,:);

    % Collect the median latency differences (2nd to 1st stimulus)
    second_spikes_t(cell) = median(lat_relative_t(2,:),'omitnan');
    second_spikes_e(cell) = median(lat_relative_e(2,:),'omitnan');
    % (3rd to 1st stimulus)
    third_spikes_t(cell) = median(lat_relative_t(3,:),'omitnan');
    third_spikes_e(cell) = median(lat_relative_e(3,:),'omitnan');
    % (4th to 1st stimulus)
    fourth_spikes_t(cell) = median(lat_relative_t(4,:),'omitnan');
    fourth_spikes_e(cell) = median(lat_relative_e(4,:),'omitnan');
end


% Calculate the corrected alpha value (Bonferroni)
alpha_korr = 0.05/6;                        % 0.0083
% ------
                                            % RT            % Cold
signrank(second_spikes_t)                   % 0.0117        % 1.22e-4  *
signrank(third_spikes_t)                    % 0.0137        % 1.22e-4  *
signrank(fourth_spikes_t)                   % 0.0273        % 1.22e-4  *
ranksum(second_spikes_t,third_spikes_t)     % 0.6499        % 0.1543
ranksum(second_spikes_t,fourth_spikes_t)    % 0.6499        % 0.0180
ranksum(third_spikes_t,fourth_spikes_t)     % 0.8798        % 0.2322
% ------
signrank(second_spikes_e)                   % 0.002    *    % 1.22e-4  *
signrank(third_spikes_e)                    % 0.0098        % 1.22e-4  *
signrank(fourth_spikes_e)                   % 0.0059   *    % 1.22e-4  *
ranksum(second_spikes_e,third_spikes_e)     % 0.7913        % 0.7304
ranksum(second_spikes_e,fourth_spikes_e)    % 0.3847        % 0.6960
ranksum(third_spikes_e,fourth_spikes_e)     % 0.6232        % 0.7131


% Combine all spikes in columns 
spikes_t = [second_spikes_t;third_spikes_t;fourth_spikes_t]';
spikes_e = [second_spikes_e;third_spikes_e;fourth_spikes_e]';

% Create figure
fig3 = figure('Position', [218,316,950,431]);
tiledlayout(1,2, 'TileSpacing', 'compact')

% Boxplot for the tactile spike latencies
nexttile; boxplot(spikes_t); ax1 = gca; bxplt1 = ax1.Children;
xlabel("number of separate stimulus"); ax1.XTickLabel = ["2","3","4"];
ylabel(["latency diffenence to 1st", "separate stimulus [ms]"])
title("A", 'FontWeight', 'bold', 'FontSize', 18)
ax1.LineWidth = 2; ax1.FontSize = 14; 
ax1.TitleHorizontalAlignment = 'left'; ax1.TitleFontSizeMultiplier = 1.5;
ax1.Title.Units = 'normalized'; ax1.Title.Position = [-0.24,1,0];


% Boxplot for the electrical spike latencies
nexttile; boxplot(spikes_e); ax2 = gca; bxplt2 = ax2.Children; 
xlabel("number of separate stimulus"); ax2.XTickLabel = ["2","3","4"];
ylabel(["latency diffenence to 1st", "separate stimulus [ms]"])
title("B", 'FontWeight', 'bold')
ax2.LineWidth = 2; ax2.FontSize = 14; 
ax2.TitleHorizontalAlignment = 'left'; ax2.TitleFontSizeMultiplier = 1.5;
ax2.Title.Units = 'normalized'; ax2.Title.Position = [-0.24,1,0];

% Modify boxplot appearence
set(bxplt1.Children(:), 'LineWidth', 1)
set(bxplt1.Children(10:21), 'Color', 'k')
set(bxplt1.Children(1:3), 'Color', [0.5, 0.5, 0.2])
set(bxplt1.Children(7:9), 'Color', [0.4     0.75    0.0])
set(bxplt1.Children(4:6), 'Color', 'k')

set(bxplt2.Children(:), 'LineWidth', 1)
set(bxplt2.Children(10:21), 'Color', 'k')
set(bxplt2.Children(1:3), 'Color', [0.75    0.0     0.35])
set(bxplt2.Children(7:9), 'Color', [0.0     0.6     0.9])
set(bxplt2.Children(4:6), 'Color', 'k')

% Save figure
exportgraphics(fig3, ...
    sprintf("Report/figures/separate_diff_%s.eps",group), ...
    'ContentType','vector','BackgroundColor','white')
