function plot_stimulation()
%% Plot the stimulation protocol
% Plot the complete elec-tac collide protocol (3x1 tiled layout). Whole 
% protocol is plotted as an overlay of the tactile and electrical 
% stimulation protocols. Additionally, both protocols are shown separately. 
% A detail of the protocol (1 elec and 1 tac pulse) are plotted in a second
% figure next to the measured stimuli.
% -------------------------------------------------------------------------
% Author: Bjarne Schultze   [last modified 26.07.2023]
% -------------------------------------------------------------------------

% Clear workspace
clear

% Load stimulus and create time vector
load("Stimuli\ELEC_collide_v4_2nA_100kHz.mat","stimulus")
stim_elec = stimulus;
load("Stimuli\TAC_collide_v4_100kHz.mat","stimulus")
stim_tac = stimulus*50; 
sampleRate = 100000;
timeVec = (1/sampleRate:1/sampleRate:length(stimulus)/sampleRate)';

% Store the onsets of the separate stimuli
onlyElec = [2,10,18,27];
onlyTac = [3,11,19,28];

% Calculations for the alignment of the y axes
l1 = -1.5; u1 = 2.5;
l2 = 0; u2 = 10;
p = (0 - l2 + ((l2 - u2)* ...
    (l1 - 0))/(l1 - u1))/ ...
    ((l1 - 0)/(l1 - u1) - 1);

% Plot stimulus protocol
fig1 = figure('Position', [218,115,950,632]);
tiledlayout(3,1, 'TileSpacing','compact')

% --- Plot the whole stimulation protocol
nexttile
% Electrcial stimulation
yyaxis left
plot(timeVec, stim_elec, 'LineWidth', 2, 'Color', [135/255,200/255,245/255])
ylim([l1, u1]); ylabel("current [nA]"); yticks([0,2])

% Highlight the separate stimuli
hold on
for i = 1:4
    % Find the indeces
    onlyElec_stim_ind = timeVec>=onlyElec(i)-2*(1/sampleRate) & ...
        timeVec<onlyElec(i)+0.005+2*(1/sampleRate);
    % Extract the data
    onlyElec_stim = stim_elec(onlyElec_stim_ind);
    onlyElec_time = timeVec(onlyElec_stim_ind);

    % Add to plot
    plot(onlyElec_time, onlyElec_stim, 'LineWidth', 2, ...
        'Color', [0/255,80/255,190/255], 'LineStyle', '-', 'Marker', 'none')

end
% Tactile protocol
yyaxis right
plot(timeVec, stim_tac, 'Color', [240/255,161/255,119/255])
ylim([l2-p, u2]); ylabel("force [mN]"), yticks([0,5])
yticklabels(["0","5"])

for i = 1:4
    % Find the indeces
    onlyTac_stim_ind = timeVec>=onlyTac(i)-2*(1/sampleRate) & ...
        timeVec<onlyTac(i)+0.005+2*(1/sampleRate);
    % Extract the data
    onlyTac_stim = stim_tac(onlyTac_stim_ind);
    onlyTac_time = timeVec(onlyTac_stim_ind);

    % Add to plot
    plot(onlyTac_time, onlyTac_stim, 'LineWidth', 2, ...
        'Color', [219/255,91/255,0/255], 'LineStyle', '-', 'Marker', 'none')

end
hold off

% Set the y axes colors
ax1 = gca; ax1.YAxis(1).Color = [0/255,80/255,190/255];
ax1.YAxis(2).Color = [219/255,91/255,0/255];
% Adjust further appearance
bold(gca)

% Add a box for the detail
hold on
plot([3.5, 4.6, 4.6, 3.5, 3.5], [-1.5, -1.5, 9, 9, -1.5], ...
    'LineWidth', 2, 'Color', [48/255,159/255,0/255], 'Marker', 'none')
hold off

% --- Plot the electrical protocol
nexttile
plot(timeVec, stim_elec, 'Color', [135/255,200/255,245/255])
hold on
for i = 1:4
    onlyElec_stim_ind = timeVec>=onlyElec(i)-2*(1/sampleRate) & ...
        timeVec<onlyElec(i)+0.005+2*(1/sampleRate);
    onlyElec_stim = stim_elec(onlyElec_stim_ind);
    onlyElec_time = timeVec(onlyElec_stim_ind);

    plot(onlyElec_time, onlyElec_stim, ...
        'Color', [0/255,80/255,190/255])

end
hold off
ylim([l1, u1]); yticks([0,2]); ylabel("current [nA]")
% Adjust the axes and the color of the y axis
ax2 = gca; ax2.YAxis.Color = [0/255,80/255,190/255];
bold(ax2)

% --- Plot the tactile protocol
nexttile
plot(timeVec, stim_tac, 'Color', [240/255,161/255,119/255])
hold on
for i = 1:4
    onlyTac_stim_ind = timeVec>=onlyTac(i)-2*(1/sampleRate) & ...
        timeVec<onlyTac(i)+0.005+2*(1/sampleRate);

    onlyTac_stim = stim_tac(onlyTac_stim_ind);
    onlyTac_time = timeVec(onlyTac_stim_ind);

    plot(onlyTac_time, onlyTac_stim, ...
        'Color', [219/255,91/255,0/255])
end
hold off
xlabel("time [s]")
ylim([l2-p, u2]); ylabel("force [mN]"), yticks([0,5])
yticklabels(["0","5"])
% Change color of y axis and adjust axes
ax3 = gca; ax3.YAxis.Color = [219/255,91/255,0/255];
bold(ax3)

% Save figure
exportgraphics(fig1, "Report/figures/stimulation.png", ...
    'ContentType','vector','BackgroundColor','white')
%  print("Report/figures/stimulation.png", '-dpng')

%% Detail of the stimulation protocol

% Load an example data file and extract the necessary data
load("Cooled_raw\22-11-22-13-14-31\tac_cold-22-11-22-13-14-31.mat", ...
    "tac_experiment")
timeVecTrace = tac_experiment(2).result.timeVector;
stim_tac = tac_experiment(2).result.stimulus(:,4) * 50;
stim_tac_rec = tac_experiment(2).result.recording(:,6);
stim_elec = tac_experiment(2).result.stimulus(:,3);
stim_elec_rec = tac_experiment(2).result.recording(:,2);

% Set up figure
fig2 = figure('Position', [218,115,950,632]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact')
% Calculations for y axes alignment
l1 = -0.7; u1 = 2.8;
l2 = -3; u2 = 10;
p = (0 - l2 + ((l2 - u2)* ...
    (l1 - 0))/(l1 - u1))/ ...
    ((l1 - 0)/(l1 - u1) - 1);

% --- requested stimulus
nexttile()
% 1st y-axis: electrical stimulus over time
yyaxis('left')
plot(timeVecTrace(timeVecTrace > 3.965 & timeVecTrace < 3.985), ...
    stim_elec(timeVecTrace > 3.965 & timeVecTrace < 3.985,1), ...
    'LineWidth', 2, 'Color', [135/255,200/255,245/255])
% Adjust the axes limits and label the y-axis
ylim([l1, u1]); yticks([0,2])
xlim([3.96,4.04])
ylabel("current [nA]")

% 2nd y-axis: tactile stimulus over time
yyaxis('right')
plot(timeVecTrace(timeVecTrace > 3.995 & timeVecTrace < 4.035), ...
    stim_tac(timeVecTrace > 3.995 & timeVecTrace < 4.035), ...
    'LineWidth', 2, 'Color', [240/255,161/255,119/255])
% Adjust the y-axis limits and label the axis
ylim([l2-p, u2]); yticks([0,4,8])
ylabel("force [mN]")

% Add a line for delta t
hold on
plot([3.97, 4], [8, 8], 'Marker', '|', 'Color', 'k', 'LineStyle','-')
text(mean([3.97, 4]), 8.7, "$\Delta$t = -30 ms", 'FontSize', 14, ...
    'HorizontalAlignment','center', 'Interpreter','latex')
hold off

% Set the y axes colors
ax1 = gca; ax1.YAxis(1).Color = [0/255,80/255,190/255];
ax1.YAxis(2).Color = [219/255,91/255,0/255];
% Adjust axes appearance
ax1 = gca; bold(ax1)

% --- actual stimulus
nexttile
% 1st y-axis: tactile stimulus over time
yyaxis('left')
plot(timeVecTrace(timeVecTrace > 3.965 & timeVecTrace < 3.985), ...
    stim_elec_rec(timeVecTrace > 3.965 & timeVecTrace < 3.985), ...
    'LineWidth', 2, 'Color', [135/255,200/255,245/255])
% Adjust the axes limits and label the y-axis
ylim([l1, u1]); yticks([0,2])
xlim([3.96,4.04])
ylabel("current [nA]")

% 2nd y-axis: electrical stimulus over time
yyaxis('right')
plot(timeVecTrace(timeVecTrace > 3.995 & timeVecTrace < 4.035), ...
    stim_tac_rec(timeVecTrace > 3.995 & timeVecTrace < 4.035), ...
    'LineWidth', 2, 'Color', [240/255,161/255,119/255])
% Adjust the y-axis limits and label the axis
ylim([l2-p, u2]); yticks([0,4,8])
ylabel("force [mN]")
xlabel("time [s]")

% Set the y axes colors
ax2 = gca; ax2.YAxis(1).Color = [0/255,80/255,190/255];
ax2.YAxis(2).Color = [219/255,91/255,0/255];
% Adjust axes appearance
bold(ax2)

% Save figure 
exportgraphics(fig2, "Report/figures/stimulation_detail.png", ...
    'ContentType','vector', 'BackgroundColor', 'white')
% print("Report/figures/stimulation_detail.png", '-dpng')

% Print success message to console
disp("Successfully saved figures: 'stimulation.eps' and " + ...
    "'stimulaiton_detail.eps to Report/figures'!")
