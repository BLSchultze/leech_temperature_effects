function varargout = plotResultsSingleStim(group)
%% Visualize the results for the single pulses
% This funciton creates various graphical representations of the results of
% the single pulse analysis. The script takes either the results of the
% cold or the room temperature experiments. The plots that are created are
% listed below.  
% -----------------
% Table of contents
% -----------------
% 1.a) Single stimuli latencies as a function of temperature
% 1.b) Single stimuli latencies as a function of temperature 
%       (colorcoded for temperature raise/fall)
% 1.c) Single stimuli latencies as a function of temperature
%       (colorcoded for the amplitude of the stimulation current)
% 2.   Single stimuli latencies as a function of the trial
% 3.a) Single pulse latencies relative to warmest trial (separately)
% 3.b) Single pulse latencies relative to warmest trial (one plot)
% 4.   Single pulse latencies as a function of the pulse index (all trials)
% 5.a) Histograms of the Q10 values (based on the single stimuli)
% 5.b) Boxplot for the Q10 values
% 6.   Number of spikes as a function of temperature
% 7.   Number of spikes as a function of the trial index
% 8.   Number of spikes for each pulse over pulse index
% 
% -----
% Input
%   group       indicating the experimental group, either 'Cold' or
%               'RoomTemp'
% ------
% Output
%   varargout   optional, figure handles for the figures from parts 1-5,
%               number in the list above gives the output index
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 28.02.2023
% -------------------------------------------------------------------------

%% Preparations

% Change working directory
cd(['C:\Users\schul\Nextcloud\Studium\2. Master-Semester\' ...
    'ResearchModule_Leech\'])

if isequal(group,"Cold")
    % Load the file with the analysis results (cold data)
    load("analysisResultsCold.mat",'singleStimResults', ...
        'relLat_elec','relLat_tac')
elseif isequal(group,"RoomTemp")
    % Load the file with the analysis results (room temperature data)
    load("analysisResultsRoomTemp.mat",'singleStimResults', ...
        'relLat_elec','relLat_tac')
end

% Define a set of colors for a consistent visualization
masterColors = [0.0     0.6     0.9; ...
                0.4     0.75    0.0; ...
                0.9     0.6     0.0];

% Define position values for new figure windows
masterPosition = [167,77,1161,702];

%% 1.a) Single stimuli latencies as a function of temperature 
% Create a new figure window and set hold on
fig1 = figure('Position',masterPosition,'Name',"Tac latencies ~ temp");
hold on

% Iterate through the results of all data files and plot the mean
% latencies as a function of the temperature
for res = 1:size(singleStimResults,1)
    plot(singleStimResults{res,7},singleStimResults{res,1}(2,:), ...
        '.-', 'MarkerSize', 14)
end
% Label the plot
title("Tactile stimulation")
xlabel("Temperature [°C]")
ylabel("Latency [ms]")


% Create a new figure window and set hold on
fig2 = figure('Position',masterPosition,'Name',"Elec latencies ~ temp");
hold on
% Label the plot
for res = 1:length(singleStimResults)
    plot(singleStimResults{res,7},singleStimResults{res,2}(2,:), ...
        '.-', 'MarkerSize', 14)
end
% Label the plot
title("Electrical stimulation")
xlabel("Temperature [°C]")
ylabel("Latency [ms]")

% Set hold off
hold off

% Catch figure handles
varargout{1} = fig1; varargout{1}(2) = fig2;

%% 1.b) Single stimuli latencies as a function of temperature
% --- Curves colorcoded according to temperature raise or fall during the 
%     experiment ---

if isequal(group,"Cold")
    % Create a new figure window and set hold on
    fig3 = figure('Position',masterPosition,'Name',"Tac latencies ~ temp");
    hold on
    
    % Iterate through the results of all data files and plot the mean
    % latencies as a function of the temperature
    for res = 1:size(singleStimResults,1)
        % Distinguish (color) between temperature raise and fall during an 
        % experiment
        if singleStimResults{res,10}(1) < singleStimResults{res,10}(end)
            plt_raise = plot(singleStimResults{res,7}, ...
                singleStimResults{res,1}(2,:),'.-', ...
                'Color',[0.74,0.25,0.57], 'MarkerSize', 14);
        else
            plt_fall = plot(singleStimResults{res,7}, ...
                singleStimResults{res,1}(2,:),'.-', ...
                'Color',[0.28,0.64,0.71], 'MarkerSize', 14);
        end
    end
    % Label the plot
    title("Tactile stimulation")
    xlabel("Temperature [°C]")
    ylabel("Latency [ms]")
    legend([plt_raise,plt_fall],"Ascending temperature", ...
        "Descending temperature",'Location','northeast')
    
    
    % Create a new figure window and set hold on
    fig4 = figure('Position',masterPosition,'Name',"Elec latencies ~ temp");
    hold on
    % Label the plot
    for res = 1:length(singleStimResults)
        % Distinguish (color) between temperature raise and fall
        if singleStimResults{res,10}(1) < singleStimResults{res,10}(end)
            plt_raise = plot(singleStimResults{res,7}, ...
                singleStimResults{res,2}(2,:),'.-', ...
                'Color',[0.74,0.25,0.57], 'MarkerSize', 14);
        else
            plt_fall = plot(singleStimResults{res,7}, ...
                singleStimResults{res,2}(2,:),'.-', ...
                'Color',[0.28,0.64,0.71], 'MarkerSize', 14);
        end
    end
    % Label the plot
    title("Electrical stimulation")
    xlabel("Temperature [°C]")
    ylabel("Latency [ms]")
    legend([plt_raise,plt_fall],"Ascending temperature", ...
        "Descending temperature",'Location','northeast')
    
    % Set hold off
    hold off
    
    % Catch axes handles
    varargout{1}(3) = fig3; varargout{1}(4) = fig4;
end

%% 1.c) Single stimuli latencies as a function of temperature
% --- colorcoded for the amplitude of the stimulation current ---

% Give the information, in which experiments which stimulus was used
if isequal(group,"Cold")
    % File numbers for the cold group
    stim1nA = [1,3,4,5,6,7,9,11,15,16];
    % Set difference gives the group with the other stimulus
    stim2nA = setdiff(1:length(relLat_tac),stim1nA);
elseif isequal(group,"RoomTemp")
    % File numbers for the control group
    stim1nA = [1,2,3,4,5,6,7,8,9,10,11];
    % Use set difference to get the file numbers for the other stimulus
    stim2nA = setdiff(1:length(relLat_tac),stim1nA);
end

% Create a new figure and set hold on
fig5 = figure('Position',masterPosition,'Name',"Tac latencies ~ temp");
hold on

% Iterate through the results of all data files and plot the mean
% 'tactile' latencies as a function of the temperature
for res = 1:size(singleStimResults,1)
    % For the 1 nA group
    if ismember(res,stim1nA)
        p1= plot(singleStimResults{res,7},singleStimResults{res,1}(2,:), ...
            '.-','Color',[0.91,0.73,0.00], 'MarkerSize', 14);
    % For the 2 nA group (1.5 not existing in this dataset)
    elseif ismember(res, stim2nA)
        p2 = plot(singleStimResults{res,7},singleStimResults{res,1}(2,:), ...
            '.-','Color',[0.55,0.23,0.95], 'MarkerSize', 14);
    end
    
end
% Label the plot
title("Tactile stimulation")
xlabel("Temperature [°C]")
ylabel("Latency [ms]")
if exist("p1","var") && exist("p2","var")
    legend([p1,p2],["1 nA stimuli","2 nA stimuli"],'Location','northeast')
end

% Create a new figure window and set hold on
fig6 = figure('Position',masterPosition,'Name',"Elec latencies ~ temp");
hold on

% Iterate through the results of all data files and plot the mean
% 'electrical' latencies as a function of the temperature
for res = 1:length(singleStimResults)
    % For the 1 nA group
    if ismember(res,stim1nA)
        p1 = plot(singleStimResults{res,7},singleStimResults{res,2}(2,:), ...
            '.-','Color',[0.91,0.73,0.00], 'MarkerSize', 14);
    % For the 2 nA group (1.5 nA not in the dataset)
    else
        p2 = plot(singleStimResults{res,7},singleStimResults{res,2}(2,:), ...
            '.-','Color',[0.55,0.23,0.95], 'MarkerSize', 14);
    end
end

% Set hold off
hold off

% Label the plot
title("Electrical stimulation")
xlabel("Temperature [°C]")
ylabel("Latency [ms]")
if exist("p1","var") && exist("p2","var")
    legend([p1,p2],["1 nA stimuli","2 nA stimuli"],'Location','northeast')
end

% Catch figure handle
varargout{1}(5) = fig5; varargout{1}(6) = fig6;

%% 2. Single stimuli latencies as a function of the trial

% Create a new figure window and set hold on
figure('Position',masterPosition,'Name',"Tac latencies ~ trial")
hold on

% Iterate through the results of all data files and plot the mean
% latencies as a function of the temperature
for res = 1:size(singleStimResults,1)
    % Re-order the dataset (back to sorting by trials)
    [~,reOrder] = sort(singleStimResults{res,10});
    sStim_reOrder = singleStimResults{res,1}(2,reOrder);
    % Plot as a funtion of the index
    plot(sStim_reOrder,".-", 'MarkerSize', 14)
end
% Label the plot
title("Tactile stimulation")
xlabel("Trial index (recording order)")
ylabel("Latency [ms]")

% Catch the figure handle
varargout{2}(1) = gcf;

% Create a new figure window and set hold on
figure('Position',masterPosition,'Name',"Elec latencies ~ trial")
hold on
% Label the plot
for res = 1:length(singleStimResults)
    % Re-order the dataset (back to sorting by trials)
    [~,reOrder] = sort(singleStimResults{res,10});
    sStim_reOrder = singleStimResults{res,2}(2,reOrder);
    % Plot as a function of the index
    plot(sStim_reOrder,".-", 'MarkerSize', 14)
end
% Label the plot
title("Electrical stimulation")
xlabel("Trial index (recording order)")
ylabel("Latency [ms]")

% Set hold off
hold off

% Catch the figure handle
varargout{2}(2) = gcf;


%% 3.a) Single pulse latencies relative to warmest trial (separately)
% Set up a new figure window 
figure('Position',masterPosition,'Name',"Relative latencies ~ temp")
tiledlayout(1,3,'TileSpacing','compact')

% --- Tactile stimulaiton
% Re-shape the data into to long column vectors
xTac = cell2mat(relLat_tac(:,2)')'; yTac = cell2mat(relLat_tac(:,1)')';
% Exclude missing values
xTac = xTac(~isnan(yTac)); yTac = yTac(~isnan(yTac));
% Calculate a linear fit and correlation coefficient
[fitTac,~] = fit(xTac,yTac,'poly1');
[correlationTac, pval_corr_tac] = corr(xTac,yTac,'type','Spearman');

% Plot the latency over temperature relative to warmest trial, for all
% datasets (tactile stimulation)
nexttile; hold on;
for dataset=1:length(relLat_tac)
    plot(relLat_tac{dataset,2},relLat_tac{dataset,1},'o', ...
        'Color',masterColors(2,:))
end

% Add the line for the linear fit
plt_fit_tac = plot(fitTac);
plt_fit_tac.Color = [0.3,0.6,0]; plt_fit_tac.LineWidth = 2;

% Insert a reference line at 0 ms
yline(0,'Color','k','LineWidth',2)
% Label the plot
title("Tactile stimulation")
xlabel("Temperature diffenence to warmest trial [°C]")
ylabel("Latency difference to warmest trial [ms]")
legend(plt_fit_tac,'Linear fit')
% Adjust the limits for the y-axis
ylim([-6,14])

% --- Electrical stimulation
% Re-shape the data into to long column vectors
xElec = cell2mat(relLat_elec(:,2)')'; yElec = cell2mat(relLat_elec(:,1)')';
% Exclude missing values
xElec = xElec(~isnan(yElec)); yElec = yElec(~isnan(yElec));
% Calculate a linear fit and correlation coefficient
[fitElec,~] = fit(xElec,yElec,'poly1');
[correlationElec,pval_corr_elec] = corr(xElec,yElec, 'type', 'Spearman');

% Plot the latency over temperature relative to warmest trial, for all
% datasets (electrical stimulaiton)
nexttile; hold on
for dataset=1:length(relLat_elec)
    plot(relLat_elec{dataset,2},relLat_elec{dataset,1},'o', ...
        'Color',masterColors(1,:))
end

% Add the line for the linear fit
plt_fit_elec = plot(fitElec);
plt_fit_elec.Color = [0,0.45,0.65]; plt_fit_elec.LineWidth = 2;

% Insert a reference line at 0 ms
yline(0,'Color','k','LineWidth',2)
% Label the plot
title("Electrical stimulation")
xlabel("Temperature difference to warmest trial [°C]")
ylabel("Latency difference to warmest trial [ms]")
legend(plt_fit_elec,"Linear fit")
% Adjust the limits for the y-axis
ylim([-6,14])

% --- Both stimulations combined
% Plot the latency over temperature relative to warmest trial, for all
% datasets (electrical stimulaiton)
nexttile; hold on
for dataset=1:length(relLat_tac)
    plot(relLat_elec{dataset,2},relLat_elec{dataset,1},'o', ...
        'Color',masterColors(1,:))
    plot(relLat_tac{dataset,2},relLat_tac{dataset,1},'o', ...
        'Color',masterColors(2,:))
end

% Add the line for the linear fit (tactile stimulation)
plt_fit_tac = plot(fitTac);
plt_fit_tac.Color = [0.3,0.6,0]; plt_fit_tac.LineWidth = 2;
% Add the line for the linear fit (electrical stimulation)
plt_fit_elec = plot(fitElec);
plt_fit_elec.Color = [0,0.45,0.65]; plt_fit_elec.LineWidth = 2;

% Insert a reference line at 0 ms
yline(0,'Color','k','LineWidth',2)
% Label the plot
title("Single tactile and electrical stimuli")
xlabel("Temperature difference to warmest trial [°C]")
ylabel("Latency difference to warmest trial [ms]")
legend([plt_fit_elec,plt_fit_tac],"Fit electrical stimulation", ...
    "Fit tactile stimulaiton")
% Adjust the limits for the y-axis
ylim([-6,14])

% Add the R^2 values
text(-7,-4, "R_s = "+round(correlationElec,2),'Color',[0,0.45,0.65])
text(-7,-5, "R_s = "+round(correlationTac,2),'Color',[0.3,0.6,0])

fprintf("\nP-value for temperature-latency correlation: \n" + ...
    "electrical stimulation: %f \n" + ...
    "for tactile stimulation: %f \n", pval_corr_elec,pval_corr_tac)


%% 3.b) Single pulse latencies relative to warmest trial (one plot)
% Set up a new figure window and put hold on
figure('Position',[321,75,841,702],'Name',"Relative latencies ~ temp")
hold on

% Plot the latency over temperature relative to warmest trial, for all
% datasets (electrical stimulaiton)
for dataset=1:length(relLat_tac)
    % Change plot symbol according to ascending or descending temperature
    % during the recording
    if singleStimResults{dataset,10}(1) < ...
            singleStimResults{dataset,10}(end)
        plot(relLat_elec{dataset,2},relLat_elec{dataset,1},'x', ...
            'Color',masterColors(1,:))
        plot(relLat_tac{dataset,2},relLat_tac{dataset,1},'x', ...
            'Color',masterColors(2,:))
    else
        plot(relLat_elec{dataset,2},relLat_elec{dataset,1},'.', ...
            'Color',masterColors(1,:),'MarkerSize',10)
        plot(relLat_tac{dataset,2},relLat_tac{dataset,1},'.', ...
            'Color',masterColors(2,:),'MarkerSize',10)
    end
    
end

% Add the line for the linear fit (tactile stimulation)
plt_fit_tac = plot(fitTac);
plt_fit_tac.Color = [0.3,0.6,0]; plt_fit_tac.LineWidth = 2;
% Add the line for the linear fit (electrical stimulation)
plt_fit_elec = plot(fitElec);
plt_fit_elec.Color = [0,0.45,0.65]; plt_fit_elec.LineWidth = 2;

% Add two black symbols outside the visible area to have black symbols in
% the legend (I don't know a more elegant way)
pltDesc = plot(1,25,'.k','MarkerSize',8); 
pltAsc = plot(1,25,'xk','MarkerSize',8);

% Label the plot
title("Single tactile and electrical stimuli")
xlabel("Temperature difference to warmest trial [°C]")
ylabel("Latency difference to warmest trial [ms]")
legend([plt_fit_elec,plt_fit_tac,pltDesc,pltAsc], ...
    ["Electrical"+newline+"stimulation","Tactile"+newline+"stimulation", ...
    "Descending"+newline+"temperature","Ascending"+newline+"temperature"])

% Add the R^2 values
text(-6.8,-0.8, "r_s = "+round(correlationElec,2),'Color',[0,0.45,0.65], ...
    'FontSize', 14)
text(-6.8,-1.5, "r_s = "+round(correlationTac,2),'Color',[0.3,0.6,0], ...
    'FontSize', 14)

% Adjust the limits for the y-axis
ylim([-2,14])
% Catch figure handle
varargout{3} = gcf;

%% 4. Single pulse latencies as a function of the pulse index (all trials)

% Set up a new figure window and set hold on
figure('Position',masterPosition,'Name',"Elec latencies ~ index")
hold on

% Iterate through all datasets
for dataset = 1:size(singleStimResults,1)
    % Re-order the dataset (back to sorting by trials)
    [~,reOrder] = sort(singleStimResults{dataset,10});
    sStim_reOrder = singleStimResults{dataset,2}(:,reOrder);

    % Re-shape the matrix with all single pulse latencies to a vector
    allElecLat = reshape(sStim_reOrder,1,[]);
    % Plot the latencies over their vector indices
    plot(allElecLat,'.-') 
end

% Set hold off 
hold off
% Catch the axes handle
ax = gca;
% Insert a vertical line at each position of a new trial
xl = xline(5:4:ax.XLim(2));

% Insert a legend for the vertical lines
legend(xl(1),"1st pulse of each trial")


% Set up a new figure window and set hold on
figure('Position',masterPosition,'Name',"Tac latencies ~ index")
hold on

% Iterate through all datasets and plot all single pulse latencies
for dataset = 1:size(singleStimResults,1)
    % Re-order the dataset (back to sorting by trials)
    [~,reOrder] = sort(singleStimResults{dataset,10});
    sStim_reOrder = singleStimResults{dataset,1}(:,reOrder);

    % Re-shape the matrix with all single pulse latencies to a vector
    allTacLat = reshape(sStim_reOrder,1,[]);
    plot(allTacLat,'.-') 
end
% Set hold off, catch the axes handle and insert vertical lines
hold off
ax2 = gca;
xl2 = xline(5:4:ax2.XLim(2));

% Insert a legend for the vertical lines
legend(xl2(1),"1st pulse of each trial")

% Label the axes in both figure windows
xlabel([ax,ax2],"Pulse index (over all trials, recording order)")
ylabel([ax,ax2],"Latency [ms]")
title(ax,"Electrical stimulation only")
title(ax2,"Tactile stimulation only")
% Store figure handle
varargout{4} = gcf;

%% 5.a) Histograms of the 10 °C change rates (based on the single stimuli)
if isequal(group,"Cold")
    % Set uo a new figure window
    figure('Position',masterPosition,'Name',"Q_10 values"); 
    tiledlayout(2,1,'TileSpacing','compact')
    
    % Define limits for x axis
    xlimits = [round(min(cell2mat(singleStimResults(:,9))),1)-0.1, ...
        round(max(cell2mat(singleStimResults(:,9))),1)+0.1];
    
    % 1st subplot: histogram of the values for electrical stimulation
    nexttile
    histogram(cell2mat(singleStimResults(:,9)),'BinWidth',0.1, ...
        'Normalization','probability')
    xline(1,'Color',[0.7,0.7,0.7])
    % Insert a line for the median
    medLineE = xline(median(cell2mat(singleStimResults(:,9))),'Color', ...
        masterColors(3,:),'LineWidth',2);
    xlim(xlimits)
    % Label the subplot
    title("Electrical stimulation")
    xlabel("10 °C change rate")
    ylabel("probability")
    legend(medLineE,"Median")
    
    % 2nd subplot: histogram of the values for tactile stimulation
    nexttile
    histogram(cell2mat(singleStimResults(:,8)),'BinWidth',0.1, ...
        'Normalization','probability')
    xline(1,'Color',[0.7,0.7,0.7])
    % Insert a line for the median
    medLineT = xline(median(cell2mat(singleStimResults(:,8))),'Color', ...
        masterColors(3,:),'LineWidth',2);
    xlim(xlimits)
    % Label the subplot
    title("Tactile stimulation")
    xlabel("10 °C change rate")
    ylabel("probability")
    legend(medLineT,"Median")
end

%% 5.b) Boxplot for the Q10 values
if isequal(group,"Cold")
    % New figure window
    figure 
    
    % Create boxplot for the Q10 values
    boxplot([cell2mat(singleStimResults(:,9)), ...
        cell2mat(singleStimResults(:,8))],'Labels', ...
        ["Electrical stimulation","Tactile stimulation"], 'Positions', [1,1.5])
    ylabel("Q_{10} value")
    
    % Store the figure handle
    varargout{5} = gcf;
end

%% 6. Number of spikes as a function of temperature 
COrder = [0	0.447	0.741;
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

% Set up a new figure window
figure('Position',masterPosition,'Name','SpikeCountTemp')
tiledlayout(2,1,'TileSpacing','compact')
tileE = nexttile; tileT = nexttile;
hold([tileE,tileT],'on')

% Iterate through all dataset
for dataset = 1:size(singleStimResults,1)
    % Plot the data for electrical stimulation
    plot(tileE,singleStimResults{dataset,7}, ...
        singleStimResults{dataset,6}(2,:),'.-', 'LineWidth',2, ...
        'MarkerSize',17);
    % Plot the data for tactile stimulation
    plot(tileT,singleStimResults{dataset,7}, ...
        singleStimResults{dataset,5}(2,:),'.-', 'LineWidth',2, ...
        'MarkerSize',17);
end

% Label the subplots
xlabel(tileT,"Temperature [°C]")
ylabel([tileE,tileT],"Number of spikes per pulse")
title(tileE,"A")
title(tileT,"B")


%% 7. Number of spikes as a function of the trial index
% Set up a new figure window
figure('Position',masterPosition,'Name','SpikeCountTrial')
tiledlayout(2,1,'TileSpacing','compact')
tileE = nexttile; tileT = nexttile;
hold([tileE,tileT],'on')

for dataset = 1:size(singleStimResults,1)
    % Re-order the dataset (back to sorting by trials)
    [~,reOrder] = sort(singleStimResults{dataset,10});
    nSpkE_reOrder = singleStimResults{dataset,6}(:,reOrder);
    nSpkT_reOrder = singleStimResults{dataset,5}(:,reOrder);

    % Plot as a funtion of the index    
    plot(tileE,nSpkE_reOrder(2,:),'.-','LineWidth',2, ...
        'MarkerSize',17);
    plot(tileT,nSpkT_reOrder(2,:),'.-','LineWidth',2, ...
        'MarkerSize',17);
end

% Label the subplots
xlabel(tileT,"Trial index (recording order)")
ylabel([tileE,tileT],"Number of spikes per pulse")
title(tileE,"A")
title(tileT,"B")

%% 8. Number of spikes for each pulse over pulse index
% Set up a new figure window
figure('Position',masterPosition)
tiledlayout(2,1,'TileSpacing','compact')
tileE = nexttile; tileT = nexttile;
hold([tileE,tileT],'on')

for dataset = 1:size(singleStimResults,1)
    % Re-order the dataset (back to sorting by trials)
    [~,reOrder] = sort(singleStimResults{dataset,10});
    nSpkE_reOrder = singleStimResults{dataset,6}(:,reOrder);
    nSpkT_reOrder = singleStimResults{dataset,5}(:,reOrder);

    % Re-shape the matrix with all single pulse latencies to a vector
    nSpkE = reshape(nSpkE_reOrder,1,[]);
    nSpkT = reshape(nSpkT_reOrder,1,[]);

    % Plot as a funtion of the index    
    plot(tileE,nSpkE,'.-');
    plot(tileT,nSpkT,'.-');
end

% Label the subplots
xlabel([tileE,tileT],"Pulse index (recording order)")
ylabel([tileE,tileT],"Number of spikes per pulse")
title(tileE,"Electrical stimulation")
title(tileT,"Tactile stimulation")

%% Spike count change over temperature change 
figure
hold on

% Iterate through all datasets
for dataset = 1:size(singleStimResults,1)
    % Calculate differences in temperature and spike count to the warmest
    % trial
    spkDiff_t = relativize(median(singleStimResults{dataset,5}));
    spkDiff_e = relativize(median(singleStimResults{dataset,6}));
    tempDiff = relativize(singleStimResults{dataset,7});
    % Plot the spike count difference over the temperature difference
    plot(tempDiff,spkDiff_e,'.','MarkerSize',8, ...
        'Color',masterColors(1,:));
    plot(tempDiff,spkDiff_t,'.','MarkerSize',8, ...
        'Color', masterColors(2,:));
end

% Label the axes
xlabel("Temperature difference [°C]")
ylabel("Difference in median spike count")
