function colors = colorGradient(cStart,cEnd,nColors,varargin)
%% Color gradient
% This function creates color codes that range between to defined colors in
% any number of given steps. The colors are linearly spaced. The result can
% be plotted if requested. 
% -----
% Input
% -----
%   cStart      1x3 vector, contains the RBG code of the color that the
%               gradient should start with
%   cEnd        1x3 vector, contains the RGB code of the color that the 
%               gradient should end with
%   nColors     1x1 doudle
%   precision   optional, positional argument, gives the number of decimal 
%               digits to round the color specification to
%   'Plotflag'  optional, name-value argument, states whether or not to
%               create a plot showing the result, default: no plot, 1 for 
%               plot
% ------
% Output
% ------
%   colors      nColors*3 matrix, contains the colors as RGB triplets
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 14.03.2023
% -------------------------------------------------------------------------

%% Check for optional input arguments (basic query) 
% If mor than the 3 standard arguments are given and the first optional
% argument is numeric ...
if nargin > 3 && isnumeric(varargin{1})
    % the first argument is interpreted as the precision value
    precision = varargin{1};
    % If there is another optional value and it reads 'PlotFlag' ...
    if nargin > 4 && isequal(varargin{2},"Plotflag")
        % The plotflag logical is set to true/1 (plot)
        plotflag = 1;
    else
        % plotflag false/0 (no plot) if no argument given 
        plotflag = 0;
    end
% If no optional arguments were given ...
else
    % The precision is set to the standard of 3
    precision = 3;
    % plotflag is set to false/0 (no plot)
    plotflag = 0;
end
    
%% Create colors
% Create linearly spaced sequences from start to end value
redSteps = linspace(cStart(1),cEnd(1),nColors)';    % red component
greenSteps = linspace(cStart(2),cEnd(2),nColors)';  % green component
blueSteps = linspace(cStart(3),cEnd(3),nColors)';   % blue component

% Combine the sequences for red, green, and blue component to a RGB matrix
colors = [redSteps,greenSteps,blueSteps];
% Round the color codes to the requested number of decimal digits
% (precision)
colors = round(colors,precision);

% If requested, create a figure that shows the created colors
if plotflag
    % Set up a figure window
    figure('Name','Color Gradient Function','NumberTitle','off', ...
        'Color', 'w','MenuBar','none'); 

    % Use the left y axis for plotting
    yyaxis("left"); hold on 
    % Plot a line with squares for each color that was created
    for i = 1:nColors
        plot(0:0.1:1,repmat(i,11,1),'s-','Color',colors(i,:),'LineWidth',3)
    end
    % Label the y axis
    ylabel("Color index")
    % Adapt the axes limits
    xlim([-0.1,1.1]); ylim([0.8,nColors+0.2])
    % Catch the axes handle 
    ax = gca; 
    % Change the plot position
    ax.Position = [0.1,0.05,0.65,0.9];
    % Adjust the ticks on both axes
    ax.XTick = []; ax.YTick = 1:nColors; ax.YColor = 'k';

    % Activate the right y axis
    yyaxis("right")
    % Label the axis
    ylabel("RBG color code")
    % Create text for the ticks containing the RGB color code
    yText = round(colors,2);
    yText = mat2cell(yText,ones(1,nColors),3);
    yText = cellfun(@mat2str,yText,'UniformOutput',false);

    % Adjust the axis limits and the color
    ylim([0.8,nColors+0.2]); ax.YColor = 'k';
    % Change the tick values and the labels
    ax.YTick = 1:nColors;
    ax.YTickLabel = yText;
end

% End of function definition 
end