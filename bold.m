%% Function to change figure appearance for usage in texts
% The function changes the line width of all plotted lines, increases the
% marker size, and enlarges the axes labels.
% -----
% Input
%   - ax    axes handle for the axes to adapt
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 24.01.2024
% -------------------------------------------------------------------------
function bold(ax)

% Get children of axis (lines)
chld = ax.Children;

for child = 1:length(chld)
    % Set line width of plot lines
    chld(child).LineWidth = 2;
    chld(child).MarkerSize = 17;
end

% Set line width of axes
ax.LineWidth = 2;
% Adjust font size
ax.FontSize = 14;

% Turn off the box around the plot
ax.Box = "off";
