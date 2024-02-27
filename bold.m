%% Function to change figure appearance for usage in texts
% The function changes the line width of all plotted lines, increases the
% marker size, and enlarges the axes labels.
% -----
% Input
%   - ax            axes handle for the axes to adapt
%   - markerSize    (optional, positional) float, gives the marker size
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 27.02.2024
% -------------------------------------------------------------------------
function bold(ax,varargin)

% Check for optional input
if nargin > 1
    markerSize = varargin{1};
else
    markerSize = 17;
end

% Get children of axis (lines)
chld = ax.Children;

% Iterate over axis children, distinguish line and text
for child = 1:length(chld)
    if isequal(chld(child).Type,'line')
        % Set line width of plot lines
        chld(child).LineWidth = 2;
        chld(child).MarkerSize = markerSize;
    elseif isequal(chld(child).Type,'text')
        % Set font size of text
        chld(child).FontSize = 14;
    end
end

% Set line width of axes
ax.LineWidth = 2;
% Adjust font size
ax.FontSize = 14;

% Turn off the box around the plot
ax.Box = "off";
