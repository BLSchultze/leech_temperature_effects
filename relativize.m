function xr = relativize(x)
%% Relativize vector items to the first item in the vector
% -----
% Input
%   x       1xn vector to relativize
% -----
% Output
%   xr      1xn vector with elements of x relative to the last element of 
%           x
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 17.03.2023
% -------------------------------------------------------------------------

% If the vector is not empty ...
if ~isempty(x)
    % Subtract the last vector element from all elements
    xr = x-x(end);
else
    % Otherwise return NaN
    xr = nan;
end