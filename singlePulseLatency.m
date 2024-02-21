function [elecLatencies,tacLatencies,nElec,nTac] = singlePulseLatency( ...
    recTrace,onlyElec,onlyTac,searchLen,timeVector)
%% First spike latency for the single stimulus pulses
% This function determines the first spike latencies for the single tactile
% and electrical stimulus pulses (not the stimulus pairs) in the collide
% stimulation protocol (Steffen, v4). This function directly operates on
% the recording trace. The Matlab 'findpeaks' function is used to search
% for spikes in a defined search window after the stimulus onset. The
% function handles one recording trace at a time. 
% -----
% Input
% -----
%   recTrace        1*n vector, recorded membrane potential trace
%   onlyElec        1*k vector, onsets of the k single electrical stimuli
%                   in ms
%   onlyTac         1*k vector, onsets of the k single tactile stimuli,
%                   in ms
%   searchLen       1*1 double, time window after the stimulus onset in
%                   which it should be searched for spikes, in ms
%   timeVector      1*n vector, time vector corresponding to the recording
%                   trace
% ------
% Output
% ------
%   elecLatencies   1*k vector, first spike latencies for the k electrical
%                   stimuli
%   tacLatencies    1*k vector, first spike latencies for the k tactile
%                   stimuli
%   Note: No spikes are encoded and NaN, therefore the dimension of the
%         output vectors always match the dimension of stimulus onset
%         vectors. 
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 18.03.2023
% -------------------------------------------------------------------------

%% Find the spikes
% Set the plotflag
plotflag = 0;

% Calculate the sampling rate
samplingRate = length(timeVector)/30;

% Find all spikes in the recording trace
[peakValues,spikeIndices] = findpeaks(recTrace, ...
    'MinPeakProminence',30,'MinPeakDistance',0.002*samplingRate);

% Initialize a spike counter/index
spikecount = 0;
% Create a variable to collect the spike times
spikeTimes = [];

% Iterate through all spikes
for spike = 1:length(spikeIndices)
    % Only execute if value 1.2 ms before spike is still in the recording
    if spikeIndices(spike) - 0.0012*samplingRate > 1
        % Get the membrane potential 1.2 ms before and after the peak
        l = recTrace(spikeIndices(spike) - 0.0012*samplingRate);
        r = recTrace(spikeIndices(spike) + 0.0012*samplingRate);
        % Calculate the difference to the peak value
        peakDiffl = peakValues(spike) - l;
        peakDiffr = peakValues(spike) - r;

        peakDiff_l(spike) = peakValues(spike) - l;
        peakDiff_r(spike) = peakValues(spike) - r;
        
        % Identify spikes by a membrane potential change of 4 mV on the
        % left and 4 mV on the right side of the peak
        if peakDiffl > 4 && peakDiffr > 4
            % Increase spike counter/index variable
            spikecount = spikecount + 1;
            % Store spike time
            spikeTimes(spikecount) = timeVector(spikeIndices(spike));
        end
    end
end

% If requested, create a plot showing the identified spikes
if plotflag
    figure
    plot(timeVector,recTrace)
    xline(spikeTimes,'Color','g')
    % Label the axes
    xlabel("Time [s]"); ylabel("Membrane potential [mV]")
end

%% Latency calculations
% Pre-define vectors to collect the spike times for the stimuli
elecSpikes = zeros(length(onlyElec),1);
tacSpikes = zeros(length(onlyTac),1);
% Pre-define vectors to collect the spike counts for the stimuli
nElec = zeros(length(onlyTac),1); nTac = zeros(length(onlyTac),1);

% Iterate through all single electrical stimuli
for stim = 1:length(onlyElec)
    % Search for spikes following the current stimulus in a time window
    % defined by 'searchLen'
    spikes = spikeTimes((spikeTimes >= onlyElec(stim)) & ...
        (spikeTimes < onlyElec(stim)+searchLen));
    % If spikes were found, save the time of the first spike
    if isempty(spikes)
        elecSpikes(stim,1) = nan;
    else
        elecSpikes(stim,1) = spikes(1);
    end
    % Determine the number of spikes following an electrical pulse
    nElec(stim,1) = numel(spikes);
end

% Iterate through all single tactile stimuli
for stim = 1:length(onlyTac)
    % Search for spikes following the current stimulus in a time frame
    % defined by 'searchLen'
    spikes = spikeTimes((spikeTimes >= onlyTac(stim)) & ...
        (spikeTimes < onlyTac(stim)+searchLen));
    % If spikes were found, save the time of the first spike
    if isempty(spikes)
        tacSpikes(stim,1) = nan;
    else
        tacSpikes(stim,1) = spikes(1);
    end
    % Determine number of spikes following a tactile pulse
    nTac(stim,1) = numel(spikes);
end

% Calculate the latencies for all stimuli and convert from s to ms
tacLatencies = (tacSpikes - onlyTac)*1000;
elecLatencies = (elecSpikes - onlyElec)*1000;

% End of function definition
end