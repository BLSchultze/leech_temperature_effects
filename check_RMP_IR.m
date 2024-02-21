function qualityCheck = check_RMP_IR(files,group,plotflag)
%% Calculate resting potential and input resistance for all trials
% This function calculates the resting membrane potential (RMP) for all
% files requested and for all trails stored within a file. The RMP is
% always calculated in the first and last 500 ms. The input resistance is
% calculated as well. If requested, the results can be visualized,
% seperately for each file. For details on the calculations, see the
% description of the helper function below. 
% -----
% Input
% -----
%   files       1*n vector, the number of the files to be analyzed
%   group       string, stating which experimental group should be 
%               analyzed, possible options: "Cold", "RoomTemp"
%   plotflag    1*1 logical, indicate whether results should be visualized
% ------
% Output
% ------
%   qualityCheck    n*3 cell array, containing the results for all files
%                   and all trials
%                       1. column - resting membrane potential, first 
%                                   500 ms [mV]
%                       2. column - resting membrane potential, last 
%                                   500 ms [mV]
%                       3. column - input resistance [MOhm]
% -------------------------------------------------------------------------
% Author: Bjarne Schultze       last modified: 22.02.2023
% -------------------------------------------------------------------------

% Store the master-path where the data folders are stored
masterPath = ['C:\Users\schul\Nextcloud\Studium\2. Master-Semester\' ...
    'ResearchModule\'];

% Paste search path depending on the dataset
if isequal(group,"Cold")
    dataPath = masterPath+"Cooled\";
elseif isequal(group,"RoomTemp")
    dataPath = masterPath+"RoomTemp\";
% Throw error if non-existing group name was given
else
    % Error with appropriate message
    error("Experiment group must be one of the following: 'Cold', " + ...
        "'RoomTemp'")
end   

% List all files in the folder
fileList = dir(dataPath+"tacData*.mat");
% Load a file with the trials to select
load(dataPath+'tempSelect.mat','tempSelect')

% Pre-define a cell array to collect the results
qualityCheck = cell(length(fileList),3);

% Iterate through all data files
for file = files
    % Load the next file and display name in command window
    load(dataPath+fileList(file).name,'tacData','metadata')
    fprintf("Currently working on file: %s!\n",fileList(file).name)

    % Reduce the recording data according to the 'tempSelect' variable
    tacData = tacData(tempSelect(file,1):tempSelect(file,2),:);

    % Extract the sampling rate
    spRate = metadata.samplingRate;
    
    % Calculate the width of the window to average the membrane potential
    avgWidth = 0.5*spRate;

    % Calculate resting potentials and input resistances for all trials
    [rmpStart,rmpEnd,inputRes] = cellfun(@(x)rmp_ir(x,spRate,avgWidth), ...
        tacData(:,2),'UniformOutput',true);

    % Store the results of the current file
    qualityCheck{file,1} = rmpStart;
    qualityCheck{file,2} = rmpEnd;
    qualityCheck{file,3} = inputRes; 

    % If requested, visualize the results
    if plotflag
        % New figure window with a tiled layout
        figure('Position',[167,77,1161,702])
        tiledlayout(2,1,'TileSpacing','compact')

        % First subplot: resting potential over trial
        nexttile
        plot(rmpStart,'.-'); hold on; plot(rmpEnd,'.-'); hold off
        % Label axes
        xlabel("Trial number"); ylabel("Resting membrane potential [mV]")
        title("File "+fileList(file).name,'Interpreter','none')
        legend("First 500 ms","Last 500 ms")

        % Second subplot: input resistance over trials
        nexttile
        plot(inputRes,'.-')
        % Label the axes
        xlabel("Trial number"); ylabel("Input resistance [M\Omega]")
    end
% End of file iteration loop
end
% End of function definition
end


%% Definition of local helper function

function [rmpStart,rmpEnd,inputRes] = rmp_ir(recData,spRate,avgWidth)
%% Calculate resting membrane potential and input resistance
% The resting membrane potential (RMP) is determiened as the average 
% membrane potential in a given time window. The RMP is determined at the
% beginning at the end of the recoding trace. The input resistance is
% determined by dividing amplitude of the membrane potential during the 
% 0.5 nA stimulation by the stimulation current (R = U/I). For calculating
% the amplitude, the second half of the stimulus (250 ms) is used to avoid 
% taking slow potential changes into account. The membrane potential is 
% averaged. The average membrane potential in the 250 ms before the
% stimulus onset are used as a reference.
% -----
% Input
% -----
%   recData     1*n double, recorded voltage trace [mV]
%   spRate      1*1 double, sampling rate for the given recData 
%               [data points/s]
%   avgWidth    1*1 double, duration for obtaining the resting membrane
%               potential [s]
% ------
% Output
% ------
%   rmpStart    1*1 double, resting membrane potential (average membrane 
%               potential from beginning until avgWidth) [mV]
%   rmpEnd      1*1 double, resting membrane potential (average membrane 
%               potential from end-avgWidth until the end) [mV]
%   inputRes    1*1 double, input resistance [MOhm]
% -------------------------------------------------------------------------

% Calculate resting membrane potential at the beginning and at the end
rmpStart = mean(recData(1:avgWidth));
rmpEnd = mean(recData(end-avgWidth:end));

% Calculate amplitude during the hyperpolarization (second half of the
% stimulus relative to 250 ms before the stimulus onset)
hypAmp = mean(recData(0.75*spRate:spRate)) - ...
    mean(recData(0.25*spRate:0.5*spRate));
% Calculate the input resistance (R = U/I)
inputRes = hypAmp / -0.5;

% End of function definition
end
