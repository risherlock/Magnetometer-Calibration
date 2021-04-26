% Add all subfolders to path

% Plot setting
set(groot,'defaulttextinterprete','latex');  
set(groot, 'defaultAxesTickLabelInterprete','latex');  
set(groot, 'defaultLegendInterprete','latex');

path = genpath('D:\GitHub\Compass-Calibration\3-axis');
addpath(path);