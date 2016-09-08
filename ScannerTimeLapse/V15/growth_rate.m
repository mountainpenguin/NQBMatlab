function growth_rate(varargin)
%% growth_rate(DirNames, ['method', method], ['debug', debug])
% -----------------------------------------------------------------------
% Purpose: Plot various histograms including colony growth rates
%
% Description: Determines appearance times, growth rates, and growth times
%               for each of the inputted data sets, and plots them. 
%
% Arguments: DirNames (optional)
%               If not provided, defaults to current working directory
%               If a single string, only processes that path
%               If a cell array (curly brackets), processes each directory and
%               combines them into a single data set
%            method (optional)
%               Defaults to 1: time for sixfold increase from first appearance
%               If set to 2: averaged time for sixfold increase from all
%               timepoints
%            debug (optional)
%               Defaults to 0
%               If set to 1, will create a figure to demonstrate how growth
%               rates are calculated for each colony (can create a lot...)
%
% Examples:
%   Single dataset:
%       growth_rate('/full/path/to/data');
%
%   Multiple datasets:
%       growth_rate({'/full/path/to/data1', '/full/path/to/data2'});
%
%   Single dataset, with debugging:
%       growth_rate('/full/path/to/data', 'debug', 1);
%
%   Multiple datasets, with debugging:
%       growth_rate({'/full/path/to/data1', '/full/path/to/data2'}, 'debug', 1);
%
% Notes:
%   Ignores any colonies that are not present in the last frame.
%   Ignores any colonies that aren't in existance for less that
%   eight frames.
% 
% Appearance Time: Time at which a colony is first detected (minutes after
% experiment start)
%
% Growth Rate: Linear growth rate fitted to colony areas, using a linear
% regression method that is robust to outliers, see: fitlm, 'RobustOps'. Areas
% are only fitted after 8 frames after appearance to capture the linear
% portion. The gradient of the linear regression line is the growth rate of the
% colony
%
% Growth Time: As specified in Levin-Reisman et al., this is the time taken for
% a six-fold increase in area of a colony. This is taken to mean from the time
% of first appearance. A future improvement could be to calculate a moving
% average of this time over the entire experiment and average this as the
% 'growth time'.
% -----------------------------------------------------------------------

% parse arguments
options = struct('debug', 0, 'method', 1);
optionNames = fieldnames(options);

% if odd number of arguments: DirNames is specified
if round(nargin / 2) ~= nargin / 2
    DirNames = varargin(1);
    varargin = varargin(2:end);
else
    DirNames = pwd;
end

for pair = reshape(varargin, 2, [])
    inpName = lower(pair{1});
    if any(strcmp(inpName, optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name', inpName);
    end
end

if options.debug > 1
    disp('Debug has been set; you can type "close all" to close all figures');
end

if ischar(DirNames)
    DirNames = {DirNames};
end

GrowthRates = [];
AppTimes = [];

Levin_AppTimes = [];
Levin_Rate = [];
for DirName = DirNames
    DirName = char(DirName);
    disp(sprintf('Loading data from <%s>', DirName));
    datafile = fullfile(DirName, 'Results', 'VecArea.mat');
    data = load(datafile);
    data = data.VecArea;
    timefile = fullfile(DirName, 'Results', 'TimeAxis.mat');
    timeaxis = load(timefile);
    timeaxis = timeaxis.TimeAxis;

    %subplot(211);
    disp('Fitting areas to determine growth rates');
    for colony = data'
        if colony(end) == 0
        else
            fidx = find(colony > 0, 1, 'first');
            if fidx == 0
                fidx = 1
            end
            area = colony(fidx-1:end);
            area_lim = area(7:end);
            time = timeaxis(fidx-1:end)';
            time_lim = time(7:end);
            % calculate time for area to increase by 6-fold
            if options.method == 1
                lindex = find(area > area(2) * 6, 1, 'first');
                if length(lindex) > 0
                    t0 = time(2);
                    t1 = time(lindex);
                    tdiff = t1 - t0;
                    Levin_Rate = [Levin_Rate, tdiff];
                    Levin_AppTimes = [Levin_AppTimes, time(2)];
                end
            else
                Levin_Rate = [Levin_Rate, get_growth_times(time, area)];
                Levin_AppTimes = [Levin_AppTimes, time(2)];
            end
            if length(area) < 10 
            else
                % fit linear regression
                f2 = fitlm(time_lim, area_lim, 'RobustOpts', 'on');
                B = table2array(f2.Coefficients(:, 1));
                intercept_r = B(1);
                gradient_r = B(2);
                GrowthRates = [GrowthRates, gradient_r * 60];
                AppTimes = [AppTimes, time(2)];

                if options.debug > 0
                    figure;
                    hold on
                    dpoints = scatter(time, area, 70, ...
                        'MarkerFaceColor', [44 162 95] / 256, ...
                        'MarkerEdgeColor', [100 100 100] / 256, ...
                        'LineWidth', 1);
                    fitline_r = plot(time_lim, (gradient_r * time_lim) + intercept_r, ...
                        'LineWidth', 3', ...
                        'Color', 'black');
                    xlabel('Time (min)');
                    ylabel('Area');
                    legend([dpoints, ...
                        fitline_r], ...
                        'Data', ...
                        'Fit', ...
                        'Location', 'Best');
                    hold off
                end
            end
        end
    end
    if options.debug > 0
        input('Press "enter" when you are happy with the debugging output; this will close all figures');
    end
    close all;
end

disp('Plotting growth rates');
figure('units', 'inches', 'pos', [0 0 14 6]);
subplot(131);
[n_at, x_at] = hist(AppTimes);
h_at = bar(x_at, n_at, 1.0);
set(h_at, ...
    'facecolor', [228 26 28] / 256, ...
    'edgecolor', 'k');
xlabel('Appearance Time (min)');
ylabel('Frequency');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
xticks = get(gca, 'XTick');
set(gca, 'XTick', unique(round(xticks / 50) * 50));

ax = subplot(132);
[n_gr, x_gr] = hist(GrowthRates);
h_gr = bar(x_gr, n_gr, 1.0);
set(h_gr, ...
    'facecolor', [55 126 184] / 256, ...
    'edgecolor', 'k');

xlabel('Growth Rate (px^2 / h)');
ylabel('Frequency');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');

subplot(133);
sc = scatter(AppTimes, GrowthRates);
set(sc, ...
    'markerfacecolor', [77 175 74] / 256, ...
    'markeredgecolor', 'k');
ylabel('Growth Rate (px^2 / h)');
xlabel('Appearance Time (min)');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperPosition', [0 0 6 3]);
set(gcf, 'PaperOrientation', 'landscape');
msg = sprintf('Saved to %s', fullfile(pwd, 'out.pdf'));
disp(msg);

print('-dpdf', '-r0', 'out');
%nhist(GrowthRates, 'xlabel', 'Growth Rate');

disp('Plotting using Levin-Reisman et al. method');
% time taken for 6-fold increase
figure('units', 'inches', 'pos', [0 0 14 6]);
subplot(131);
h_at2 = bar(x_at, n_at, 1.0);
set(h_at2, ...
    'facecolor', [228 26 28] / 256, ...
    'edgecolor', 'k');
xlabel('Appearance Time (min)');
ylabel('Frequency');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
xticks = get(gca, 'XTick');
set(gca, 'XTick', unique(round(xticks / 50) * 50));

ax = subplot(132);
[n_lev, x_lev] = hist(Levin_Rate);
h_lev = bar(x_lev, n_lev, 1.0);
set(h_lev, ...
    'facecolor', [152 78 163] / 256, ...
    'edgecolor', 'k');

xlabel('Growth Time (min)');
ylabel('Frequency');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');

subplot(133);
sc = scatter(Levin_AppTimes, Levin_Rate);
set(sc, ...
    'markerfacecolor', [77 175 74] / 256, ...
    'markeredgecolor', 'k');
ylabel('Growth Time (min)');
xlabel('Appearance Time (min)');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperPosition', [0 0 6 3]);
set(gcf, 'PaperOrientation', 'landscape');

msg = sprintf('Saved to %s', fullfile(pwd, 'out2.pdf'));
disp(msg);

print('-dpdf', '-r0', 'out2');
end

%nhist(GrowthRates, 'xlabel', 'Growth Rate');

function growth_time = get_growth_times(time, area)
    fold_factor = 6;
    aindexes = [];
    lindexes = [];
    growth_times = [];
    for idx = 2:length(area)
        lindex = find(area > area(idx) * 6, 1, 'first');
        if length(lindex) > 0
            lindexes = [lindexes lindex];
            aindexes = [aindexes idx];
            t0 = time(idx);
            t1 = time(lindex);
            tdiff = t1 - t0;
            growth_times = [growth_times tdiff];
        end
    end
    growth_time = mean(growth_times);
end
