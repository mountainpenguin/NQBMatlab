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
%            merge_method (optional)
%               Defaults to 1: take region before the merge point
%               If set to 2: ignore any colonies that merge
%               If set to 3: ignore merge events
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

% switch off warnings
warning('off', 'all');

% parse arguments
options = struct('debug', 0, 'method', 1, 'merge_method', 1);
optionNames = fieldnames(options);

% if odd number of arguments: DirNames is specified
if round(nargin / 2) ~= nargin / 2
    DirNames = varargin(1);
    DirNames = DirNames{1};
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

if ischar(DirNames)
    DirNames = {DirNames};
end

disp('Directories to be processed:')
for k = 1:length(DirNames)
    fprintf(' %s\n', char(DirNames{k}));
end
fprintf('\n');
disp('Options:')
disp(' Method:')
if options.method == 1
    disp('  [1] time taken for six-fold increase in area from time of first appearance');
    area_min = inputdlg('Area to start from is?');
    area_max = inputdlg('Area to end at?');
    area_min = str2num(area_min{1});
    area_max = str2num(area_max{1});
elseif options.method == 2
    disp('  [2] averaged time taken for six-fold increase in area for all timepoints');
elseif options.method == 3
    disp('  [3] not implemented');
else
    fprintf('  [%i] unknown method\n', options.method);
end
disp(' Merge Method:')
if options.merge_method == 2
    disp('  [2] ignore all merged colonies');
elseif options.merge_method == 1
    disp('  [1] use region before merge event');
elseif options.merge_method == 3
    disp('  [3] ignore any merge events');
else
    fprintf('  [%i] unknown merge method\n', options.merge_method);
end
disp('  Debug:')
if options.debug == 0
    disp('  [0] No debugging');
else
    disp('  [1] Debug growth rate fitting');
end
fprintf('\n');

if options.debug > 1
    disp('Debug has been set; you can type "close all" to close all figures');
end

GrowthRates = [];
AppTimes = [];

Levin_AppTimes = [];
Levin_Rate = [];

GrowthAreas = {};
GrowthTimes = {};

Levin_GrowthAreas = {};
Levin_GrowthTimes = {};

Colony_Nums = [];
Levin_Colony_Nums = [];

for k = 1:length(DirNames)
    DirName = char(DirNames{k});
    disp(sprintf('Loading data from <%s>', DirName));
    datafile = fullfile(DirName, 'Results', 'VecArea.mat');
    data = load(datafile);
    data = data.VecArea;

    merged = getMergedColonies(fullfile(DirName, 'Results'));
    excluded = load(fullfile(DirName, 'Results', 'ExcludedBacteria.txt'));
    notclosetoborder = FindColoniesInWorkingArea(DirName);

    timefile = fullfile(DirName, 'Results', 'TimeAxis.mat');
    timeaxis = load(timefile);
    timeaxis = timeaxis.TimeAxis;

    %subplot(211);
    disp('Fitting areas to determine growth rates');
    for colony_num = 1:size(data, 1)
        %perc = sprintf('%.2f%%', 100 * colony_num / size(data, 1));
        %fprintf('Dataset %3i of %3i: Colony %3i of %3i (%7s) - ', k, length(DirNames), colony_num, size(data, 1), perc);
        fprintf('Dataset %3i of %3i: Colony %3i of %3i (%7s) - ', k, length(DirNames), colony_num, size(data, 1));
        border_status = find(notclosetoborder == colony_num, 1, 'first');
        exclude_status = find(excluded == colony_num, 1, 'first');
        if options.merge_method == 2
            merge_status = merged(colony_num);
        elseif options.merge_method == 1
            merge_status = 0;
        elseif options.merge_method == 3
            merge_status = 0;
        end
        colony = data(colony_num, :);

        if (colony(end) == 0)
            fprintf('rejected (does not exist in final frame)\n');
        elseif (merge_status ~= 0) 
            fprintf('rejected (colony merges)\n');
        elseif (isempty(border_status)) 
            fprintf('rejected (colony is on plate border)\n');
        elseif (~isempty(exclude_status))
            fprintf('rejected (colony is on exclusion list)\n');
        else
            fidx = find(colony > 0, 1, 'first');
            if fidx == 0
                fidx = 1
            end
            if options.merge_method == 1
                merge_val = merged(colony_num);
                if merge_val ~= 0
                    area = colony(fidx-1:merge_val-1);
                    time = timeaxis(fidx-1:merge_val-1);
                else
                    area = colony(fidx-1:end);
                    time = timeaxis(fidx-1:end);
                end
            else
                area = colony(fidx-1:end)';
                %area_lim = area(7:end);
                time = timeaxis(fidx-1:end)';
                %time_lim = time(7:end);
            end
            if length(area) >= 10 
%                % fit linear regression
%                f2 = fitlm(time_lim, area_lim, 'RobustOpts', 'on');
%                B = table2array(f2.Coefficients(:, 1));
%                intercept_r = B(1);
%                gradient_r = B(2);
%                GrowthRates = [GrowthRates, gradient_r * 60];  % px^2 per hour
%                AppTimes = [AppTimes, time(2)];
                
                % fit polynomial
                %model = @(B, t)B(1) + B(2).*t + B(3).*t.^2;
                %B0 = [1 1 1];  % initial estimate

                % fit logistic function
                model = @(B,t) B(1) + B(2) ./ (1 + exp(-B(3) .* (t-B(4))));
                B0 = [0 1000 0.001 300];
                coeffnames = {'yoffset', 'maxY', 'gamma', 'alpha'};
                f3 = fitnlm(time, area, model, B0, 'CoefficientNames', coeffnames);

%                model = @(B, t) B(1) + B(2) ./ (B(5) + B(4) .* exp(-B(3) .* (t - B(6))));
%                B0 = [0 1000 0.001 1 1 300];
%                f3 = fitnlm(time, area, model, B0);

                coeff = table2array(f3.Coefficients(:, 1));
                rmse = f3.RMSE;
                colony_discard = 0;
                if (rmse < 30) && (coeff(3) > 0) && (coeff(3) < 0.02)
                    %% midpoint of linear portion of logistic curve
                    midX = coeff(4);
                    mididx = find(abs(time - midX) == min(abs(time - midX)), 1, 'first');
                    if (mididx > 3)
    %                    % use data region
    %                    if time(mididx) < midX
    %                        linear_time = time(mididx - 2:mididx + 3);
    %                        linear_area = area(mididx - 2:mididx + 3);
    %                    else
    %                        linear_time = time(mididx - 3:mididx + 2);
    %                        linear_area = area(mididx - 3:mididx + 2);
    %                    end

                        % use modelled region
                        linear_time = time(mididx - 3:mididx);
                        linear_area = model(coeff, time);
                        linear_area = linear_area(mididx - 3:mididx);
                        % fit linear regression to linear modelled region
                        f4 = fitlm(linear_time, linear_area, 'RobustOpts', 'on');
                        lin_coeff = table2array(f4.Coefficients(:, 1));
                        GrowthRates = [GrowthRates, lin_coeff(2) * 60]; 
                        GrowthAreas{length(GrowthAreas) + 1} = area;
                        GrowthTimes{length(GrowthTimes) + 1} = time;
                        Colony_Nums = [Colony_Nums, colony_num];

                        %GrowthRates = [GrowthRates, coeff(3) * 60];  % hr^{-1}
                        AppTimes = [AppTimes, time(2)];
                        colony_discard = 1;

                        fprintf('accepted\n');
                    else
                        fprintf('rejected (linear region too early)\n');
                    end
                else
                    fprintf('rejected (poor fit)\n');
                end

                if (options.debug > 0) && (colony_discard == 1)% & length(GrowthRates) == 8
                    %f3
%                    grad = gradient(gradient(model(B, time)));
%                    lims = 0.3;
%                    linear_x1 = find(abs(grad) < lims, 1, 'first');
%                    linear_x2 = linear_x1 + find(abs(grad(linear_x1:end)) > lims, 1, 'first') - 1;
%                    linear_portion = model(B, time);
%                    linear_portion = linear_portion(linear_x1:linear_x2);
%                    linear_time = time(linear_x1:linear_x2);
%                    %linear = plot(linear_time, linear_portion, 'LineWidth', 5);
                    figure;
                    hold on
                    dpoints = scatter(time, area, 70, ...
                        'MarkerFaceColor', [44 162 95] / 256, ...
                        'MarkerEdgeColor', [100 100 100] / 256, ...
                        'LineWidth', 1);
%                    fitline_r = plot(time, (gradient_r * time) + intercept_r, ...
%                        'LineWidth', 3, ...
%                        'Color', 'black');
                    fitline_r2 = plot(time, model(coeff, time), ...
                        'LineWidth', 3, ...
                        'Color', 'red');

                    f4_model = @(B, x) B(1) + B(2) .* x;
                    ylim = get(gca, 'YLim');
                    plot(time, f4_model(lin_coeff, time), '--', 'Color', 'k', 'LineWidth', 2);
                    set(gca, 'YLim', ylim);
 
                    xlabel('Time (min)');
                    ylabel('Area');
                    t = sprintf('Colony %i', colony_num);
                    title(t);
                    legend([dpoints, ...
                        fitline_r2], ...
                        'Data', ...
                        'Fit2', ...
                        'Location', 'Best');
                    hold off
%                    assignin('base', 'time', time');
%                    assignin('base', 'area', area);
                    input('next? ')
                end
            else
                fprintf('rejected (not enough datapoints)\n');
            end

            % calculate time for area to increase by 6-fold
            if options.method == 1
                index1 = find(area > area_min, 1, 'first');
                index2 = find(area > area_max, 1, 'first');
                if (length(index1) > 0) + (length(index2) > 0) == 2
                    t0 = time(index1);
                    t1 = time(index2);
                    tdiff = t1 - t0;
                    Levin_Rate = [Levin_Rate, tdiff];
                    Levin_AppTimes = [Levin_AppTimes, time(2)];
                    Levin_GrowthAreas{length(Levin_GrowthAreas) + 1} = area;
                    Levin_GrowthTimes{length(Levin_GrowthTimes) + 1} = time;
                    Levin_Colony_Nums = [Levin_Colony_Nums, colony_num];
                end
            elseif options.method == 2
                Levin_Rate = [Levin_Rate, get_growth_times(time, area)];
                Levin_AppTimes = [Levin_AppTimes, time(2)];
                Levin_GrowthAreas{length(Levin_GrowthAreas) + 1} = area;
                Levin_GrowthTimes{length(Levin_GrowthTimes) + 1} = time;
                Levin_Colony_Nums = [Levin_Colony_Nums, colony_num];
            elseif options.method == 3
                Levin_Rate = [Levin_Rate, get_linear_growth_times(time, area, f3)];
                Levin_AppTimes = [Levin_AppTimes, time(2)];
                Levin_GrowthAreas{length(Levin_GrowthAreas) + 1} = area;
                Levin_GrowthTimes{length(Levin_GrowthTimes) + 1} = time;
                Levin_Colony_Nums = [Levin_Colony_Nums, colony_num];
            end
        end
    end
    fprintf('\n');
    if options.debug > 0
        input('Press "enter" when you are happy with the debugging output; this will close all figures');
        close all;
    end
end

disp('Plotting area graph for colonies used for each plot');
figure('Units', 'Inches');
subplot(1, 2, 1);
cc = lines(length(GrowthAreas));
hold on;
for k = 1:length(GrowthAreas)
    ga = GrowthAreas{k};
    gt = GrowthTimes{k};
    plot(gt, ga, 'color', cc(k, :), 'LineWidth', 3);
end
hold off;
title('Colonies used for determining growth rates');
xlabel('Time (min)');
ylabel('Area (px^2)');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);
set(gca, 'XTick', unique(round(timeaxis / 500) * 500));

subplot(1, 2, 2);
cc = lines(length(Levin_GrowthAreas));
hold on;
for k = 1:length(Levin_GrowthAreas)
    ga = Levin_GrowthAreas{k};
    gt = Levin_GrowthTimes{k};
    plot(gt, ga, 'color', cc(k, :), 'LineWidth', 3);
end
hold off;
title('Colonies used for determining growth times');
xlabel('Time (min)');
ylabel('Area (px^2)');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);
set(gca, 'XTick', unique(round(timeaxis / 500) * 500));

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperSize', [pos(3), pos(4)]);

msg = sprintf('Saved to %s', fullfile(pwd, 'areas.pdf'));
disp(msg);

print('-dpdf', '-r0', 'areas');

disp('Plotting growth rates');
figure('Units', 'Inches');
subplot(2, 2, 1);
[n_at, x_at] = hist(AppTimes);
h_at = bar(x_at, n_at, 1.0);
set(h_at, ...
    'facecolor', [228 26 28] / 256, ...
    'edgecolor', 'k');
xlabel('Appearance Time (min)');
ylabel('Count');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);
set(gca, 'XTick', unique(round(timeaxis / 500) * 500));

ax = subplot(2, 2, 2);
[n_gr, x_gr] = hist(GrowthRates);
h_gr = bar(x_gr, n_gr, 1.0);
set(h_gr, ...
    'facecolor', [55 126 184] / 256, ...
    'edgecolor', 'k');

xlabel('Growth Rate (px^2 / h)');
ylabel('Count');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
%set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);

subplot(2, 2, [3, 4]);
sc = scatter(AppTimes, GrowthRates);
set(sc, ...
    'markerfacecolor', [77 175 74] / 256, ...
    'markeredgecolor', 'k');
ylabel('Growth Rate (px^2 / h)');
xlabel('Appearance Time (min)');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);
set(gca, 'XTick', unique(round(timeaxis / 500) * 500));
t = sprintf('n = %i', length(AppTimes));
title(t);

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperSize', [pos(3), pos(4)]);
msg = sprintf('Saved to %s', fullfile(pwd, 'out.pdf'));
disp(msg);

print('-dpdf', '-r0', 'out');
%nhist(GrowthRates, 'xlabel', 'Growth Rate');

disp('Plotting using Levin-Reisman et al. method');
% time taken for 6-fold increase
figure('Units', 'Inches');
subplot(2, 2, 1);
[n_at, x_at] = hist(Levin_AppTimes);
h_at2 = bar(x_at, n_at, 1.0);
set(h_at2, ...
    'facecolor', [228 26 28] / 256, ...
    'edgecolor', 'k');
xlabel('Appearance Time (min)');
ylabel('Count');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);
set(gca, 'XTick', unique(round(timeaxis / 500) * 500));

ax = subplot(2, 2, 2);
[n_lev, x_lev] = hist(Levin_Rate);
h_lev = bar(x_lev, n_lev, 1.0);
set(h_lev, ...
    'facecolor', [152 78 163] / 256, ...
    'edgecolor', 'k');

xlabel('Growth Time (min)');
ylabel('Count');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
xticks = get(gca, 'XTick');
set(gca, 'XTick', unique(round(xticks / 50) * 50));

subplot(2, 2, [3, 4]);
sc = scatter(Levin_AppTimes, Levin_Rate);
set(sc, ...
    'markerfacecolor', [77 175 74] / 256, ...
    'markeredgecolor', 'k');

ylabel('Growth Time (min)');
xlabel('Appearance Time (min)');
set(gca, 'box', 'off');
set(gca, 'Color', 'none');
set(gca, 'XLim', [min(timeaxis), max(timeaxis)]);
set(gca, 'XTick', unique(round(timeaxis / 500) * 500));
t = sprintf('n = %i', length(Levin_AppTimes));
title(t);

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperSize', [pos(3), pos(4)]);

msg = sprintf('Saved to %s', fullfile(pwd, 'out2.pdf'));
disp(msg);

print('-dpdf', '-r0', 'out2');

% save data to Excel
sheet1 = table([1:length(AppTimes)]', AppTimes', GrowthRates', 'VariableNames', {'Colony', 'Appearance_Time', 'Growth_Rate'});
writetable(sheet1, 'output_data.xlsx', 'Sheet', 1);
sheet2 = table([1:length(Levin_AppTimes)]', Levin_AppTimes', Levin_Rate', 'VariableNames', {'Colony', 'Appearance_Time', 'Time_to_6fold_increase'}); 
writetable(sheet2, 'output_data.xlsx', 'Sheet', 2);

mega_table = zeros(length(timeaxis), length(GrowthTimes));
mega_table(:, 1) = timeaxis;
for k = 1:length(GrowthAreas)
    this_times = GrowthTimes{k};
    this_areas = GrowthAreas{k}';

    first_time = find(timeaxis == this_times(1), 1, 'first');
    for i = 1:length(this_times)
        mega_table(i + first_time - 1, k + 1) = this_areas(i);
    end
end
column_names = {'Time'};
for a = 1:length(Colony_Nums)
    column_names{a + 1} = sprintf('Colony %d', Colony_Nums(a));
end
column_names = table(column_names);
mega_table = table(mega_table);
writetable(mega_table, 'output_data.xlsx', 'Sheet', 3);
writetable(column_names, 'output_data.xlsx', 'Sheet', 4);

mega_table = zeros(length(timeaxis), length(Levin_GrowthTimes));
mega_table(:, 1) = timeaxis;
for k = 1:length(Levin_GrowthAreas)
    this_times = Levin_GrowthTimes{k};
    this_areas = Levin_GrowthAreas{k}';

    first_time = find(timeaxis == this_times(1), 1, 'first');
    for i = 1:length(this_times)
        mega_table(i + first_time - 1, k + 1) = this_areas(i);
    end
end
column_names = {'Time'};
for a = 1:length(Levin_Colony_Nums)
    column_names{a + 1} = sprintf('Colony %d', Levin_Colony_Nums(a));
end
column_names = table(column_names);
mega_table = table(mega_table);
writetable(mega_table, 'output_data.xlsx', 'Sheet', 5);
writetable(column_names, 'output_data.xlsx', 'Sheet', 6);

end

%nhist(GrowthRates, 'xlabel', 'Growth Rate');

function growth_time = get_growth_times(time, area)
    fold_factor = 6;
    growth_times = [];
    for idx = 2:length(area)
        lindex = find(area > area(idx) * 6, 1, 'first');
        if length(lindex) > 0
            t0 = time(idx);
            t1 = time(lindex);
            tdiff = t1 - t0;
            growth_times = [growth_times tdiff];
        end
    end
    growth_time = mean(growth_times);
end

function growth_time = get_linear_growth_times(time, area, nlmmodel)
    fold_factor = 6;
    coeff = table2array(nlmmodel.Coefficients(:, 1));
    yoffset = coeff(1)  % Y-offset
    maxY = coeff(2)  % max Y
    gamma = coeff(3)  % steepness
    alpha = coeff(4)  % midpoint X

    error('stop here');
    growth_time = 0;
end
