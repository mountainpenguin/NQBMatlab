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
for k = 1:length(DirNames)
    DirName = char(DirNames{k});
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

                        %GrowthRates = [GrowthRates, coeff(3) * 60];  % hr^{-1}
                        AppTimes = [AppTimes, time(2)];
                        colony_discard = 1;
                    end
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
            end

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
            elseif options.method == 2
                Levin_Rate = [Levin_Rate, get_growth_times(time, area)];
                Levin_AppTimes = [Levin_AppTimes, time(2)];
            end
        end
    end
    if options.debug > 0
        input('Press "enter" when you are happy with the debugging output; this will close all figures');
        close all;
    end
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
