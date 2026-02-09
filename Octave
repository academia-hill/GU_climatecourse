clear; clc; close all;

%% 1. Startup and loading
% Configuration of data files and plot titles
datasets = {
    'EastAsia.txt',       3, 1, 'East Asia';       
    'Tibet.txt',          3, 1, 'Tibetan Plateau';           
    'China.txt',          3, 1, 'North Central China';          
    'West_mongolia.txt',  3, 1, 'West Mongolia';    
    'North_mongolia.txt', 3, 1, 'North Mongolia';  
    'Japan.txt',          4, 3, 'North Japan';     
};

% Define global x-axis limits
common_xlim = [1500 2005]; 

% Initialize figure window
figure('Name', 'Climate Analysis', 'NumberTitle', 'off', 'Position', [50, 50, 1400, 1000]);

num_files = size(datasets, 1);
h_leg_objects = []; 

%% 2. Loop to conduct research for all data sets, and analysis
for i = 1:num_files
    filename = datasets{i, 1};
    z_col_idx = datasets{i, 2};     
    header_skip = datasets{i, 3};   
    plot_title = datasets{i, 4};
    
    % --- Data Import ---
    if ~exist(filename, 'file'), continue; end
    try
        file_content = fileread(filename);
        file_content = strrep(file_content, ',', '.');
        temp_file = ['temp_' filename];
        fid = fopen(temp_file, 'w');
        fprintf(fid, '%s', file_content);
        fclose(fid);
        data = dlmread(temp_file, '', header_skip, 0);
        delete(temp_file);
    catch
        continue;
    end

    years = data(:, 1);
    if size(data, 2) < z_col_idx, continue; end
    temp_raw = data(:, z_col_idx); 

    % Outlier Removal (3-sigma threshold)
    mu = mean(temp_raw); sd = std(temp_raw);
    temp = temp_raw;
    temp(abs(temp_raw - mu) > 3 * sd) = median(temp_raw);

    % Baseline calculations (1720-1800)
    ref_idx = find(years >= 1720 & years <= 1800);
    if isempty(ref_idx)
        b_mean = mean(temp); b_2sigma = 2 * std(temp);
    else
        b_mean = mean(temp(ref_idx)); b_2sigma = 2 * std(temp(ref_idx));
    end
    threshold = b_mean + b_2sigma;

    % Trend analysis through Gaussian filters
    all_onsets = []; all_toes = [];
    search_start = find(years == 1801);
    if isempty(search_start), search_start = 1; end
    limit = length(years) - 30;
    
    if search_start <= limit
        for fw = 15:5:50
            sigma = fw / 6; k = -ceil(3*sigma):ceil(3*sigma);
            kernel = exp(-k.^2 / (2*sigma^2)); kernel = kernel / sum(kernel);
            sm = conv(temp, kernel, 'same') ./ conv(ones(size(temp)), kernel, 'same');
            
            % Regressive Change-Point Analysis (Onset)
            best_cp = NaN; min_rss = inf;
            for k = search_start:limit
                x2 = years(k:end); y2 = sm(k:end);
                X = [ones(length(x2),1), x2(:)]; b = X \ y2(:);
                if b(2) > 0 && (b(2)/sqrt((sum((y2-X*b).^2)/(length(x2)-2))/sum((x2-mean(x2)).^2))) > 1.65 ...
                   && sm(end) > sm(k) && sm(k) > b_mean
                    rss = sum((sm(1:k)-mean(sm(1:k))).^2) + sum((y2-X*b).^2);
                    if rss < min_rss, min_rss = rss; best_cp = years(k); end
                end
            end
            if ~isnan(best_cp), all_onsets = [all_onsets, best_cp]; end
            
            % Time of Emergence (ToE) detection
            idx_toe = find(all(sm(search_start:end) > threshold, 2), 1, 'first');
            if ~isempty(idx_toe), all_toes = [all_toes, years(search_start + idx_toe - 1)]; end
        end
    end
    med_onset = median(all_onsets);
    med_toe = median(all_toes);

    % --- PLOTTING ---
    ax = subplot(num_files, 1, i);
    hold on; box on;

    % Subplot lettering (a, b, c...)
    labels = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)'};
    text(-0.05, 1.12, labels{i}, 'Units', 'normalized', ...
         'FontSize', 14, 'FontWeight', 'bold', 'Clipping', 'off'); 

    set(gca, 'FontSize', 12);

    % Subplot positioning
    pos = get(ax, 'Position'); 
    set(ax, 'Position', [0.10, pos(2), 0.60, pos(4)]); 

    % Smoothing filters for visualization (15 and 50 years)
    fw = 15; sig = fw/6; k = exp(-(-3*sig:3*sig).^2/(2*sig^2)); k=k/sum(k);
    s15 = conv(temp, k, 'same') ./ conv(ones(size(temp)), k, 'same');
    fw = 50; sig = fw/6; k = exp(-(-3*sig:3*sig).^2/(2*sig^2)); k=k/sum(k);
    s50 = conv(temp, k, 'same') ./ conv(ones(size(temp)), k, 'same');

    % Graphical rendering
    if ~isnan(med_toe)
        fill([med_toe max(years) max(years) med_toe], [-4 -4 4 4], [1 0.9 0.9], 'EdgeColor', 'none');
    end
    h_base = fill([1720 1800 1800 1720], [b_mean-b_2sigma b_mean-b_2sigma threshold threshold], ...
                  [0.9 0.9 0.95], 'EdgeColor', 'none');
    h_raw = plot(years, temp_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    h_s15 = plot(years, s15, 'k-', 'LineWidth', 1.0);
    h_s50 = plot(years, s50, 'b-', 'LineWidth', 2.5);
    h_thr = plot(years, repmat(threshold, size(years)), 'r--', 'LineWidth', 1.5);

    % Markers for Onset and ToE
    if ~isnan(med_onset)
        h_onset = line([med_onset med_onset], [-4 4], 'Color', [0.9 0.4 0], 'LineWidth', 2.0);
        text(med_onset, 3.2, sprintf('%d', round(med_onset)), ...
             'Color', [0.9 0.4 0], 'FontWeight', 'bold', 'FontSize', 10);
    else
        h_onset = plot(NaN, NaN, 'Color', [0.9 0.4 0], 'LineWidth', 2.0);
    end
    
    if ~isnan(med_toe)
        h_toe = line([med_toe med_toe], [-4 4], 'Color', [0.6 0 0.8], 'LineWidth', 2.0, 'LineStyle', '--');
        text(med_toe, 2.5, sprintf('%d', round(med_toe)), ...
             'Color', [0.6 0 0.8], 'FontWeight', 'bold', 'FontSize', 10);
    else
        h_toe = plot(NaN, NaN, 'Color', [0.6 0 0.8], 'LineWidth', 2.0, 'LineStyle', '--');
    end

    % Legend objects capture (from first subplot)
    if i == 1
        h_leg_objects = [h_base, h_raw, h_s15, h_s50, h_thr, h_onset, h_toe];
    end

    title(plot_title, 'FontSize', 12, 'FontWeight', 'bold');
    ylim([-4 4]); 
    xlim(common_xlim);
    ylabel('Z-score', 'FontSize', 12);
    
    if i < num_files
        set(gca, 'XTickLabel', []); 
    else
        xlabel('Year', 'FontSize', 12); 
    end
end

%% 3. Legend Configuration
if ~isempty(h_leg_objects)
    lgd = legend(h_leg_objects, ...
           'Baseline (1720-1800)', 'Annual Z-scores', '15-yr Trend', ...
           '50-yr Trend', 'Threshold (+2\sigma)', 'Year of Onset', 'Time of Emergence');
    
    set(lgd, 'FontSize', 11, 'Location', 'best'); 
    drawnow; 
    pos_tight = get(lgd, 'Position');
    set(lgd, 'Position', [0.72, 0.5 - (pos_tight(4)/2), pos_tight(3), pos_tight(4)]);
end

hold off;
