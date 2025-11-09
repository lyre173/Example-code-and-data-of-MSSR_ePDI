% clc;
% clear;

% Assume intensityVals and intensityVals_mssr already exist in workspace
% Calculate percentage change relative to baseline (first 155 samples)
data_pdi = 100 .* (intensityVals - mean(intensityVals(1:155))) ./ mean(intensityVals(1:155));
data_mssrpdi = 100 .* (intensityVals_mssr - mean(intensityVals_mssr(1:155))) ./ mean(intensityVals_mssr(1:155));

% ===== Group data =====
a = [];
for i = 1:5
    a = cat(1, a, data_mssrpdi((155*(2*i-2)+1):(155*(2*i-1))));
end

b = [];
for i = 1:5
    b = cat(1, b, data_mssrpdi((155*(2*i-1)+1):(155*(2*i))));
end

% Plot data
figure;
plot(a);
hold on;
plot(b);
title('MSSR Data');

% Run analysis
analyze_ttest(a, b, 2); % Threshold set to 2%

% ===== PDI Data =====
a = [];
for i = 1:5
    a = cat(1, a, data_pdi((155*(2*i-2)+1):(155*(2*i-1))));
end

b = [];
for i = 1:5
    b = cat(1, b, data_pdi((155*(2*i-1)+1):(155*(2*i))));
end

figure;
plot(a);
hold on;
plot(b);
title('PDI Data');

% Run analysis
analyze_ttest(a, b, 2); % Threshold set to 2%


% ===== Function: automatic paired t-test + practical significance assessment =====
function analyze_ttest(data1, data2, threshold)
    % Paired t-test
    [h, p, ci, stats] = ttest(data1, data2);
    
    % Mean difference
    meanDiff = mean(data1 - data2);
    
    % Compute Cohen's d
    pooledStd = sqrt((std(data1)^2 + std(data2)^2) / 2);
    cohen_d = meanDiff / pooledStd;
    
    % Display results
    fprintf('\n==== Paired t-test Results ====\n');
    fprintf('h = %d (1 indicates significant difference, 0 indicates not significant)\n', h);
    fprintf('p-value = %.6f\n', p);
    fprintf('t-statistic = %.4f\n', stats.tstat);
    fprintf('Degrees of freedom = %d\n', stats.df);
    fprintf('95%% Confidence interval = [%.4f, %.4f]\n', ci(1), ci(2));
    fprintf('Mean difference = %.4f %%\n', meanDiff);
    fprintf('Cohen''s d = %.4f\n', cohen_d);
    
    % Dual criteria evaluation (statistical significance + practical relevance)
    if h == 1 && abs(meanDiff) >= threshold && abs(cohen_d) >= 0.5
        fprintf('Conclusion: Statistically significant and practically meaningful (difference >= %.2f%% and effect size >= 0.5)\n', threshold);
    elseif h == 1
        fprintf('Conclusion: Statistically significant but practically insufficient (difference or effect size too small)\n');
    else
        fprintf('Conclusion: No significant difference\n');
    end
end
