% Analyze phase distribution
coarse_range = 60:5:120;
fine_range = (105 - 5):0.1:(105 + 5);
ultra_range = (105 - 0.1):0.001:(105 + 0.1);
all_phases = [coarse_range, fine_range, ultra_range];
unique_phases = unique(all_phases);

fprintf('Phase distribution analysis:\n');
fprintf('Total unique phases: %d\n', length(unique_phases));
fprintf('Min phase: %.3f, Max phase: %.3f\n', min(unique_phases), max(unique_phases));

% Count distribution
n_60_100 = sum(unique_phases >= 60 & unique_phases < 100);
n_100_110 = sum(unique_phases >= 100 & unique_phases < 110);
n_104_106 = sum(unique_phases >= 104 & unique_phases < 106);
n_104_9_105_1 = sum(unique_phases >= 104.9 & unique_phases <= 105.1);

fprintf('\nPoint distribution:\n');
fprintf('  60-100°: %d points (%.1f%% of total)\n', n_60_100, 100*n_60_100/length(unique_phases));
fprintf('  100-110°: %d points (%.1f%% of total)\n', n_100_110, 100*n_100_110/length(unique_phases));
fprintf('  104-106°: %d points (%.1f%% of total)\n', n_104_106, 100*n_104_106/length(unique_phases));
fprintf('  104.9-105.1°: %d points (%.1f%% of total)\n', n_104_9_105_1, 100*n_104_9_105_1/length(unique_phases));

fprintf('\nThis explains the flat top - most points are concentrated near 105°!\n');