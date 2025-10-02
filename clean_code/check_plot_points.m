coarse_range = 60:5:120;
fine_range = (105 - 5):0.1:(105 + 5);
plot_phases = unique([coarse_range, fine_range]);
fprintf('Plotting points: %d\n', length(plot_phases));
fprintf('Range: %.1f to %.1f degrees\n', min(plot_phases), max(plot_phases));