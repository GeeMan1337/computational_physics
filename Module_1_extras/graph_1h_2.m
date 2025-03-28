%graph for question 1h_2
graph_data = readtable('E:\computational_physics\Module_1_out\graph_data_1h_1.dat');
x = table2array(graph_data(:, 1));
y = table2array(graph_data(:, 7));
gauss_fit = fit(x, y, 'gauss1');

figure;
grid on;
hold on;

% Scatter plot
scatter(x, y, 'red', 'Marker', '.');

% Plotting Gaussian
plot(gauss_fit, 'green');

xlabel('Bin Centre of Sums');
ylabel('Normalized Frequency');
title('10^4 sums of 10^4 random numbers in [0,1] with bin size 0.5 (normalised)');
legend('Data points for 10^4 sums','Gaussian fit');

mean = gauss_fit.b1;
std = gauss_fit.c1/sqrt(2);

text(max(x)-50, max(y)*0.8, ['Mean = ', num2str(mean, '%.2f'), ', Std = ', num2str(std, '%.2f')]);

hold off;