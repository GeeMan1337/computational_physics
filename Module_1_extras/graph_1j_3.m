%graph for question 1j_3
graph_data = readtable('E:\computational_physics\Module_1_out\graph_data_1i_5.dat');

x = table2array(graph_data(:, 1));
y = table2array(graph_data(:, 5));
gauss_fit = fit(x, y, 'gauss1');

figure;
grid on;
hold on;

% Scatter plot
scatter(x, y, 'blue', 'Marker', '.');

% Plotting Gaussian
plot(gauss_fit, 'red');

xlabel('Bin Centre of Sums');
ylabel('Normalized Frequency');
title('10^4 sums of 10^4 random numbers (either -1 or 1) with bin size 10 (normalised)');
legend('Data points for 10^4 sums','Gaussian fit');

mean = gauss_fit.b1;
std = gauss_fit.c1/sqrt(2);

text(max(x)-200, max(y)*0.8, {['Mean = ', num2str(mean, '%.2f'), ', Std = ', num2str(std, '%.2f')]...
    ['a1 = ',num2str(gauss_fit.a1)]});


hold off;