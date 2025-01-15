%graph for question 1h_7
graph_data_1 = readtable('E:\computational_physics\Module_1_out\graph_data_1h_4.dat');
graph_data_2 = readtable('E:\computational_physics\Module_1_out\graph_data_1h_5.dat');

x1 = table2array(graph_data_1(:, 1));
y1 = table2array(graph_data_1(:, 5));

x2 = table2array(graph_data_2(:, 1));
y2 = table2array(graph_data_2(:, 5));

gauss_fit_1 = fit(x1, y1, 'gauss1');
gauss_fit_2 = fit(x2, y2, 'gauss1');


figure;
grid on;
hold on;

% Scatter plot
scatter(x1, y1, 'blue', 'Marker', 'o');
scatter(x2, y2, 'red', 'Marker','x')

xlabel('Bin Centre of Sums');
ylabel('Normalized Frequency');
title('10^4 random number sums in [-1,1] with bin size 2 (normalised)');
legend('Data points for 10^4 sums','Data points for 10^5 sums');

mean_1 = gauss_fit_1.b1;
std_1 = gauss_fit_1.c1/sqrt(2);

mean_2 = gauss_fit_2.b1;
std_2 = gauss_fit_2.c1/sqrt(2);

text(max(x1)-90, max(y1)*0.8, {['Mean for 10^4 sums= ', num2str(mean_1, '%.2f')],...
    ['Std for 10^4 sums= ', num2str(std_1, '%.2f')],...
    ['Mean for 10^5 sums= ', num2str(mean_2, '%.2f')],...
    ['Std for 10^5 sums= ', num2str(std_2, '%.2f')]});

hold off;