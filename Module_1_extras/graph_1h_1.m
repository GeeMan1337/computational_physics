%graph for question 1h_1
graph_data = readtable('E:\computational_physics\Module_1_out\graph_data_1h_1.dat');
x = table2array(graph_data(:, 1));
y = table2array(graph_data(:, 3));

figure;
grid on;
hold on;

% Scatter plot
scatter(x, y, 'red', 'Marker', '.');
xlabel('Bin Centre of Sums');
ylabel('Frequency');
title('10^4 sums of 10^4 random numbers in [0,1] with bin size 0.5');
legend('Data points for 10^4 sums');

hold off;