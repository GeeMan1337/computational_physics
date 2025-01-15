%graph for question 1j_1
graph_data_1 = readtable('E:\computational_physics\Module_1_out\graph_data_1i_2.dat');
graph_data_2 = readtable('E:\computational_physics\Module_1_out\graph_data_1i_3.dat');

x1 = table2array(graph_data_1(:, 1));
y1 = table2array(graph_data_1(:, 7));

x2 = table2array(graph_data_2(:, 1));
y2 = table2array(graph_data_2(:, 7));

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
title('10^4 sums of 10^4 random numbers (either -1 or 1) with bin size 2 (normalised)');
legend('Data points for even bins','Data points for odd bins');

mean_1 = gauss_fit_1.b1;
std_1 = gauss_fit_1.c1/sqrt(2);
mean_2 = gauss_fit_2.b1;
std_2 = gauss_fit_2.c1/sqrt(2);

function_1 = @(x) gauss_fit_1.a1*exp(-((x-gauss_fit_1.b1)/gauss_fit_1.c1).^2);
function_2 = @(x) gauss_fit_2.a1*exp(-((x-gauss_fit_2.b1)/gauss_fit_2.c1).^2);

area_1 = integral(function_1, -inf, inf);
area_2 = integral(function_2, -inf, inf);

text(max(x1)-180, max(y1)*0.8, {['Mean for even bins= ', num2str(mean_1, '%.2f')],...
    ['Std for even bins bins= ', num2str(std_1, '%.2f')],...
    ['Mean for odd bins= ', num2str(mean_2, '%.2f')],...
    ['Std for odd bins= ', num2str(std_2, '%.2f')]});

text(max(x1)-180, max(y1)*0.5,{['area 1 is= ',num2str(area_1)],['area 2 is= ',num2str(area_2)]})

hold off;