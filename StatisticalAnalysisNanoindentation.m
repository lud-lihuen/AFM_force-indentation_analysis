% Analisis estadistico de resultados de nanoindentacion
close all
clear
clc

%% Leer resultados - Modulo de Young

% Grupo control:
E_control = dlmread('E_control'); % modulo de Young en kPa
E_control_mean = mean(E_control); % media
E_control_deviation = std(E_control); % desviacion estandar

% Grupo modificado:
E_modified = dlmread('E_modified'); % modulo de Young en kPa
E_modified_mean = mean(E_modified); % media
E_modified_deviation = std(E_modified); % desviacion estandar

% Resultados:
disp('--- Resultados para módulo de Young ---');
disp(' ');
disp(['Grupo control: ' num2str(length(E_control)) ' curvas de fuerza analizadas.']);
disp(['Módulo de Young = ' num2str(E_control_mean,'%.2f') ' +/- ' num2str(E_control_deviation,'%.2f') ' kPa']);
disp(' ');
disp(['Grupo modificado: ' num2str(length(E_modified)) ' curvas de fuerza analizadas.']);
disp(['Módulo de Young = ' num2str(E_modified_mean,'%.2f') ' +/- ' num2str(E_modified_deviation,'%.2f') ' kPa']);
disp(' ');

%% Leer resultados - Trabajo de Adhesion

% Grupo control:
A_control = dlmread('A_control'); % trabajo de adhesion en J/m2
A_control_mean = mean(A_control); % media
A_control_deviation = std(A_control); % desviacion estandar

% Grupo modificado:
A_modified = dlmread('A_modified'); % trabajo de adhesion en J/m2
A_modified_mean = mean(A_modified); % media
A_modified_deviation = std(A_modified); % desviacion estandar

% Resultados:
disp('--- Resultados para trabajo de adhesión ---');
disp(' ');
disp(['Grupo control: ' num2str(length(A_control)) ' curvas de fuerza analizadas.']);
disp(['Trabajo de Adhesión = ' num2str(A_control_mean,'%.2e') ' +/- ' num2str(A_control_deviation,'%.2e') ' J/m2']);
disp(' ');
disp(['Grupo modificado: ' num2str(length(A_modified)) ' curvas de fuerza analizadas.']);
disp(['Trabajo de Adhesión = ' num2str(A_modified_mean,'%.2e') ' +/- ' num2str(A_modified_deviation,'%.2e') ' J/m2']);
disp(' ');

%% Histograma - Modulo de Young

figure;
histogram(E_control, 100, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, 50]);
hold on;
histogram(E_modified, 100, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, 50]);
hold off;
xlabel("Young's Modulus (kPa)");
ylabel("Force curves");
legend("Control", "Modified");
title("Young's Modulus");
grid on;

%% Histograma - Trabajo de Adhesion

figure;
histogram(A_control, 100, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, 0.01]);
hold on;
histogram(A_modified, 100, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, 0.01]);
hold off;
xlabel("Work of Adhesion (J/m2)");
ylabel("Force curves");
legend("Control", "Modified");
title("Work of Adhesion");
grid on;

%% Boxplot - Modulo de Young

% Preparo datos para hacer el boxplot:
E_data = [E_control; E_modified];
E_boxplot = [repmat({'Control'}, length(E_control), 1); repmat({'Modified'}, length(E_modified), 1)];

% Grafico el boxplot:
figure;
boxplot(E_data, E_boxplot, 'Labels', {'Control', 'Modified'});
ylim([0, 50]);
title("Young's Modulus (kPa)");
grid on;

%% Boxplot - Trabajo de Adhesion

% Preparo datos para hacer el boxplot:
A_data = [A_control; A_modified];
A_boxplot = [repmat({'Control'}, length(A_control), 1); repmat({'Modified'}, length(A_modified), 1)];

% Grafico el boxplot:
figure;
boxplot(A_data, A_boxplot, 'Labels', {'Control', 'Modified'});
ylim([0, 0.01]);
title("Work of Adhesion (J/m2)");
grid on;

%% Test de Student - Modulo de Young

[hE, pE] = ttest2(E_control, E_modified);

disp('--- Resultados del test de Student para módulo de Young ---');
disp(' ');
if hE == 1
  disp('Las muestras son significativamente diferentes.');
end
if hE == 0
  disp('Las muestras NO son significativamente diferentes.');
end
disp(['p-valor = ' num2str(pE)]);
disp(' ');

%% Test de Student - Trabajo de Adhesion

[hA, pA] = ttest2(A_control, A_modified);

disp('--- Resultados del test de Student para trabajo de adhesión ---');
disp(' ');
if hE == 1
  disp('Las muestras son significativamente diferentes.');
end
if hE == 0
  disp('Las muestras NO son significativamente diferentes.');
end
disp(['p-valor = ' num2str(pA)]);
disp(' ');