% Analisis estadistico de resultados de mapas QNM
close all
clear
clc

%% Leer resultados - Modulo de Young

% Grupo control:
nE_control = 1; % numero de histogramas procesados (modificar)
E_control = readmatrix('E_control.txt');
E_control = repelem(E_control(:, 1), ceil((E_control(:, 2)/(100*nE_control))*size(E_control, 1))); % modulo de Young en kPa
E_control_mean = mean(E_control); % media
E_control_deviation = std(E_control); % desviacion estandar

% Grupo modificado:
nE_modified = 1; % numero de histogramas procesados (modificar)
E_modified = readmatrix('E_modified.txt');
E_modified = repelem(E_modified(:, 1), ceil((E_modified(:, 2)/(100*nE_modified))*size(E_modified, 1))); % modulo de Young en kPa
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

%% Leer resultados - Fuerza de Adhesion

% Grupo control:
nA_control = 1; % numero de histogramas procesados (modificar)
A_control = readmatrix('E_control.txt');
A_control = repelem(A_control(:, 1), ceil((A_control(:, 2)/(100*nE_control))*size(A_control, 1))); % modulo de Young en kPa
A_control_mean = mean(A_control); % media
A_control_deviation = std(A_control); % desviacion estandar

% Grupo modificado:
nA_modified = 1; % numero de histogramas procesados (modificar)
A_modified = readmatrix('E_modified.txt');
A_modified = repelem(A_modified(:, 1), ceil((A_modified(:, 2)/(100*nE_modified))*size(A_modified, 1))); % modulo de Young en kPa
A_modified_mean = mean(A_modified); % media
A_modified_deviation = std(A_modified); % desviacion estandar

% Resultados:
disp('--- Resultados para fuerza máxima de adhesión ---');
disp(' ');
disp(['Grupo control: ' num2str(length(A_control)) ' curvas de fuerza analizadas.']);
disp(['Fuerza de Adhesión = ' num2str(A_control_mean,'%.2e') ' +/- ' num2str(A_control_deviation,'%.2e') ' nN']);
disp(' ');
disp(['Grupo modificado: ' num2str(length(A_modified)) ' curvas de fuerza analizadas.']);
disp(['Fuerza de Adhesión = ' num2str(A_modified_mean,'%.2e') ' +/- ' num2str(A_modified_deviation,'%.2e') ' nN']);
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

%% Histograma - Fuerza de Adhesion

figure;
histogram(A_control, 100, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, 0.01]);
hold on;
histogram(A_modified, 100, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, 0.01]);
hold off;
xlabel("Adhesion Force (nN)");
ylabel("Force curves");
legend("Control", "Modified");
title("Adhesion Force");
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

%% Boxplot - Fuerza de Adhesion

% Preparo datos para hacer el boxplot:
A_data = [A_control; A_modified];
A_boxplot = [repmat({'Control'}, length(A_control), 1); repmat({'Modified'}, length(A_modified), 1)];

% Grafico el boxplot:
figure;
boxplot(A_data, A_boxplot, 'Labels', {'Control', 'Modified'});
ylim([0, 0.01]);
title("Adhesion Force (nN)");
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

%% Test de Student - Fuerza de Adhesion

[hA, pA] = ttest2(A_control, A_modified);

disp('--- Resultados del test de Student para fuerza máxima de adhesión ---');
disp(' ');
if hE == 1
  disp('Las muestras son significativamente diferentes.');
end
if hE == 0
  disp('Las muestras NO son significativamente diferentes.');
end
disp(['p-valor = ' num2str(pA)]);
disp(' ');