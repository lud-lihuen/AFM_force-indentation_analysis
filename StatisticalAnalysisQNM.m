%% Analisis estadistico de resultados de mapas QNM

close all
clear
clc

%% --- CONFIGURAR PARAMETROS ---

% Defino directorio donde estan los archivos a concatenar:
directorio = [uigetdir('C:\ruta_archivos') '\']; % Modificar ruta
archivos = dir([directorio 'Archivo0*']); % Modificar nombre generico de los archivos

% Defino nombre del archivo concatenado a guardar:
nombre = 'E_control.txt'; % Modificar nombre del archivo

% Defino numero de archivos exportados para cada grupo de datos:
nEc = 1; % Numero de archivos exportados de modulo de Young del grupo control
nEm = 1; % Numero de archivos exportados de modulo de Young del grupo modificado
nAc = 1; % Numero de archivos exportados de fuerza maxima de adhesion del grupo control
nAm = 1; % Numero de archivos exportados de fuerza maxima de adhesion del grupo modificado

% Defino filtro limite superior de valores de modulo de Young y fuerza maxima de adhesion:
Eclim = 100; % Limite de corte de modulo de Young del grupo control
Emlim = 100; % Limite de corte de modulo de Young del grupo modificado
Aclim = 100; % Limite de corte de fuerza maxima de adhesion del grupo control
Amlim = 100; % Limite de corte de fuerza maxima de adhesion del grupo modificado

% Defino limite superior de ejes de histograma y boxplot:
Elim = 100; % Modulo de Young
Alim = 0.01; % Fuerza maxima de adhesion

% Defino numero de bins de los histogramas:
nbinsEc = 100; % Numero de bins de modulo de Young del grupo control
nbinsEm = 100; % Numero de bins de modulo de Young del grupo modificado
nbinsAc = 100; % Numero de bins de fuerza maxima de adhesion del grupo control
nbinsAm = 100; % Numero de bins de fuerza maxima de adhesion del grupo modificado

% Defino nombres de los grupos para leyendas de graficos:
control = 'Control';
modified = 'Modified';

%% --- CONCATENAR ARCHIVOS ---

% Concateno los archivos exportados de NanoScope Analysis:
datos = [];
for i = 1:length(archivos)
    archivo = [directorio getfield(archivos(i),'name')]; % Obtengo la ruta de cada archivo
    datos_i = readmatrix(archivo); % Leo cada archivo
    datos = cat(1, datos, datos_i); % Concateno en matriz de resultados
end

% Guardo la matriz de resultados en formato ascii:
save ([directorio nombre],'datos','-ascii');

%% --- ELASTICIDAD ---

% --- Leer resultados ---

% Primero copiar archivos E_control.txt y E_modified.txt a carpeta del script

% Grupo control:
E_control = readmatrix('E_control.txt'); % Leo archivo
E_control = repelem(E_control(:, 1), ceil((E_control(:, 2)/(100*nEc))*size(E_control, 1))); % Extraigo vector de modulo de Young en kPa
E_control = E_control(E_control>0&E_control<Eclim); % Filtro datos con limites de corte
E_control_mean = mean(E_control); % Media
E_control_deviation = std(E_control); % Desviacion estandar

% Grupo modificado:
E_modified = readmatrix('E_modified.txt'); % Leo archivo
E_modified = repelem(E_modified(:, 1), ceil((E_modified(:, 2)/(100*nEm))*size(E_modified, 1))); % Extraigo vector de modulo de Young en kPa
E_modified = E_modified(E_modified>0&E_modified<Emlim); % Filtro datos con limites de corte
E_modified_mean = mean(E_modified); % Media
E_modified_deviation = std(E_modified); % Desviacion estandar

% Resultados:
disp('--- Resultados para módulo de Young ---');
disp(' ');
disp(['Grupo control: ' num2str(length(E_control)) ' curvas de fuerza analizadas.']);
disp(['Módulo de Young = ' num2str(E_control_mean,'%.2f') ' +/- ' num2str(E_control_deviation,'%.2f') ' kPa']);
disp(' ');
disp(['Grupo modificado: ' num2str(length(E_modified)) ' curvas de fuerza analizadas.']);
disp(['Módulo de Young = ' num2str(E_modified_mean,'%.2f') ' +/- ' num2str(E_modified_deviation,'%.2f') ' kPa']);
disp(' ');

% --- Histograma ---

% Preparo el ajuste gaussiano para superponer al histograma:
x = linspace(0, Elim);
fit_control = normpdf(x, E_control_mean, E_control_deviation);
fit_modified = normpdf(x, E_modified_mean, E_modified_deviation);

% Realizo el histograma:
figure;
hE_control = histogram(E_control, nbinsEc, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, Elim]);
hold on;
hE_modified = histogram(E_modified, nbinsEm, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, Elim]);

% Superpongo el ajuste gaussiano a la grafica del histograma:
plot(x, fit_control * numel(E_control) * diff(hE_control.BinEdges(1:2)), 'b', 'LineWidth', 2);
plot(x, fit_modified * numel(E_modified) * diff(hE_modified.BinEdges(1:2)), 'r', 'LineWidth', 2);
hold off;

% Etiquetas y leyendas:
title("Young's Modulus");
xlabel("Young's Modulus (kPa)");
ylabel("Force curves");
legend(control, modified);
grid on;

% --- Boxplot ---

% Preparo datos para hacer el boxplot:
E_data = [E_control; E_modified];
E_boxplot = [repmat({'Control'}, length(E_control), 1); repmat({'Modified'}, length(E_modified), 1)];

% Grafico el boxplot:
figure;
boxplot(E_data, E_boxplot, 'Labels', {control, modified}); % Etiquetas
ylim([0, Elim]); % Limites
title("Young's Modulus (kPa)"); % Titulo
grid on;

% --- Prueba de normalidad de Kolmogorov-Smirnov ---

% Grupo control:
[hEksc, pEksc] = kstest(E_control, 'Alpha', 0.05);

% Grupo modificado:
[hEksm, pEksm] = kstest(E_modified, 'Alpha', 0.05);

% Resultados:
disp('--- Resultados de la prueba de normalidad de Kolmogorov-Smirnov para módulo de Young ---');
disp(' ');
if hEksc == 0
  disp('Grupo control: La muestra sigue una distribución normal.');
end
if hEksc == 1
  disp('Grupo control: La muestra NO sigue una distribución normal.');
end
disp(' ');
if hEksm == 0
  disp('Grupo modificado: La muestra sigue una distribución normal.');
end
if hEksm == 1
  disp('Grupo modificado: La muestra NO sigue una distribución normal.');
end
disp(' ');

% --- Prueba T de Student (si la muestra sigue una distribucion normal) ---

if hEksc == 0 && hEksm == 0
    % Prueba T:
    [hE, pE] = ttest2(E_control, E_modified);
    
    % Resultados:
    disp('--- Resultados de la prueba T de Student para módulo de Young ---');
    disp(' ');
    if hE == 1
      disp('Las muestras son significativamente diferentes.');
    end
    if hE == 0
      disp('Las muestras NO son significativamente diferentes.');
    end
    disp(['p-valor = ' num2str(pE)]);
    disp(' ');
end

% --- Prueba U de Mann-Whitney (si la muestra NO sigue una distribucion normal) ---

if hEksc == 1 || hEksm == 1
    % Prueba U:
    [pE, hE] = ranksum(E_control, E_modified);

    % Resultados:
    disp('--- Resultados de la prueba U de Mann-Whitney para módulo de Young ---');
    disp(' ');
    if hE == 1
      disp('Las muestras son significativamente diferentes.');
    end
    if hE == 0
      disp('Las muestras NO son significativamente diferentes.');
    end
    disp(['p-valor = ' num2str(pE)]);
    disp(' ');
end
    
%% --- ADHESION ---

% --- Leer resultados ---

% Primero copiar archivos A_control.txt y A_modified.txt a carpeta del script

% Grupo control:
A_control = readmatrix('A_control.txt'); % Leo archivo
A_control = repelem(A_control(:, 1), ceil((A_control(:, 2)/(100*nAc))*size(A_control, 1))); % Extraigo vector de fuerza maxima de adhesion en nN
A_control = A_control(A_control>0&A_control<Aclim); % Filtro datos con limites de corte
A_control_mean = mean(A_control); % Media
A_control_deviation = std(A_control); % Desviacion estandar

% Grupo modificado:
A_modified = readmatrix('A_modified.txt'); % Leo archivo
A_modified = repelem(A_modified(:, 1), ceil((A_modified(:, 2)/(100*nAm))*size(A_modified, 1))); % Extraigo vector de fuerza maxima de adhesion en nN
A_modified = A_modified(A_modified>0&A_modified<Amlim); % Filtro datos con limites de corte
A_modified_mean = mean(A_modified); % Media
A_modified_deviation = std(A_modified); % Desviacion estandar

% Resultados:
disp('--- Resultados para fuerza máxima de adhesión ---');
disp(' ');
disp(['Grupo control: ' num2str(length(A_control)) ' curvas de fuerza analizadas.']);
disp(['Fuerza de Adhesión = ' num2str(A_control_mean,'%.2e') ' +/- ' num2str(A_control_deviation,'%.2e') ' nN']);
disp(' ');
disp(['Grupo modificado: ' num2str(length(A_modified)) ' curvas de fuerza analizadas.']);
disp(['Fuerza de Adhesión = ' num2str(A_modified_mean,'%.2e') ' +/- ' num2str(A_modified_deviation,'%.2e') ' nN']);
disp(' ');

% --- Histograma ---

% Preparo el ajuste gaussiano para superponer al histograma:
x = linspace(0, Alim);
fit_control = normpdf(x, A_control_mean, A_control_deviation);
fit_modified = normpdf(x, A_modified_mean, A_modified_deviation);

% Realizo el histograma:
figure;
hA_control = histogram(A_control, nbinsAc, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, Alim]);
hold on;
hA_modified = histogram(A_modified, nbinsAm, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, Alim]);

% Superpongo el ajuste gaussiano a la grafica del histograma:
plot(x, fit_control * numel(A_control) * diff(hA_control.BinEdges(1:2)), 'b', 'LineWidth', 2);
plot(x, fit_modified * numel(A_modified) * diff(hA_modified.BinEdges(1:2)), 'r', 'LineWidth', 2);
hold off;

% Etiquetas y leyendas:
title("Adhesion Force");
xlabel("Adhesion Force (nN)");
ylabel("Force curves");
legend(control, modified);
grid on;

% --- Boxplot ---

% Preparo datos para hacer el boxplot:
A_data = [A_control; A_modified];
A_boxplot = [repmat({'Control'}, length(A_control), 1); repmat({'Modified'}, length(A_modified), 1)];

% Grafico el boxplot:
figure;
boxplot(A_data, A_boxplot, 'Labels', {control, modified}); % Etiquetas
ylim([0, Alim]); % Limites
title("Adhesion Force (nN)"); % Titulo
grid on;

% --- Prueba de normalidad de Kolmogorov-Smirnov ---

% Grupo control:
[hAksc, pAksc] = kstest(A_control, 'Alpha', 0.05);

% Grupo modificado:
[hAksm, pAksm] = kstest(A_modified, 'Alpha', 0.05);

% Resultados:
disp('--- Resultados de la prueba de normalidad de Kolmogorov-Smirnov para fuerza máxima de adhesión ---');
disp(' ');
if hAksc == 0
  disp('Grupo control: La muestra sigue una distribución normal.');
end
if hAksc == 1
  disp('Grupo control: La muestra NO sigue una distribución normal.');
end
disp(' ');
if hAksm == 0
  disp('Grupo modificado: La muestra sigue una distribución normal.');
end
if hAksm == 1
  disp('Grupo modificado: La muestra NO sigue una distribución normal.');
end
disp(' ');

% --- Prueba T de Student (si la muestra sigue una distribucion normal) ---

if hAksc == 0 && hAksm == 0
    % Prueba T:
    [hA, pA] = ttest2(A_control, A_modified);

    % Resultados:
    disp('--- Resultados de la prueba T de Student para fuerza máxima de adhesión ---');
    disp(' ');
    if hA == 1
      disp('Las muestras son significativamente diferentes.');
    end
    if hA == 0
      disp('Las muestras NO son significativamente diferentes.');
    end
    disp(['p-valor = ' num2str(pA)]);
    disp(' ');
end

% --- Prueba U de Mann-Whitney (si la muestra NO sigue una distribucion normal) ---

if hAksc == 1 || hAksm == 1
    % Prueba U:
    [pA, hA] = ranksum(A_control, A_modified);

    % Resultados:
    disp('--- Resultados de la prueba U de Mann-Whitney para fuerza máxima de adhesión ---');
    disp(' ');
    if hA == 1
      disp('Las muestras son significativamente diferentes.');
    end
    if hA == 0
      disp('Las muestras NO son significativamente diferentes.');
    end
    disp(['p-valor = ' num2str(pA)]);
    disp(' ');
end