%% Analisis estadistico de resultados

% Este programa realiza el analisis estadistico de los valores de
% elasticidad y adhesion obtenidos de
% curvas de fuerza de nanoindentacion o mapas QNM.

% Los archivos con los datos deben estar en la misma carpeta que el script.

% Limpio el espacio de trabajo:
close all
clear
clc

%% --- CONFIGURAR PARAMETROS ---

% Indico el metodo por el que obtuve los datos:
metodo = 1; % 1: Nanoindentacion / 2: QNM

% Defino la cantidad de muestras (condiciones, grupos, poblaciones):
muestras = 2; % Minimo 1, maximo 4

% Defino nombres de los grupos (poblaciones, condiciones) para leyendas de graficos:
muestra1 = 'Control';
muestra2 = 'Nombre';
muestra3 = 'Nombre';
muestra4 = 'Nombre';

% Elijo si quiero analizar adhesion o solo elasticidad:
adhesion = 0; % 0: Analiza solo elasticidad / 1: Analiza elasticidad y adhesion

% Elijo si quiero graficar boxplot o solo histograma:
boxpt = 0; % 0: Grafica solo histograma / 1: Grafica histograma y boxplot

% --- Parametros para analisis de elasticidad ---

% Indico los nombres de los archivos con los datos a analizar:
E_muestra1 = 'E_muestra1.txt';
E_muestra2 = 'E_muestra2.txt';
E_muestra3 = 'E_muestra3.txt';
E_muestra4 = 'E_muestra4.txt';

% Solo QNM: indico numero de archivos exportados para cada grupo de datos:
nE1 = 1; % Numero de archivos exportados de elasticidad del grupo 1
nE2 = 1; % Numero de archivos exportados de elasticidad del grupo 2
nE3 = 1; % Numero de archivos exportados de elasticidad del grupo 3
nE4 = 1; % Numero de archivos exportados de elasticidad del grupo 4

% Defino filtro limites de corte de valores de elasticidad:
minE1 = 0; % Valor minimo del grupo 1
maxE1 = 100; % Valor maximo del grupo 1
minE2 = 0; % Valor minimo del grupo 2
maxE2 = 100; % Valor maximo del grupo 2
minE3 = 0; % Valor minimo del grupo 3
maxE3 = 100; % Valor maximo del grupo 3
minE4 = 0; % Valor minimo del grupo 4
maxE4 = 100; % Valor maximo del grupo 4

% Defino limite superior de ejes de histograma y boxplot:
limE = 100; % Elasticidad

% Defino numero de bins de los histogramas de elasticidad:
nbinsE1 = 100; % Numero de bins del grupo 1
nbinsE2 = 100; % Numero de bins del grupo 2
nbinsE3 = 100; % Numero de bins grupo 3
nbinsE4 = 100; % Numero de bins grupo 4

% --- Parametros para analisis de adhesion ---

% Indico los nombres de los archivos con los datos a analizar:
A_muestra1 = 'A_muestra1.txt';
A_muestra2 = 'A_muestra2.txt';
A_muestra3 = 'A_muestra3.txt';
A_muestra4 = 'A_muestra4.txt';

% Solo QNM: indico numero de archivos exportados para cada grupo de datos:
nA1 = 1; % Numero de archivos exportados de adhesion del grupo 1
nA2 = 1; % Numero de archivos exportados de adhesion del grupo 2
nA3 = 1; % Numero de archivos exportados de adhesion del grupo 3
nA4 = 1; % Numero de archivos exportados de adhesion del grupo 4

% Defino filtro limites de corte de valores de adhesion:
minA1 = 0; % Valor minimo del grupo 1
maxA1 = 100; % Valor maximo del grupo 1
minA2 = 0; % Valor minimo del grupo 2
maxA2 = 100; % Valor maximo del grupo 2
minA3 = 0; % Valor minimo del grupo 3
maxA3 = 100; % Valor maximo del grupo 3
minA4 = 0; % Valor minimo del grupo 4
maxA4 = 100; % Valor maximo del grupo 4

% Defino limite superior de ejes de histograma y boxplot:
limA = 0.01; % Adhesion

% Defino numero de bins de los histogramas de adhesion:
nbinsA1 = 100; % Numero de bins del grupo 1
nbinsA2 = 100; % Numero de bins del grupo 2
nbinsA3 = 100; % Numero de bins del grupo 3
nbinsA4 = 100; % Numero de bins del grupo 4

%% --- ELASTICIDAD ---

% --- Leer, analizar y mostrar resultados ---

disp('--- Resultados para módulo de Young ---');
disp(' ');

% Grupo 1:
E_muestra1 = readmatrix(E_muestra1); % Leo archivo
% Extraigo vector de modulo de Young en kPa:
if metodo == 1
  E_muestra1 = E_muestra1(:,4); % Nanoindentacion
else
  E_muestra1 = repelem(E_muestra1(:, 1), ceil((E_muestra1(:, 2)/(100*nE1))*size(E_muestra1, 1))); % QNM
end
E_muestra1 = E_muestra1(E_muestra1>minE1 & E_muestra1<maxE1); % Filtro datos con limites de corte
E_muestra1_mean = mean(E_muestra1); % Media
E_muestra1_deviation = std(E_muestra1); % Desviacion estandar
[hEks1, pEks1] = kstest(E_muestra1, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
% Resultados:
disp([muestra1 ': ' num2str(length(E_muestra1)) ' curvas de fuerza analizadas.']);
disp(['Módulo de Young = ' num2str(E_muestra1_mean,'%.2f') ' +/- ' num2str(E_muestra1_deviation,'%.2f') ' kPa']);
if hEks1 == 0
  disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
end
if hEks1 == 1
  disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
end
disp(' ');

% Grupo 2:
if muestras >= 2
  E_muestra2 = readmatrix(E_muestra2); % Leo archivo
  % Extraigo vector de modulo de Young en kPa:
  if metodo == 1
    E_muestra2 = E_muestra2(:,4); % Nanoindentacion
  else
    E_muestra2 = repelem(E_muestra2(:, 1), ceil((E_muestra2(:, 2)/(100*nE2))*size(E_muestra2, 1))); % QNM
  end
  E_muestra2 = E_muestra2(E_muestra2>minE2 & E_muestra2<maxE2); % Filtro datos con limites de corte
  E_muestra2_mean = mean(E_muestra2); % Media
  E_muestra2_deviation = std(E_muestra2); % Desviacion estandar
  [hEks2, pEks2] = kstest(E_muestra2, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
  % Resultados:
  disp([muestra2 ': ' num2str(length(E_muestra2)) ' curvas de fuerza analizadas.']);
  disp(['Módulo de Young = ' num2str(E_muestra2_mean,'%.2f') ' +/- ' num2str(E_muestra2_deviation,'%.2f') ' kPa']);
  if hEks2 == 0
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
  end
  if hEks2 == 1
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
  end
  disp(' ');
end

% Grupo 3:
if muestras >= 3
  E_muestra3 = readmatrix(E_muestra3); % Leo archivo
  % Extraigo vector de modulo de Young en kPa:
  if metodo == 1
    E_muestra3 = E_muestra3(:,4); % Nanoindentacion
  else
    E_muestra3 = repelem(E_muestra3(:, 1), ceil((E_muestra3(:, 2)/(100*nE3))*size(E_muestra3, 1))); % QNM
  end
  E_muestra3 = E_muestra3(E_muestra3>minE3 & E_muestra3<maxE3); % Filtro datos con limites de corte
  E_muestra3_mean = mean(E_muestra3); % Media
  E_muestra3_deviation = std(E_muestra3); % Desviacion estandar
  [hEks3, pEks3] = kstest(E_muestra3, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
  % Resultados:
  disp([muestra3 ': ' num2str(length(E_muestra3)) ' curvas de fuerza analizadas.']);
  disp(['Módulo de Young = ' num2str(E_muestra3_mean,'%.2f') ' +/- ' num2str(E_muestra3_deviation,'%.2f') ' kPa']);
  if hEks3 == 0
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
  end
  if hEks3 == 1
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
  end
  disp(' ');
end

% Grupo 4:
if muestras >= 4
  E_muestra4 = readmatrix(E_muestra4); % Leo archivo
  % Extraigo vector de modulo de Young en kPa:
  if metodo == 1
    E_muestra4 = E_muestra4(:,4); % Nanoindentacion
  else
    E_muestra4 = repelem(E_muestra4(:, 1), ceil((E_muestra4(:, 2)/(100*nE4))*size(E_muestra4, 1))); % QNM
  end
  E_muestra4 = E_muestra4(E_muestra4>minE4 & E_muestra4<maxE4); % Filtro datos con limites de corte
  E_muestra4_mean = mean(E_muestra4); % Media
  E_muestra4_deviation = std(E_muestra4); % Desviacion estandar
  [hEks4, pEks4] = kstest(E_muestra4, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
  % Resultados:
  disp([muestra4 ': ' num2str(length(E_muestra4)) ' curvas de fuerza analizadas.']);
  disp(['Módulo de Young = ' num2str(E_muestra4_mean,'%.2f') ' +/- ' num2str(E_muestra4_deviation,'%.2f') ' kPa']);
  if hEks4 == 0
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
  end
  if hEks4 == 1
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
  end
  disp(' ');
end

% --- Prueba T de Student (para 2 muestras que siguen una distribucion normal) ---

if muestras >= 2 && hEks1 == 0 && hEks2 == 0
  [hE, pE] = ttest2(E_muestra1, E_muestra2); % Prueba T
  % Resultados:
  disp('--- Resultados de la prueba T de Student para módulo de Young ---');
  disp(' ');
  if hE == 1
    disp(['Las muestras ' muestra1 'y' muestra2 'son significativamente diferentes.']);
  end
  if hE == 0
    disp(['Las muestras ' muestra1 'y' muestra2 'NO son significativamente diferentes.']);
  end
  disp(['p-valor = ' num2str(pE)]);
  disp(' ');
end

% --- Prueba U de Mann-Whitney (para 2 muestras que NO siguen una distribucion normal) ---

if muestras >= 2 && hEks1 == 1 || hEks2 == 1
  [pE, hE] = ranksum(E_muestra1, E_muestra2); % Prueba U
  % Resultados:
  disp('--- Resultados de la prueba U de Mann-Whitney para módulo de Young ---');
  disp(' ');
  if hE == 1
    disp(['Las muestras ' muestra1 'y' muestra2 'son significativamente diferentes.']);
  end
  if hE == 0
    disp(['Las muestras ' muestra1 'y' muestra2 'NO son significativamente diferentes.']);
  end
  disp(['p-valor = ' num2str(pE)]);
  disp(' ');
end

% --- Histograma ---

x = linspace(0, limE);
figure;
title("Young's Modulus");
xlabel("Young's Modulus (kPa)");
ylabel("Force curves");
grid on;

% Para muestra única:
hE_muestra1 = histogram(E_muestra1, nbinsE1, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limE]);
if muestras == 1
  legend(muestra1); % Leyenda para muestra única
end
hold on;
% Superpongo el ajuste gaussiano a la grafica del histograma (si la muestra sigue una distribucion normal):
if muestras <= 2 && hEks1 == 0 % Solo para 1 o 2 muestras
  fit_muestra1 = normpdf(x, E_muestra1_mean, E_muestra1_deviation); % Ajuste gaussiano para grupo 1
  plot(x, fit_muestra1 * numel(E_muestra1) * diff(hE_muestra1.BinEdges(1:2)), 'b', 'LineWidth', 2);
end

% Para 2 muestras:
if muestras >= 2
  hE_muestra2 = histogram(E_muestra2, nbinsE2, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limE]);
  if muestras == 2
    legend(muestra1, muestra2); % Leyenda para 2 muestras
  end
  % Superpongo el ajuste gaussiano a la grafica del histograma (si la muestra sigue una distribucion normal):
  if muestras <= 2 && hEks1 == 0 && hEks2 == 0 % Solo para 1 o 2 muestras
    fit_muestra2 = normpdf(x, E_muestra2_mean, E_muestra2_deviation); % Ajuste gaussiano para grupo 2
    plot(x, fit_muestra2 * numel(E_muestra2) * diff(hE_muestra2.BinEdges(1:2)), 'r', 'LineWidth', 2);
  end
  
  % Para 3 y 4 muestras:
  if muestras >= 3
    hE_muestra3 = histogram(E_muestra3, nbinsE3, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limE]);
    if muestras == 3
      legend(muestra1, muestra2, muestra3); % Leyenda para 3 muestras
    end
    if muestras >= 4
      hE_muestra4 = histogram(E_muestra4, nbinsE4, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limE]);
      if muestras == 4
        legend(muestra1, muestra2, muestra3, muestra4); % Leyenda para 4 muestras
      end
    end
  end
end
hold off;

% --- Boxplot ---

% Para 2 muestras:
if boxpt == 1 && muestras == 2
  % Preparo datos para hacer el boxplot:
  E_data = [E_muestra1; E_muestra2];
  E_boxplot = [repmat({'muestra1'}, length(E_muestra1), 1);
    repmat({'muestra2'}, length(E_muestra2), 1)];
  % Grafico el boxplot:
  figure;
  boxplot(E_data, E_boxplot, 'Labels', {muestra1, muestra2}); % Etiquetas
  ylim([0, limE]); % Limites
  title("Young's Modulus (kPa)"); % Titulo
  grid on;
end

% Para 3 muestras:
if boxpt == 1 && muestras == 3
  % Preparo datos para hacer el boxplot:
  E_data = [E_muestra1; E_muestra2; E_muestra3];
  E_boxplot = [repmat({'muestra1'}, length(E_muestra1), 1);
    repmat({'muestra2'}, length(E_muestra2), 1);
    repmat({'muestra3'}, length(E_muestra3), 1)];
  % Grafico el boxplot:
  figure;
  boxplot(E_data, E_boxplot, 'Labels', {muestra1, muestra2, muestra3}); % Etiquetas
  ylim([0, limE]); % Limites
  title("Young's Modulus (kPa)"); % Titulo
  grid on;
end

% Para 4 muestras:
if boxpt == 1 && muestras == 4
  % Preparo datos para hacer el boxplot:
  E_data = [E_muestra1; E_muestra2; E_muestra3; E_muestra4];
  E_boxplot = [repmat({'muestra1'}, length(E_muestra1), 1);
    repmat({'muestra2'}, length(E_muestra2), 1);
    repmat({'muestra3'}, length(E_muestra3), 1);
    repmat({'muestra4'}, length(E_muestra4), 1)];
  % Grafico el boxplot:
  figure;
  boxplot(E_data, E_boxplot, 'Labels', {muestra1, muestra2, muestra3, muestra4}); % Etiquetas
  ylim([0, limE]); % Limites
  title("Young's Modulus (kPa)"); % Titulo
  grid on;
end

%% --- ADHESION ---

if adhesion == 1
  
  % --- Leer, analizar y mostrar resultados ---
  
  if metodo == 1
    disp('--- Resultados para trabajo de adhesión ---');
  else
    disp('--- Resultados para fuerza máxima de adhesión ---');
  end
  disp(' ');
  
  % Grupo 1:
  A_muestra1 = readmatrix(A_muestra1); % Leo archivo
  % Extraigo vector de valores de adhesion:
  if metodo == 1
    A_muestra1 = A_muestra1(:,4); % Nanoindentacion: trabajo de adhesion en J/m2
  else
    A_muestra1 = repelem(A_muestra1(:, 1), ceil((A_muestra1(:, 2)/(100*nA1))*size(A_muestra1, 1))); % QNM: fuerza maxima de adhesion en nN
  end
  A_muestra1 = A_muestra1(A_muestra1>minA1 & A_muestra1<maxA1); % Filtro datos con limites de corte
  A_muestra1_mean = mean(A_muestra1); % Media
  A_muestra1_deviation = std(A_muestra1); % Desviacion estandar
  [hAks1, pAks1] = kstest(A_muestra1, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
  % Resultados:
  disp([muestra1 ': ' num2str(length(A_muestra1)) ' curvas de fuerza analizadas.']);
  if metodo == 1
    disp(['Trabajo de adhesión = ' num2str(A_muestra1_mean,'%.2f') ' +/- ' num2str(A_muestra1_deviation,'%.2f') ' J/m2']);
  else
    disp(['Fuerza máxima de adhesión = ' num2str(A_muestra1_mean,'%.2f') ' +/- ' num2str(A_muestra1_deviation,'%.2f') ' nN']);
  end
  if hAks1 == 0
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
  end
  if hAks1 == 1
    disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
  end
  disp(' ');
  
  % Grupo 2:
  if muestras >= 2
    A_muestra2 = readmatrix(A_muestra2); % Leo archivo
    % Extraigo vector de valores de adhesion:
    if metodo == 1
      A_muestra2 = A_muestra2(:,4); % Nanoindentacion: trabajo de adhesion en J/m2
    else
      A_muestra2 = repelem(A_muestra2(:, 1), ceil((A_muestra2(:, 2)/(100*nA2))*size(A_muestra2, 1))); % QNM: fuerza maxima de adhesion en nN
    end
    A_muestra2 = A_muestra2(A_muestra2>minA2 & A_muestra2<maxA2); % Filtro datos con limites de corte
    A_muestra2_mean = mean(A_muestra2); % Media
    A_muestra2_deviation = std(A_muestra2); % Desviacion estandar
    [hAks2, pAks2] = kstest(A_muestra2, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
    % Resultados:
    disp([muestra2 ': ' num2str(length(A_muestra2)) ' curvas de fuerza analizadas.']);
    if metodo == 1
      disp(['Trabajo de adhesión = ' num2str(A_muestra2_mean,'%.2f') ' +/- ' num2str(A_muestra2_deviation,'%.2f') ' J/m2']);
    else
      disp(['Fuerza máxima de adhesión = ' num2str(A_muestra2_mean,'%.2f') ' +/- ' num2str(A_muestra2_deviation,'%.2f') ' nN']);
    end
    if hAks2 == 0
      disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
    end
    if hAks2 == 1
      disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
    end
    disp(' ');
  end
  
  % Grupo 3:
  if muestras >= 3
    A_muestra3 = readmatrix(A_muestra3); % Leo archivo
    % Extraigo vector de valores de adhesion:
    if metodo == 1
      A_muestra3 = A_muestra3(:,4); % Nanoindentacion: trabajo de adhesion en J/m2
    else
      A_muestra3 = repelem(A_muestra3(:, 1), ceil((A_muestra3(:, 2)/(100*nA3))*size(A_muestra3, 1))); % QNM: fuerza maxima de adhesion en nN
    end
    A_muestra3 = A_muestra3(A_muestra3>minA3 & A_muestra3<maxA3); % Filtro datos con limites de corte
    A_muestra3_mean = mean(A_muestra3); % Media
    A_muestra3_deviation = std(A_muestra3); % Desviacion estandar
    [hAks3, pAks3] = kstest(A_muestra3, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
    % Resultados:
    disp([muestra3 ': ' num2str(length(A_muestra3)) ' curvas de fuerza analizadas.']);
    if metodo == 1
      disp(['Trabajo de adhesión = ' num2str(A_muestra3_mean,'%.2f') ' +/- ' num2str(A_muestra3_deviation,'%.2f') ' J/m2']);
    else
      disp(['Fuerza máxima de adhesión = ' num2str(A_muestra3_mean,'%.2f') ' +/- ' num2str(A_muestra3_deviation,'%.2f') ' nN']);
    end
    if hAks3 == 0
      disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
    end
    if hAks3 == 1
      disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
    end
    disp(' ');
  end
  
  % Grupo 4:
  if muestras >= 4
    A_muestra4 = readmatrix(A_muestra4); % Leo archivo
    % Extraigo vector de valores de adhesion:
    if metodo == 1
      A_muestra4 = A_muestra4(:,4); % Nanoindentacion: trabajo de adhesion en J/m2
    else
      A_muestra4 = repelem(A_muestra4(:, 1), ceil((A_muestra4(:, 2)/(100*nA4))*size(A_muestra4, 1))); % QNM: fuerza maxima de adhesion en nN
    end
    A_muestra4 = A_muestra4(A_muestra4>minA4 & A_muestra4<maxA4); % Filtro datos con limites de corte
    A_muestra4_mean = mean(A_muestra4); % Media
    A_muestra4_deviation = std(A_muestra4); % Desviacion estandar
    [hAks4, pAks4] = kstest(A_muestra4, 'Alpha', 0.05); % Prueba de normalidad de Kolmogorov-Smirnov
    % Resultados:
    disp([muestra4 ': ' num2str(length(A_muestra4)) ' curvas de fuerza analizadas.']);
    if metodo == 1
      disp(['Trabajo de adhesión = ' num2str(A_muestra4_mean,'%.2f') ' +/- ' num2str(A_muestra4_deviation,'%.2f') ' J/m2']);
    else
      disp(['Fuerza máxima de adhesión = ' num2str(A_muestra4_mean,'%.2f') ' +/- ' num2str(A_muestra4_deviation,'%.2f') ' nN']);
    end
    if hAks4 == 0
      disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra sigue una distribución normal.');
    end
    if hAks4 == 1
      disp('Prueba de normalidad de Kolmogorov-Smirnov: La muestra NO sigue una distribución normal.');
    end
    disp(' ');
  end
  
  % --- Prueba T de Student (para 2 muestras que siguen una distribucion normal) ---
  
  if muestras >= 2 && hAks1 == 0 && hAks2 == 0
    [hA, pA] = ttest2(A_muestra1, A_muestra2); % Prueba T
    % Resultados:
    if metodo == 1
      disp('--- Resultados de la prueba T de Student para trabajo de adhesión ---');
    else
      disp('--- Resultados de la prueba T de Student para fuerza máxima de adhesión ---');
    end
    disp(' ');
    if hA == 1
      disp(['Las muestras ' muestra1 'y' muestra2 'son significativamente diferentes.']);
    end
    if hA == 0
      disp(['Las muestras ' muestra1 'y' muestra2 'NO son significativamente diferentes.']);
    end
    disp(['p-valor = ' num2str(pA)]);
    disp(' ');
  end
  
  % --- Prueba U de Mann-Whitney (para 2 muestras que NO siguen una distribucion normal) ---
  
  if muestras >= 2 && hAks1 == 1 || hAks2 == 1
    [pA, hA] = ranksum(A_muestra1, A_muestra2); % Prueba U
    % Resultados:
    if metodo == 1
      disp('--- Resultados de la prueba U de Mann-Whitney para trabajo de adhesión ---');
    else
      disp('--- Resultados de la prueba U de Mann-Whitney para fuerza máxima de adhesión ---');
    end
    disp(' ');
    if hA == 1
      disp(['Las muestras ' muestra1 'y' muestra2 'son significativamente diferentes.']);
    end
    if hA == 0
      disp(['Las muestras ' muestra1 'y' muestra2 'NO son significativamente diferentes.']);
    end
    disp(['p-valor = ' num2str(pA)]);
    disp(' ');
  end
  
  % --- Histograma ---
  
  x = linspace(0, limA);
  figure;
  if metodo == 1
    title("Work of Adhesion"); % Nanoindentacion
    xlabel("Work of Adhesion (J/m2)");
  else
    title("Adhesion Force"); % QNM
    xlabel("Adhesion Force (nN)");
  end
  ylabel("Force curves");
  grid on;
  
  % Para muestra única:
  hA_muestra1 = histogram(A_muestra1, nbinsA1, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limA]);
  if muestras == 1
    legend(muestra1); % Leyenda para muestra única
  end
  hold on;
  % Superpongo el ajuste gaussiano a la grafica del histograma (si la muestra sigue una distribucion normal):
  if muestras <= 2 && hAks1 == 0 % Solo para 1 o 2 muestras
    fit_muestra1 = normpdf(x, A_muestra1_mean, A_muestra1_deviation); % Ajuste gaussiano para grupo 1
    plot(x, fit_muestra1 * numel(A_muestra1) * diff(hA_muestra1.BinAdges(1:2)), 'b', 'LineWidth', 2);
  end
  
  % Para 2 muestras:
  if muestras >= 2
    hA_muestra2 = histogram(A_muestra2, nbinsA2, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limA]);
    if muestras == 2
      legend(muestra1, muestra2); % Leyenda para 2 muestras
    end
    % Superpongo el ajuste gaussiano a la grafica del histograma (si la muestra sigue una distribucion normal):
    if muestras <= 2 && hAks1 == 0 && hAks2 == 0 % Solo para 1 o 2 muestras
      fit_muestra2 = normpdf(x, A_muestra2_mean, A_muestra2_deviation); % Ajuste gaussiano para grupo 2
      plot(x, fit_muestra2 * numel(A_muestra2) * diff(hA_muestra2.BinAdges(1:2)), 'r', 'LineWidth', 2);
    end
    
    % Para 3 y 4 muestras:
    if muestras >= 3
      hA_muestra3 = histogram(A_muestra3, nbinsA3, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limA]);
      if muestras == 3
        legend(muestra1, muestra2, muestra3); % Leyenda para 3 muestras
      end
      if muestras >= 4
        hA_muestra4 = histogram(A_muestra4, nbinsA4, 'Normalization', 'count', 'EdgeColor', 'none', 'BinLimits', [0, limA]);
        if muestras == 4
          legend(muestra1, muestra2, muestra3, muestra4); % Leyenda para 4 muestras
        end
      end
    end
  end
  hold off;
  
  % --- Boxplot ---
  
  % Para 2 muestras:
  if boxpt == 1 && muestras == 2
    % Preparo datos para hacer el boxplot:
    A_data = [A_muestra1; A_muestra2];
    A_boxplot = [repmat({'muestra1'}, length(A_muestra1), 1);
      repmat({'muestra2'}, length(A_muestra2), 1)];
    % Grafico el boxplot:
    figure;
    boxplot(A_data, A_boxplot, 'Labels', {muestra1, muestra2}); % Etiquetas
    ylim([0, limA]); % Limites
    if metodo == 1
      title("Work of Adhesion (J/m2)"); % Nanoindentacion
    else
      title("Adhesion Force (nN)"); % QNM
    end
    grid on;
  end
  
  % Para 3 muestras:
  if boxpt == 1 && muestras == 3
    % Preparo datos para hacer el boxplot:
    A_data = [A_muestra1; A_muestra2; A_muestra3];
    A_boxplot = [repmat({'muestra1'}, length(A_muestra1), 1);
      repmat({'muestra2'}, length(A_muestra2), 1);
      repmat({'muestra3'}, length(A_muestra3), 1)];
    % Grafico el boxplot:
    figure;
    boxplot(A_data, A_boxplot, 'Labels', {muestra1, muestra2, muestra3}); % Etiquetas
    ylim([0, limA]); % Limites
    if metodo == 1
      title("Work of Adhesion (J/m2)"); % Nanoindentacion
    else
      title("Adhesion Force (nN)"); % QNM
    end
    grid on;
  end
  
  % Para 4 muestras:
  if boxpt == 1 && muestras == 4
    % Preparo datos para hacer el boxplot:
    A_data = [A_muestra1; A_muestra2; A_muestra3; A_muestra4];
    A_boxplot = [repmat({'muestra1'}, length(A_muestra1), 1);
      repmat({'muestra2'}, length(A_muestra2), 1);
      repmat({'muestra3'}, length(A_muestra3), 1);
      repmat({'muestra4'}, length(A_muestra4), 1)];
    % Grafico el boxplot:
    figure;
    boxplot(A_data, A_boxplot, 'Labels', {muestra1, muestra2, muestra3, muestra4}); % Etiquetas
    ylim([0, limA]); % Limites
    if metodo == 1
      title("Work of Adhesion (J/m2)"); % Nanoindentacion
    else
      title("Adhesion Force (nN)"); % QNM
    end
    grid on;
  end
  
end