%% Concatenar archivos

% Limpio el espacio de trabajo:
close all
clear
clc

%% --- CONFIGURAR PARAMETROS ---

% Defino directorio donde estan los archivos a concatenar:
directorio = [uigetdir('C:\ruta_archivos') '\']; % Modificar ruta
archivos = dir([directorio 'Archivo0*']); % Modificar nombre generico de los archivos

% Defino nombre del archivo concatenado a guardar:
nombre = 'E_control.txt'; % Modificar nombre del archivo

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