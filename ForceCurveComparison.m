%% Comparacion de curvas de fuerza

% Este programa grafica dos curvas de fuerza en un mismo par de ejes.
% Las curvas deben estar guardadas en un mismo directorio.

% Limpio el espacio de trabajo:
close all
clear
clc
warning('off', 'all'); % Inhabilito los warnings

%% --- SETTINGS ---

% Defino directorio donde estan los archivos con los datos a analizar:
directorio = [uigetdir('C:\ruta_archivos') '\']; % Modificar ruta
archivos = dir([directorio 'Archivo0*']); % Modificar nombre generico de los archivos, dejar el primer cero

% Defino umbral de tolerancia del coeficiente de correlacion para el modulo de Young:
threshold = 0.95; % Modificar entre 0.95 y 0.99 segun de mejor resultado
check = 0; % Elijo si quiero visualizar el ajuste que hizo (0: No / 1: Si)

% Defino parametros con los que fueron obtenidos los datos a analizar:
k = 0.1234; % Constante elastica del cantilever en N/m
ny = 0.49; % Modulo de Poisson de la muestra (0.5 vivas, 0.3 materiales rigidos)
maxIndent = 1000*1e-9; % Indentacion maxima a considerar (sugerido 10% de la altura de la muestra) en metros
resolucion = 2048; % Cantidad de puntos (datos) por cada curva de fuerza

% Elijo tipo de indentador utilizado y modelo matematico para el analisis:
tip = 4; % 1: Esferico / 2: Parabolico / 3: Conico / 4: Piramidal
model = 3; % 1: Hertz / 2: Sneddon / 3: Sirghi
adhesion = 0; % 0: Calcula solo elasticidad / 1: Calcula elasticidad y adhesion

% Defino angulos y/o radio del indentador según corresponda al tipo utilizado:
% Indentador piramidal:
alfa = deg2rad(15); % Front Angle de la punta piramidal (en web del fabricante)
beta = deg2rad(25); % Back Angle de la punta piramidal (en web del fabricante)
gamma = deg2rad(17.5); % Side Angle de la punta piramidal (en web del fabricante)
ang = atan((tan(alfa)+tan(beta))*tan(gamma)); % Calculo del angulo que utiliza el modelo de indentador piramidal
% Indentador conico:
ang_cono = deg2rad(0); % Semi-angulo de apertura de la punta conica
% Indentador esferico o parabolico:
Rc = 1000*1e-9; % Radio de curvatura de la punta esferica o parabolica (en metros)

%% --- ANALISIS ---

curvasEx = zeros(resolucion,4); % Defino matriz generica para guardar fuerza e indentacion para curva de bajada (extend)
curvasRt = zeros(resolucion,4); % Defino matriz generica para guardar fuerza e indentacion para curva de subida (retrace)

for i=1:2
  
  % Leo los datos de un archivo de una curva de fuerza:
  archivo = [directorio getfield(archivos(i),'name')]; % Obtengo la ruta de cada archivo
  datos = dlmread(archivo,'\t',1,0); % Reescribo los datos en una matriz
  Calc_Ramp_Ex_nm = datos(:,1); % Altura, movimiento en Z del piezoelectrico durante la bajada (extend) en nm
  Calc_Ramp_Rt_nm = flipud(datos(:,2)); % Altura, movimiento en Z del piezoelectrico durante la subida (retract) en nm
  Defl_nm_Ex = datos(:,3); % Fuerza, deflexion del cantilever durante la bajada (extend) en nm
  Defl_nm_Rt = flipud(datos(:,4)); % Fuerza, deflexion del cantilever durante la subida (retract) en nm
  
  % Busco punto de contacto ("cero") manualmente, con ginput:
  figure, plot(Calc_Ramp_Ex_nm,Defl_nm_Ex,'linewidth',2) % Grafico curva de bajada (extend)
  hold on, plot(Calc_Ramp_Rt_nm,Defl_nm_Rt,'r','linewidth',2) % Grafico curva de subida (retract)
  xlabel('Z (nm)')
  ylabel('Deflection (nm)')
  [zCero,deflCero] = ginput(1); % Selecciono cero manualmente (punto a partir del que ambas curvas tienen comportamiento parabolico)
  
  % Corrijo ejes para que el 0 sea el punto de contacto:
  z_Ex = Calc_Ramp_Ex_nm - zCero; % Curva de bajada (extend)
  defl_Ex = Defl_nm_Ex - deflCero;
  z_Rt = Calc_Ramp_Rt_nm - zCero; % Curva de subida (retract)
  defl_Rt = Defl_nm_Rt - deflCero;
  indCero = find(z_Ex>0,1); % Indice del punto de contacto
  
  % Calculo la fuerza (resultado convertido a Newtons):
  f_Ex = (defl_Ex*1e-9)*k;
  f_Rt = (defl_Rt*1e-9)*k;
  
  % Determino la indentacion (resultado convertido a metros):
  indent_Ex = (z_Ex-defl_Ex)*1e-9;
  indent_Rt = (z_Rt-defl_Rt)*1e-9;
  raizIndent_Ex = sqrt(indent_Ex); % Raiz cuadrada de la indentacion (se usa en el modelo para indentador esferico o parabolico)
  indMaxIndent = find(indent_Ex>maxIndent,1); % Indice del punto de maxima indentacion a considerar
  
  % Guardo curvas corregidas en matrices:
  curvasEx(:,(2*i)-1) = indent_Ex; % Guardo indentacion para curva de bajada (extend)
  curvasEx(:,2*i) = f_Ex; % Guardo fuerza para curva de bajada (extend)
  curvasRt(:,(2*i)-1) = indent_Rt; % Guardo indentacion para curva de subida (retrace)
  curvasRt(:,2*i) = f_Rt; % Guardo fuerza para curva de subida (retrace)
  
  % Ajuste para hallar Modulo de Young (sobre la curva de bajada):
  if tip == 1 || tip == 2 % Indentador esferico o parabolico
    PE = polyfit(raizIndent_Ex(indCero:indMaxIndent),f_Ex(indCero:indMaxIndent),3); % Coeficientes del polinomio de ajuste de grado 3 para la curva de fuerza contra raiz de indentacion
  end
  if tip == 3 || tip == 4 % Indentador conico o piramidal
    PE = polyfit(indent_Ex(indCero:indMaxIndent),f_Ex(indCero:indMaxIndent),2); % Coeficientes del polinomio de ajuste de grado 2 para la curva de fuerza contra indentacion
  end
  PEeval = polyval(PE,indent_Ex(indCero:indMaxIndent)); % Evaluo polinomio de ajuste para los valores de indentacion a considerar
  rEeval = corrcoef(PEeval,f_Ex(indCero:indMaxIndent)); % Coeficientes de correlacion entre el polinomio de ajuste y la curva de fuerza contra indentacion
  r1 = rEeval(1,2); % Coeficiente de correlacion para el Modulo de Young
  
  % Ajuste para hallar Trabajo de Adhesion (sobre la curva de subida):
  if adhesion == 1
    if tip == 3 || tip == 4 % Indentador conico o piramidal
      PA = polyfit(indent_Rt(indCero:indMaxIndent),f_Rt(indCero:indMaxIndent),2); % Coeficientes del polinomio de ajuste de grado 2 para la curva de fuerza contra indentacion
    end
    PAeval = polyval(PA,indent_Rt(indCero:indMaxIndent)); % Evaluo polinomio de ajuste para los valores de indentacion a considerar
    rAeval = corrcoef(PAeval,f_Rt(indCero:indMaxIndent)); % Coeficientes de correlacion entre el polinomio de ajuste y la curva de fuerza contra indentacion
    r2 = rAeval(1,2); % Coeficiente de correlacion para el Trabajo de Adhesion
  end
  
  % Compruebo que hace el ajuste correctamente:
  if check == 1
    figure, plot(indent_Ex,f_Ex,'linewidth',2) % Grafico fuerza contra indentacion (extend)
    xlabel('Indentation (m)')
    ylabel('Force (N)')
    hold on, plot(indent_Rt,f_Rt,'r','linewidth',2) % Grafico fuerza contra indentacion (retract)
    plot(indent_Ex(indCero:indMaxIndent),PEeval,'k','linewidth',2) % Superpongo el ajuste para la indentacion considerada (extend)
    if adhesion == 1
      plot(indent_Rt(indCero:indMaxIndent),PAeval,'k','linewidth',2) % Superpongo el ajuste para la indentacion considerada (retract)
    end
  else
    close all % Cierro graficas
  end
  
  % Calculo Modulo de Young:
  if (tip == 1 || tip == 2) && model == 1
    Ei = PE(1)*3*(1-ny^2)/(4*Rc^(1/2)); % Indentador esferico o parabolico, modelo de Hertz
  end
  if tip == 3 && model == 2
    Ei = PE(1)*pi*(1-ny^2)/(2*tan(ang_cono)); % Indentador conico, modelo de Sneddon
  end
  if tip == 4 && model == 2
    Ei = PE(1)*sqrt(2)*(1-ny^2)/(tan(ang)); % Indentador piramidal, modelo de Sneddon
  end
  if tip == 4 && model == 3
    Ei = PE(1)*pi^(3/2)*(1-ny^2)/(4*tan(ang)); % Indentador piramidal, modelo de Sirghi 2008
  end
  % Si los datos vienen en N y m, E queda en Pa (para todos los modelos anteriores)
  Ei = Ei*1e-3; % Resultado convertido de Pa a kPa
  
  % Calculo Trabajo de Adhesion:
  if adhesion == 1
    if tip == 3 || tip == 4
      Ai=-PA(2)*pi^2*cos(ang)/(32*tan(ang)); % Indentador conico o piramidal, modelo de Sirghi 2008
    end
    % Si los datos vienen en N y m, A queda en J/m2 (joules sobre metro cuadrado)
  end
  
  % Muestro informacion en consola:
  disp(['Índice del punto de contacto: ' num2str(indCero)]) % Muestra el indice del punto de contacto
  disp(['Modulo de Young: ' num2str(Ei,'%.2f') ' kPa (correlacion: ' num2str(r1,'%.3f') ')']) % Muestra resultado de E en consola
  if adhesion == 1
    disp(['Trabajo de Adhesion: ' num2str(Ai,'%.2e') ' J/m2 (correlacion: ' num2str(r2,'%.3f') ')']) % Muestra resultado de A en consola
  end
  disp(' ');
  
end

%% --- GRAFICO ---

figure, plot(curvasEx(:,1),curvasEx(:,2),'b','linewidth',3)
hold on, plot(curvasRt(:,1),curvasRt(:,2),':b','linewidth',3)
plot(curvasEx(:,3),curvasEx(:,4),'r','linewidth',3)
plot(curvasRt(:,3),curvasRt(:,4),':r','linewidth',3), hold off

xlabel({'Indentation (m)'},'FontSize',12)
ylabel({'Force (N)'},'FontSize',12)
xlim([-1.5e-6 2e-6]) % Modificar limites del eje x segun mejore la visualizacion
% Modificar leyendas del grafico:
legend({'Control (approach)','Control (withdrawal)','Modified (approach)','Modified (withdrawal)'},'FontSize',12)

warning('on', 'all'); % Vuelvo a habilitar los warnings