%% Force-indentation curves analysis - Young's Modulus - Work of Adhesion
% Este programa analiza todas las curvas de fuerza de un directorio
% Esta version (completa) obtiene Modulo de Young a partir de la curva de bajada (extend)
% y Trabajo de Adhesion a partir de la curva de subida (retract)
close all
clear
clc

%% Defino parametros de ubicacion de los datos a analizar:

dirScript = pwd; % Defino directorio actual donde esta el programa
dirDatos = uigetdir('C:\Users\ruta_del_archivo\'); % Defino directorio donde estan los datos
dirDatos = [dirDatos '\'];
cd (dirDatos)
arch = 'Archivo0*'; % Defino nombre generico de archivos a abrir
cd (dirScript)

%% Defino parametros con los que fueron obtenidos los datos a analizar:

k = 0.1422; % Constante elastica del cantilever en N/m
ny = 0.49; % Modulo de Poisson de la muestra (0.5 vivas, 0.3 materiales rigidos)
maxIndent = 1000*1e-9; % Indentacion maxima a considerar (sugerido 10% de la altura de la muestra) en metros
% Indentador piramidal:
alfa = deg2rad(15); % Front Angle de la punta piramidal (en web del fabricante)
beta = deg2rad(25); % Back Angle de la punta piramidal (en web del fabricante)
gamma = deg2rad(17.5); % Side Angle de la punta piramidal (en web del fabricante)
ang = atan((tan(alfa)+tan(beta))*tan(gamma)); % Calculo del angulo que utiliza el modelo de indentador piramidal
% Indentador conico:
%ang = deg2rad(); % Semi-angulo de apertura de la punta conica
% Indentador esferico o parabolico:
%Rc = ; % Radio de curvatura de la punta esferica o parabolica

%% Compruebo que encuentra los archivos con los datos:

D = dir([dirDatos arch]);
nD = length(D);
if nD == 0
  disp('No se encontraron archivos');
end

%% Inicializo vectores para guardar los resultados:

E = zeros(nD,1); % Defino vector generico para guardar resultados de Modulo de Young
re = zeros(nD,1); % Defino vector generico para guardar coeficientes de correlacion de ajuste para elasticidad
A = zeros(nD,1); % Defino vector generico para guardar resultados de Trabajo de Adhesion
ra = zeros(nD,1); % Defino vector generico para guardar coeficientes de correlacion de ajuste para adhesion

%% Analizo cada curva de fuerza:

threshold = 0.98; % Defino umbral de tolerancia del coeficiente de correlacion para el modulo de Young

for i = 1:nD
    
    % Leo los datos de un archivo de una curva de fuerza:
    datos = getfield(D(i),'name');
    datos = [dirDatos datos]; % Obtengo la ruta del archivo especifico
    M = dlmread(datos,'\t',1,0); % Reescribo los datos en una matriz
      Calc_Ramp_Ex_nm = M(:,1); % Altura, movimiento en Z del piezoelectrico durante la bajada (extend) en nm
      Calc_Ramp_Rt_nm = flipud(M(:,2)); % Altura, movimiento en Z del piezoelectrico durante la subida (retract) en nm
      Defl_nm_Ex = M(:,3); % Fuerza, deflexion del cantilever durante la bajada (extend) en nm
      Defl_nm_Rt = flipud(M(:,4)); % Fuerza, deflexion del cantilever durante la subida (retract) en nm
    
    % Pruebo analisis con distintos puntos de contacto hasta hallar el mejor ajuste:
    indCero = 1; % Inicializo el indice del punto de contacto para entrar al bucle
    r1 = 0; % Inicializo el coeficiente de correlacion para entrar al bucle
    while r1 < threshold
        
        % Determino el punto de contacto:
        zCero = Calc_Ramp_Ex_nm(indCero);
        deflCero = Defl_nm_Ex(indCero);

        % Corrijo ejes para que el 0 sea el punto de contacto:
        z_Ex = Calc_Ramp_Ex_nm - zCero; % Curva de bajada (extend)
        defl_Ex = Defl_nm_Ex - deflCero;
        z_Rt = Calc_Ramp_Rt_nm - zCero; % Curva de subida (retract)
        defl_Rt = Defl_nm_Rt - deflCero;
        
        % Calculo la fuerza (resultado convertido a Newtons):
        f_Ex = (defl_Ex*1e-9)*k;
        f_Rt = (defl_Rt*1e-9)*k;
        
        % Determino la indentacion (resultado convertido a metros):
        indent_Ex = (z_Ex-defl_Ex)*1e-9;
        indent_Rt = (z_Rt-defl_Rt)*1e-9;
        raizIndent_Ex = sqrt(indent_Ex); % Raiz cuadrada de la indentacion (se usa en el modelo para indentador esferico o parabolico)
        indMaxIndent = find(indent_Ex>maxIndent,1); % Indice del punto de maxima indentacion a considerar
        
        % Ajuste para hallar Modulo de Young (sobre la curva de bajada):
        % Indentador piramidal o conico:
        PE = polyfit(indent_Ex(indCero:indMaxIndent),f_Ex(indCero:indMaxIndent),2); % Coeficientes del polinomio de ajuste de grado 2 para la curva de fuerza contra indentacion
        % Indentador esferico o parabolico:
        %PE = polyfit(raizIndent_Ex(indCero:indMaxIndent),f(indCero:indMaxIndent),3); % Coeficientes del polinomio de ajuste de grado 3 para la curva de fuerza contra raiz de indentacion
        PEeval = polyval(PE,indent_Ex(indCero:indMaxIndent)); % Evaluo polinomio de ajuste para los valores de indentacion a considerar
        rEeval = corrcoef(PEeval,f_Ex(indCero:indMaxIndent)); % Coeficientes de correlacion entre el polinomio de ajuste y la curva de fuerza contra indentacion
        r1 = rEeval(1,2); % Coeficiente de correlacion para el Modulo de Young
        
        % Compruebo correlacion o busco un nuevo punto de contacto:
        if r1 < threshold
            indCero = indCero + 1;
            if indCero > length(z_Ex)
                break; % Corto el bucle si no hay mas puntos para probar
            end
        end
    end
    
    % Ajuste para hallar Trabajo de Adhesion (sobre la curva de subida):
    % Indentador piramidal o conico:
    PA = polyfit(indent_Rt(indCero:indMaxIndent),f_Rt(indCero:indMaxIndent),2); % Coeficientes del polinomio de ajuste de grado 2 para la curva de fuerza contra indentacion
    PAeval = polyval(PA,indent_Rt(indCero:indMaxIndent)); % Evaluo polinomio de ajuste para los valores de indentacion a considerar
    rAeval = corrcoef(PAeval,f_Rt(indCero:indMaxIndent)); % Coeficientes de correlacion entre el polinomio de ajuste y la curva de fuerza contra indentacion
    r2 = rAeval(1,2); % Coeficiente de correlacion para el Trabajo de Adhesion
    
    % Compruebo que hace el ajuste correctamente:
    %figure, plot(indent_Ex,f_Ex,'linewidth',2) % Grafico fuerza contra indentacion (extend)
    %xlabel('Indentation (m)')
    %ylabel('Force (N)')
    %hold on, plot(indent_Rt,f_Rt,'r','linewidth',2) % Grafico fuerza contra indentacion (retract)
    %plot(indent_Ex(indCero:indMaxIndent),PEeval,'k','linewidth',2) % Superpongo el ajuste para la indentacion considerada (extend)
    %plot(indent_Rt(indCero:indMaxIndent),PAeval,'k','linewidth',2) % Superpongo el ajuste para la indentacion considerada (retract)
    
    % Calculo Modulo de Young (elegir ecuacion del modelo que corresponda al indentador):
    Ei = PE(1)*pi^(3/2)*(1-ny^2)/(4*tan(ang)); % Indentador piramidal, modelo de Sirghi 2008
    %Ei = PE(1)*sqrt(2)*(1-ny^2)/(tan(ang)); % Indentador piramidal, modelo de Sneddon
    %Ei = PE(1)*pi*(1-ny^2)/(2*tan(ang)); % Indentador conico, modelo de Sneddon
    %Ei = PE(1)*3*(1-ny^2)/(4*Rc^(1/2)); % Indentador esferico o parabolico, modelo de Hertz
    % Si los datos vienen en N y m, E queda en Pa (para todos los modelos anteriores)
    Ei = Ei*1e-3; % Resultado convertido de Pa a kPa
    
    % Calculo Trabajo de Adhesion (modelo de Sirghi, 2008):
    Ai=-PA(2)*pi^2*cos(ang)/(32*tan(ang)); % Indentador piramidal o conico
    % Si los datos vienen en N y m, A queda en J/m2 (joules sobre metro cuadrado)
    
    % Guardo resultado en vectores, muestro resultado para la curva analizada y cierro graficas:
    E(i) = Ei; % Guardo Modulo de Young en vector
    re(i) = r1; % Guardo coeficiente de correlacion para E
    A(i) = Ai; % Guardo Trabajo de Adhesion en vector
    ra(i) = r2; % Guardo coeficiente de correlacion para A
    disp(['Modulo de Young: ' num2str(Ei,'%.2f') ' kPa (correlacion: ' num2str(r1,'%.3f') ')']) % Muestra resultado de E en consola
    disp(['Trabajo de Adhesion: ' num2str(Ai,'%.2e') ' J/m2 (correlacion: ' num2str(r2,'%.3f') ')']) % Muestra resultado de A en consola
    disp(['Quedan ' num2str(nD-i) ' curvas...']) % Muestra cuantas curvas faltan por analizar
    close all
end

%% Filtro, guardo y muestro los resultados finales:

E = E(re>=threshold&E>0); % Vector solo con los E positivos con coeficiente de correlacion mayor al umbral de tolerancia definido
A = A(ra>=0.95&A>0); % Vector solo con los A positivos con coeficiente de correlacion mayor a 0.95
save ([dirDatos 'Resultado Modulo de Young'],'E','-ascii'); % Guardo el vector de resultados de E en formato ascii
save ([dirDatos 'Resultado Trabajo de Adhesion'],'A','-ascii'); % Guardo el vector de resultados de A en formato ascii
disp(['Modulo de Young promedio para esta celula: ' num2str(mean(E),'%.2f') ' kPa'])
disp(['Trabajo de Adhesion promedio para esta celula: ' num2str(mean(A),'%.2e') ' J/m2'])