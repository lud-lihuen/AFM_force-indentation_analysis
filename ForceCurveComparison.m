%% Force curve comparison
% Este programa grafica dos curvas de fuerza de bajada (extend) en un mismo par de ejes
% Las curvas deben estar en un mismo directorio
clear all
close all
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
%alfa = ; % Semi-angulo de apertura de la punta conica
% Indentador esferico o parabolico:
%Rc = ; % Radio de curvatura de la punta esferica o parabolica

%% Compruebo que encuentra los archivos con los datos:

D = dir([dirDatos arch]);
nD = length(D);
if nD == 0
  disp('No se encontraron archivos');
end

%% Analizo cada curva de fuerza:

E = zeros(nD,1); % Defino vector generico para guardar resultados de Modulo Elastico
re = zeros(nD,1); % Defino vector generico para guardar coeficientes de correlacion
curvasDefl = zeros(512,2*nD); % Defino matriz generica para guardar deflexion y altura
curvasIndent = zeros(512,2*nD); % Defino matriz generica para guardar fuerza e indentacion

for i=1:nD
    
    % Leo los datos de un archivo de una curva de fuerza:
    datos = getfield(D(i),'name');
    datos = [dirDatos datos]; % Obtengo la ruta del archivo especifico
    M = dlmread(datos,'\t',1,0); % Reescribo los datos en una matriz
      Calc_Ramp_Ex_nm = M(:,1); % Altura, movimiento en Z del piezoelectrico durante la bajada (extend) en nm
      Defl_nm_Ex = M(:,3); % Fuerza, deflexion del cantilever durante la bajada (extend) en nm
      curvasDefl(:,(2*i)-1) = M(:,1); % Guardo altura
      curvasDefl(:,2*i) = M(:,3); % Guardo deflexion
    
    % Busco punto de contacto ("cero") en la curva de bajada, manualmente, con ginput:
    figure, plot(Calc_Ramp_Ex_nm,Defl_nm_Ex,'linewidth',2) % Grafico curva de fuerza (extend)
    xlabel('Z (nm)')
    ylabel('Deflection (nm)')
    [zCero,deflCero] = ginput(1); % Selecciono cero manualmente
    
    % Corrijo ejes para que el 0 sea el punto de contacto:
    z = Calc_Ramp_Ex_nm - zCero;
    defl_Ex = Defl_nm_Ex - deflCero;
    indCero = find(z>0,1); % Indice del punto de contacto
    
    % Determino la indentacion:
    indent = (z-defl_Ex)*1e-9; % Resultado convertido a metros
    raizIndent = sqrt(indent); % Raiz cuadrada de la indentacion (se usa en el modelo para indentador esferico o parabolico)
    curvasIndent(:,(2*i)-1) = indent; % Guardo indentacion (si la punta es piramidal o conica)
    %curvasIndent(:,(2*i)-1) = raizIndent; % Guardo raiz cuadrada de la indentacion (si la punta es esferica o parabolica)
    indMaxIndent = find(indent>maxIndent,1); % Indice del punto de maxima indentacion a considerar
    
    % Calculo la fuerza:
    f = (defl_Ex*1e-9)*k; % Resultado convertido a Newtons
    curvasIndent(:,2*i) = f; % Guardo fuerza
    
    % Ajuste cuadratico para hallar Modulo Elastico:
    % Indentador piramidal o conico:
    P = polyfit(indent(indCero:indMaxIndent),f(indCero:indMaxIndent),2); % Coeficientes del polinomio de ajuste de grado 2 para la curva de fuerza contra indentacion 
    % Indentador esferico o parabolico:
    %P = polyfit(raizIndent(indCero:indMaxIndent),f(indCero:indMaxIndent),3); % Coeficientes del polinomio de ajuste de grado 3 para la curva de fuerza contra raiz de indentacion
    pev = polyval(P,indent(indCero:indMaxIndent)); % Evaluo polinomio de ajuste para los valores de indentacion a considerar
    ro = corrcoef(pev,f(indCero:indMaxIndent)); % Coeficientes de correlacion entre el polinomio de ajuste y la curva de fuerza contra indentacion
    r1 = ro(1,2); % Coeficiente de correlacion para el Modulo Elastico
    
    % Controlo que hace el ajuste correctamente:
    figure, plot(indent,f,'linewidth',2) % Grafico fuerza contra indentacion   
    xlabel('Indentation (m)')
    ylabel('Force (N)')
    hold on, plot(indent(indCero:indMaxIndent),pev,'r','linewidth',2) % Superpongo el ajuste cuadratico para la indentacion considerada
    
    % Calculo Modulo Elastico (elegir ecuacion del modelo que corresponda al indentador):
    Ei = P(1)*pi^(3/2)*(1-ny^2)/(4*tan(ang)); % Indentador piramidal, modelo de Sirghi 2008
    %Ei = P(1)*2^(1/2)*(1-ny^2)/(tan(ang)); % Indentador piramidal, modelo de Bilodeau
    %Ei = P(1)*pi*(1-ny^2)/(2*tan(alfa)); % Indentador conico, modelo de Sneddon
    %Ei = P(1)*3*(1-ny^2)/(4*Rc^(1/2)); % Indentador esferico o parabolico, modelo de Hertz
    
    % Si los datos vienen en N y m, E queda en Pa (para todos los modelos anteriores)
    Ei = Ei*1e-3 % Resultado convertido de Pa a kPa
    
    % Guardo resultados en vectores y cierro graficas:
    E(i) = Ei; % Guardo Modulo Elastico en vector
    re(i) = r1; % Guardo coeficiente de correlacion para E
    disp(['Quedan ' num2str(nD-i) ' curvas...']) % Imprime en pantalla cuantas curvas faltan por analizar
    close all
end

%% Grafico las curvas de fuerza en el mismo par de ejes:

% Grafico deflexion contra altura:
%figure, plot(curvasDefl(:,1),curvasDefl(:,2),'linewidth',2)
%xlabel('Z (nm)')
%ylabel('Deflection (nm)')
%hold on, plot(curvasDefl(:,3),curvasDefl(:,4),'r','linewidth',2)

% Grafico fuerza contra indentacion:
figure, plot(curvasIndent(:,1),curvasIndent(:,2),'r','linewidth',2)
xlabel('Indentation (m)')
ylabel('Force (N)')
xlim([-2e-6 2e-6]) % Elijo limites del eje X que permitan visualizar mejor las curvas
hold on, plot(curvasIndent(:,3),curvasIndent(:,4),'linewidth',2)