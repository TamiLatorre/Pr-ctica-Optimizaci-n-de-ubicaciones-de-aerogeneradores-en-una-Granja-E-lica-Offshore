% Los datos de viento constan solamente de las componentes u y v de la
% velocidad a una altura de 100metros. Esto a lo largo de una grilla la
% cual debe ser filtrada para obteer aquellos puntos mas cercanos al area
% de emplazamiento. 

% Con las variables u,v se calcula la magnitud de velocidad y el angulo
% azimuth de la direccion del viento. 

% Porsteriormente al tener unos 6 puntos coordenados (longitud, latitud)
% con datos mensuales desde 1940-2022 se calcula un promedio en cada
% instante de tiempo para esos puntos obteniendo un vector de direccion y
% magnitud de velocidad con los cuales calculan los parametros de Weibull y
% el sector de frecuencia. 

% En un ultimo paso, se obtienen nuevamente parametros de Weibull y el
% sector de frecuencia para el caso del mejor 40% de los datos de viento. 

%% leyendo datos de viento

% Archivo en formato .nc 
info = ncinfo('datostotal.nc');

% Para leer una varible
u_100m = ncread('datostotal.nc', 'u100');
v_100m = ncread('datostotal.nc', 'v100');
longitud = ncread('datostotal.nc', 'longitude');
latitud = ncread('datostotal.nc', 'latitude');
tiempo = ncread('datostotal.nc', 'time');
ex = ncread('datostotal.nc', 'expver');

%% Cálculo de azimuth

% Se elimina los datos de expver (experimentales correspondientes a expver[5])
v_final=squeeze(v_100m(:,:,1,:));
u_final=squeeze(u_100m(:,:,1,:));

% Obtenemos la direccion del viento en cada punto del tiempo (1 dato por mes)
% formato en grados, 0 es el norte y aumentando en sentido horario.(direccion a la que se dirige el viento)
azimuth=NaN(size(u_final));
[indx] =find((u_final>=0) & (v_final>=0));
    azimuth(indx)= rad2deg(atan(u_final(indx)./v_final(indx)));
[indx] =find((u_final>=0) & (v_final<0));
    azimuth(indx)=rad2deg(pi+atan(u_final(indx)./v_final(indx)));
[indx] =find((u_final<0) & (v_final<0));
    azimuth(indx)=rad2deg(pi+atan(u_final(indx)./v_final(indx)));
[indx] =find((u_final<0) & (v_final>=0));
    azimuth(indx)=rad2deg(2*pi+atan(u_final(indx)./v_final(indx)));
    
%% filtrando cerca del area de emplazamiento 

coordenadas =  [longitud, latitud];

% mayor -36.540058, menor -36.714829
% mayor -73.463581, menor -73.555744 

% Encontramos los valores que esten en este rango
lat_f= find(coordenadas(:,2)>=-36.75 & coordenadas(:,2)<=-36.5);
lon_f = find(coordenadas(:,1)>=-73.75 & coordenadas(:,1)<=-73.25);

% Crea la matriz u_filtrado con los valores dentro de los rangos de longitud y latitud
u_filtrado = u_final(lon_f, lat_f, :);
v_filtrado = v_final(lon_f, lat_f, :);
azimuth_filtrado = azimuth(lon_f, lat_f, :);

% Calcular la magnitud de la velocidad en cada instante de tiempo 
magnitud_vel =sqrt(u_filtrado.^2 + v_filtrado.^2);


% Datos promedio para cada instante tiempo
promedio_vel_tiempo = mean(magnitud_vel,[1 2],'omitnan');
promedio_azi_tiempo = mean(azimuth_filtrado,[1 2],'omitnan');

% La funcion indica hacia adonde apunta la direccion del viento, en este
% caso le podemos sumar 180 si queremos ver desde donde proviene.
wind_rose(promedio_azi_tiempo+180,promedio_vel_tiempo);

% Promedio para todos los puntos en el a traves del tiempo
promediovel=mean(promedio_vel_tiempo,'all','omitnan'); 
% Promedio de cada punto coordenado para todos los tiempos
promedio_puntos = squeeze(mean(magnitud_vel, 3,'omitnan'));

wd= squeeze(promedio_azi_tiempo(1,1,:)); % datos de magnitud de velocidad
ws= squeeze(promedio_vel_tiempo(1,1,:)); % datos de direccion 

%% Parametros de Weibull caso 1


% Se utiliza la funcion wblfit de MATLAB

[paramEsts2, paramCIs] = wblfit(ws(1:1006));
% paramEsts contiene el parametro de escala (a), y el segundo es el
% parametro de forma (k)
scaleParameter2 = paramEsts2(1);
shapeParameter2 = paramEsts2(2);

wd2=wd(1:1006)+180 % desde donde proviene el viento y ultimos dos datos que son NaN

numBins = 24; 
histogram(wd2, numBins);
numSectors = 24;  
sectorCounts = histcounts(wd2, linspace(0, 360, numSectors+1));
% Para visualizar lo obtenido 
x = linspace(min(ws(1:1006)), max(ws(1:1006)), 100);
y = wblpdf(x, scaleParameter2, shapeParameter2);
histogram(ws(1:1006), numBins, 'Normalization', 'pdf');
hold on;
legend('Datos','pdf')
plot(x, y, 'r', 'LineWidth', 2);
xlabel('Velocidad [m/s]')
ylabel('Probabilidad')
title('Ajuste de distribucion de Weibull')
hold off;


%% sector de frecuencia

numSectors = 36;
% Calcula los límites de los sectores
sectorLimits = linspace(0, 360, numSectors+1);

% Calcula la frecuencia de ocurrencia en cada sector
[f, edges] = histcounts(wd2, sectorLimits, 'Normalization', 'probability');
% 'f' es el vector de sector de frecuencias que se utilizara en Python

% Visualiza el vector de frecuencias normalizado junto con los límites de los sectores
disp('Probabilidad en cada sector y límites:');
for i = 1:numSectors
    fprintf('Sector %d: Probabilidad = %.2f, Límites: %.0f° a %.0f°\n', i, f(i), edges(i), edges(i+1));
end

% Visualiza el histograma de dirección del viento con sectores normalizado
figure;
histogram(wd2, sectorLimits, 'Normalization', 'probability');
xlabel('Dirección del Viento (grados)');
ylabel('Probabilidad');
title('Histograma de Dirección del Viento con Sectores Normalizado');


% Serie de timpo de la velocidad 
%plot(squeeze(promedio_vel_tiempo(1,1,:)))


%% caso 2: 50% mejor de los datos de viento 

velocidad_caso2 = mean(magnitud_vel,[1 2],'omitnan');
c2_vel =squeeze(velocidad_caso2(1,1,1:1006));
azimuth_caso2 = mean(azimuth_filtrado,[1 2],'omitnan');
c_azi =squeeze(azimuth_caso2(1,1,1:1006));
c2_azi = c_azi+180 ;% desde donde viene el viento

% Filtrar el 40% mejor de mag de vel, encontrar esas ubicaciones y usarlas
% y tambien su direccion

% Calcula el percentil 60 
percentil_60 = prctile(c2_vel, 60);
% Selecciona las velocidades que son mayores o iguales al percentil 50
c2_40_v = c2_vel(c2_vel >= percentil_60);
c2_40_a = c2_azi(c2_vel >= percentil_60);
% El 40%mejor del viento

% Rosa de los vientos
wind_rose(c2_40_a,c2_40_v)

% Parametros de Weibull
[paramEsts2, paramCIs] = wblfit(c2_40_v);
scaleParameter40 = paramEsts2(1);
shapeParameter40 = paramEsts2(2);

%
numBins = 24; 
histogram(c2_40_a, numBins);
numSectors = 24;  
sectorCounts = histcounts(c2_40_a, linspace(0, 360, numSectors+1));
% Visualizamos 
x = linspace(min(c2_40_v), max(c2_40_v), 100);
y = wblpdf(x, scaleParameter40, shapeParameter40);
histogram(c2_40_v, numBins, 'Normalization', 'pdf');
hold on;
plot(x, y, 'r', 'LineWidth', 2);
xlabel('Velocidad [m/s]')
ylabel('Probabilidad')
title('Ajuste de distribucion de Weibull')

% Sector de frecuencias
numSectors = 36;
sectorLimits = linspace(0, 360, numSectors+1);
% Calcula la frecuencia de ocurrencia en cada sector
[f, edges] = histcounts(c2_40_a, sectorLimits, 'Normalization', 'probability');
% Visualiza el vector de frecuencias normalizado junto con los límites de los sectores
disp('Probabilidad en cada sector y límites:');
for i = 1:numSectors
    fprintf('Sector %d: Probabilidad = %.2f, Límites: %.0f° a %.0f°\n', i, f(i), edges(i), edges(i+1));
end
% Visualiza el histograma de dirección del viento con sectores normalizado
figure;
histogram(c2_40_a, sectorLimits, 'Normalization', 'probability');
xlabel('Dirección del Viento (grados)');
ylabel('Probabilidad');
title('Histograma de Dirección del Viento con Sectores Normalizado');
