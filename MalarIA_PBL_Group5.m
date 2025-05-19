clear; close all; clc;

% --- Cambia solo estos nombres ---
image_filename = '20170829_174225.jpg';
txt_filename = '20170829_174225.txt';
% ---------------------------------

% Leer la imagen
IPF = imread(image_filename);
img = IPF;

% Leer todas las líneas del archivo de texto
fileID = fopen(txt_filename, 'r');
lines = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
lines = lines{1};

% Eliminar la primera línea (no contiene datos relevantes)
lines(1) = [];

% Preparar estructuras para los datos
numLines = length(lines);
data = cell(numLines, 7);  % Solo 7 columnas para la tabla
x2 = NaN(numLines, 1);     % X2 para círculos
y2 = NaN(numLines, 1);     % Y2 para círculos

% Procesar cada línea
for i = 1:numLines
    tokens = strsplit(lines{i}, ',');
    % Llenar columnas fijas
    for j = 1:7
        if j <= length(tokens)
            data{i, j} = tokens{j};
        else
            data{i, j} = '';
        end
    end
    % Leer X2/Y2 si existen
    if length(tokens) >= 9
        x2(i) = str2double(tokens{8});
        y2(i) = str2double(tokens{9});
    end
end

% Convertir a tabla con solo 7 columnas
dataTable = cell2table(data, ...
    'VariableNames', {'ID', 'Label', 'Comment', 'Shape', 'NumPoints', 'X1', 'Y1'});

% Convertir columnas numéricas
dataTable.NumPoints = str2double(dataTable.NumPoints);
dataTable.X1 = str2double(dataTable.X1);
dataTable.Y1 = str2double(dataTable.Y1);

% Mostrar cantidad de cada tipo
parasites = strcmpi(dataTable.Label, 'Parasite') | strcmpi(dataTable.Label, 'Parasitized');
white_cells = strcmpi(dataTable.Label, 'White_Blood_Cell');

fprintf('Número de parásitos: %d\n', sum(parasites));
fprintf('Número de glóbulos blancos: %d\n', sum(white_cells));

% Mostrar imagen
figure;
imshow(IPF);
hold on;

% Dibujar los elementos
for i = 1:height(dataTable)
    shape = lower(dataTable.Shape{i});
    x1 = dataTable.X1(i); y1 = dataTable.Y1(i);

    if strcmp(shape, 'circle') && all(~isnan([x2(i), y2(i)]))
        c = ([x1, y1] + [x2(i), y2(i)]) / 2;
        r = norm([x2(i)-x1, y2(i)-y1]) / 2;
        viscircles(c, r, 'Color', 'r');
    elseif strcmp(shape, 'point')
        plot(x1, y1, 'y+', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
end

hold off;
title('Anotaciones en la imagen');
%%
% Quitar el ruido gaussiano
IPF_gauss_filt=imgaussfilt(IPF,2,'FilterSize',5);

% Quitar el ruido salt&pepper
IPF_salt_pepper_filt=zeros(size(IPF_gauss_filt), 'uint8');

    for k = 1:3
        IPF_salt_pepper_filt(:,:,k) = medfilt2(IPF_gauss_filt(:,:,k), [3 3]);
    end

Imagen_FILT = im2double(IPF_salt_pepper_filt);
Imagen_FILT(Imagen_FILT==0)=eps;

% Convertir a OD (densidad óptica)
IPF_double = im2double(IPF);
OD = -log(IPF_double + 0.01);
OD = reshape(OD, [], 3);  % Matriz N x 3

% Definir manualmente la matriz de canales W (3 vectores de tinción):
W_manual = [
    0.650, 0.072, 0.707;   % Canal 1 (núcleo)
    0.704, 0.990, 0.707;   % Canal 2 (citoplasma)
    0.286, 0.117, 0.000    % Canal 3 (residual)
];

% Normalizar los vectores para que sean unitarios
W = W_manual ./ vecnorm(W_manual, 2, 2);

% Transponer para usar como base de tinción
W = W';  % Ahora es 3x3, como se espera

% Proyectar OD a los nuevos canales definidos
C = OD * inv(W);  % Inversa de la matriz base

% Reconstruir imagen deconvolucionada
deconvolved_manual = reshape(C, size(IPF));

% Visualizar cada canal
figure;
for i = 1:3
    subplot(1, 3, i);
    imshow(mat2gray(deconvolved_manual(:,:,i)));
    switch i
        case 1
            title('Canal 1 - Residual');
        case 2
            title('Canal 2 - Citoplasma');
        case 3
            title('Canal 3 - Núcleo');
    end
end
subplot;

Imagen_nucleos=deconvolved_manual(:,:,3);
Imagen_mascara=cat(3,IPF(:,:,1), IPF(:,:,2));
fov_mask=segmentCellsStainBased(Imagen_mascara);

bw = imbinarize(Imagen_nucleos);
Imagen_masked=bw.*fov_mask;
imshow(Imagen_masked);title('Imagen de núcleos solo con circulo de interes')
%%
% Preprocesamiento
nucleo_bin = imbinarize(imadjust(deconvolved_manual(:,:,3))); % mejorar contraste y binarizar

% Quitar objetos pequeños
nucleo_clean = bwareaopen(nucleo_bin, 20);

% Extraer propiedades de regiones candidatas
props = regionprops(nucleo_clean, 'Centroid', 'Area', 'Eccentricity');

% Filtrar posibles parásitos (regiones pequeñas, redondas, etc.)
parasite_candidates = [];
for i = 1:length(props)
    if props(i).Area < 200 && props(i).Eccentricity < 0.8
        parasite_candidates(end+1,:) = props(i).Centroid;
    end
end

% Suponiendo que el archivo tiene las coordenadas y etiquetas
dos = readtable(txt_filename); % Asegúrate de que el nombre y formato son correctos
% Ver etiquetas únicas
disp(unique(dataTable{:,2}));

% Filtrar parásitos reales
solo_parasitos = dataTable(strcmp(dataTable{:,2}, 'Parasite'), :);

% Mostrar la imagen original
imshow(IPF); title('Parásitos detectados');
hold on;

% Dibujar candidatos automáticos
if ~isempty(parasite_candidates)
plot(parasite_candidates(:,1), parasite_candidates(:,2), 'r+', 'MarkerSize', 5);
end

% Dibujar parásitos reales si existen y coordenadas son válidas
if ~isempty(solo_parasitos) && all(size(solo_parasitos,2) >= 7)
x_real = solo_parasitos{:,6};
y_real = solo_parasitos{:,7};

% Verificar que estén dentro del rango de la imagen
[h, w, ~] = size(IPF);
validos = x_real > 0 & x_real < w & y_real > 0 & y_real < h;

plot(x_real(validos), y_real(validos), 'go', 'MarkerSize', 5, 'LineWidth', 1.5);
end
hold off;
%%
nucleo = mat2gray(Imagen_nucleos);
figure;
subplot(1,2,1); imshow(nucleo); title('Canal de núcleos');
subplot(1,2,2); imhist(nucleo); title('Histograma');
%%
% --- DETECTAR GLOBULOS BLANCOS Y CREAR MÁSCARA ---
% Convertir a escala de grises
I_gray = rgb2gray(Imagen_FILT);

% Mejorar contraste y suavizar
I_eq = adapthisteq(I_gray);
I_filt = medfilt2(I_eq, [3 3]);

% Aplicar máscara del portaobjetos (si es necesario)
BW_circle = imbinarize(I_filt, graythresh(I_filt));
BW_circle = imfill(BW_circle, 'holes');
BW_circle = bwareafilt(BW_circle, 1);
I_filt(~BW_circle) = 0;

% Visualizar imagen filtrada e histograma
figure;
imshow(I_filt);
title('Imagen con fondo eliminado');

figure;
imhist(I_filt);
title('Histograma de intensidad');
disp(['Tipo de datos de I_filt: ', class(I_filt)]);

% Umbral por intensidad
umbral_bajo = 0.35;  
mascara_oscuros = I_filt < umbral_bajo;

% Ajustable: tamaño mínimo de glóbulos blancos
tam_minimo = 1500;  
mascara_final = bwareafilt(mascara_oscuros, [tam_minimo Inf]);

% Mostrar resultado
figure;
imshow(mascara_final);
title(['Glóbulos blancos (área mínima: ', num2str(tam_minimo), ' px)']);
%%
% Globulo zurixak dilata
% Maskari aplika umbralekuai
% Behin globuko zurixak kendutakun parasituk detekta.
%%
% DESAROLLO DE 4 MÉTODOS DE SEGMENTACIÓN

% Aplicar máscara para eliminar glóbulos blancos de la imagen
Imagen_masked_sin_GB = Imagen_masked;
Imagen_masked_sin_GB(mascara_final) = 0;  % Rellenar glóbulos blancos con negro (0)

% Mostrar imagen sin glóbulos blancos
figure;
imshow(Imagen_masked_sin_GB);
title('Imagen sin glóbulos blancos');

% Asegurar que la imagen está en escala de grises
if size(Imagen_masked_sin_GB, 3) == 3
    Imagen_masked_sin_GB = rgb2gray(Imagen_masked_sin_GB);
end

imagen_segmentable = uint8(Imagen_masked_sin_GB);

% MODELO 1: LOCAL THRESHOLDING
bw_local = imbinarize(Imagen_masked_sin_GB, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.6); 
figure;
imshow(bw_local); title('Umbral Local (Adaptativo)');

% Filtrado morfológico para eliminar objetos pequeños (glóbulos blancos residuales, ruido)
umbral_tamano = 400;
bw_filtrado = filtrar_por_tamano(bw_local, umbral_tamano);

% Visualizar resultado del filtrado
figure;
imshow(bw_filtrado);
title('Parásitos filtrados (sin glóbulos blancos)');

% Aplicar dilatación morfológica para engrosar los objetos
se = strel('disk', 1);  
mascara_dilatada = imdilate(bw_filtrado, se);

% Visualizar resultado del filtrado y dilatado
figure;
imshow(mascara_dilatada);
title('Parásitos filtrados y dilatados');

% Detección circular de parásitos
[centros_parasitos, radios_parasitos] = imfindcircles(mascara_dilatada, [2 5], 'Sensitivity', 0.92, 'EdgeThreshold', 0.5);

% Contar parásitos detectados automáticamente y reales
num_auto = size(centros_parasitos, 1);
num_reales = size(solo_parasitos, 1);
num_total = num_auto + num_reales;

% Mostrar resultados sobre imagen original
figure; imshow(IPF);
title('Detección de parásitos');
hold on;

% Dibujar círculos detectados (automáticos)
viscircles(centros_parasitos, radios_parasitos, 'EdgeColor', 'r');

% Dibujar parásitos reales (desde archivo)
plot(solo_parasitos{:,6}, solo_parasitos{:,7}, 'go', 'MarkerSize', 5, 'LineWidth', 1.5);

% Mostrar cantidades en texto
text(10, 20, ['Detectados (rojo): ' num2str(num_auto)], 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
text(10, 40, ['Reales (verde): ' num2str(num_reales)], 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold');
text(10, 60, ['Total: ' num2str(num_total)], 'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'black');

hold off;

% También mostrar en consola
fprintf('Detectados: %d\nReales: %d\nTotal: %d\n', num_auto, num_reales, num_total);

%----------------------------------

% Parámetro para ampliar el radio de los círculos si es necesario
radio_extra = 15;  % Aumentar la tolerancia del radio (puedes ajustarlo)

% Verificar cuántos parásitos reales caen dentro de algún círculo detectado
aciertos = 0;
for i = 1:size(solo_parasitos, 1)
    punto_real = [solo_parasitos{i,6}, solo_parasitos{i,7}];
    
    for j = 1:size(centros_parasitos, 1)
        centro_detectado = centros_parasitos(j,:);
        radio_detectado = radios_parasitos(j) + radio_extra;  % Ampliamos el radio de detección

        % Calcular distancia entre el punto real y el centro del círculo
        distancia = norm(punto_real - centro_detectado);

        % Verificar si el punto real está dentro del círculo ampliado
        if distancia <= radio_detectado
            aciertos = aciertos + 1; % Incrementar aciertos si está dentro
            break;  % Si se encuentra uno, no es necesario buscar más círculos para este punto real
        end
    end
end

% Mostrar en consola
fprintf('Parásitos correctamente detectados: %d de %d reales\n', aciertos, size(solo_parasitos,1));

% Mostrar en la imagen también (opcional)
text(10, 80, ['Aciertos: ' num2str(aciertos)], 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');

%----------------------------------

%%
% MODELO 2: GLOBAL THRESHOLDING (OTSU)

level_global = graythresh(Imagen_masked_sin_GB);
bw_otsu = imbinarize(Imagen_masked_sin_GB, level_global);
figure;
imshow(bw_otsu); title('Umbral Global (Otsu)');

% Filtrado morfológico para eliminar objetos pequeños
umbral_tamano = 400;
bw_filtrado = filtrar_por_tamano(bw_otsu, umbral_tamano);

% Visualizar resultado del filtrado
figure;
imshow(bw_filtrado);
title('Parásitos filtrados (sin glóbulos blancos)');

% Aplicar dilatación morfológica para engrosar los objetos
se = strel('disk', 1);  % Puedes ajustar el tamaño del disco
mascara_dilatada = imdilate(bw_filtrado, se);

% Visualizar resultado del filtrado y dilatado
figure;
imshow(mascara_dilatada);
title('Parásitos filtrados y dilatados');

% Detección circular de parásitos
[centros_parasitos, radios_parasitos] = imfindcircles(mascara_dilatada, [2 5], 'Sensitivity', 0.92, 'EdgeThreshold', 0.5);

% Contar parásitos detectados automáticamente y reales
num_auto = size(centros_parasitos, 1);
num_reales = size(solo_parasitos, 1);
num_total = num_auto + num_reales;

% Mostrar resultados sobre imagen original
figure; imshow(IPF);
title('Detección de parásitos (Otsu)');
hold on;

% Dibujar círculos detectados (automáticos)
viscircles(centros_parasitos, radios_parasitos, 'EdgeColor', 'r');

% Dibujar parásitos reales (desde archivo)
plot(solo_parasitos{:,6}, solo_parasitos{:,7}, 'go', 'MarkerSize', 5, 'LineWidth', 1.5);

% Mostrar cantidades en texto
text(10, 20, ['Detectados (rojo): ' num2str(num_auto)], 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
text(10, 40, ['Reales (verde): ' num2str(num_reales)], 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold');
text(10, 60, ['Total: ' num2str(num_total)], 'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'black');

hold off;

% También mostrar en consola
fprintf('Detectados: %d\nReales: %d\nTotal: %d\n', num_auto, num_reales, num_total);

%----------------------------------

% Parámetro para ampliar el radio de los círculos si es necesario
radio_extra = 15;  % Aumentar la tolerancia del radio (puedes ajustarlo)

% Verificar cuántos parásitos reales caen dentro de algún círculo detectado
aciertos = 0;
for i = 1:size(solo_parasitos, 1)
    punto_real = [solo_parasitos{i,6}, solo_parasitos{i,7}];
    
    for j = 1:size(centros_parasitos, 1)
        centro_detectado = centros_parasitos(j,:);
        radio_detectado = radios_parasitos(j) + radio_extra;  % Ampliamos el radio de detección

        % Calcular distancia entre el punto real y el centro del círculo
        distancia = norm(punto_real - centro_detectado);

        % Verificar si el punto real está dentro del círculo ampliado
        if distancia <= radio_detectado
            aciertos = aciertos + 1; % Incrementar aciertos si está dentro
            break;  % Si se encuentra uno, no es necesario buscar más círculos para este punto real
        end
    end
end

% Mostrar en consola
fprintf('Parásitos correctamente detectados: %d de %d reales\n', aciertos, size(solo_parasitos,1));

% Mostrar en la imagen también (opcional)
text(10, 80, ['Aciertos: ' num2str(aciertos)], 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');

%----------------------------------
%%
% MODELO 3: K-MEANS

[num_clases, centros] = imsegkmeans(imagen_segmentable, 3);

% Crear visualización con superposición
imagen_agrupada = labeloverlay(imagen_segmentable, num_clases);

figure;
imshow(imagen_agrupada);
title('Segmentación por K-means (k = 3)');

% Extraer máscaras binarias para cada clase
mascara_1 = num_clases == 1;
mascara_2 = num_clases == 2;
mascara_3 = num_clases == 3;

% Mostrar clusters individualmente
figure;
montage({mascara_1, mascara_2, mascara_3}, 'Size', [1 3]);
title('Clusters individuales');

% Se considera que el cluster 2 (del citoplasma) contiene los parásitos
mascara_parasitos = mascara_2;

% Filtrado morfológico para eliminar objetos pequeños
umbral_tamano = 400;
mascara_filtrada = filtrar_por_tamano(mascara_parasitos, umbral_tamano);

% Visualizar resultado del filtrado
figure;
imshow(mascara_filtrada);
title('Parásitos filtrados (sin glóbulos blancos)');

% Detección circular de parásitos
[centros_parasitos, radios_parasitos] = imfindcircles(mascara_filtrada, [2 5], 'Sensitivity', 0.92, 'EdgeThreshold', 0.5);

% Contar parásitos detectados automáticamente y reales
num_auto = size(centros_parasitos, 1);
num_reales = size(solo_parasitos, 1);
num_total = num_auto + num_reales;

% Mostrar resultados sobre imagen original
figure; imshow(IPF);
title('Detección de parásitos (K-means)');
hold on;

% Dibujar círculos detectados (automáticos)
viscircles(centros_parasitos, radios_parasitos, 'EdgeColor', 'r');

% Dibujar parásitos reales (desde archivo)
plot(solo_parasitos{:,6}, solo_parasitos{:,7}, 'go', 'MarkerSize', 5, 'LineWidth', 1.5);

% Mostrar cantidades en texto
text(10, 20, ['Detectados (rojo): ' num2str(num_auto)], 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
text(10, 40, ['Reales (verde): ' num2str(num_reales)], 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold');
text(10, 60, ['Total: ' num2str(num_total)], 'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'black');

hold off;

% También mostrar en consola
fprintf('Detectados: %d\nReales: %d\nTotal: %d\n', num_auto, num_reales, num_total);

%----------------------------------

% Parámetro para ampliar el radio de los círculos si es necesario
radio_extra = 15;  % Aumentar la tolerancia del radio (puedes ajustarlo)

% Verificar cuántos parásitos reales caen dentro de algún círculo detectado
aciertos = 0;
for i = 1:size(solo_parasitos, 1)
    punto_real = [solo_parasitos{i,6}, solo_parasitos{i,7}];
    
    for j = 1:size(centros_parasitos, 1)
        centro_detectado = centros_parasitos(j,:);
        radio_detectado = radios_parasitos(j) + radio_extra;  % Ampliamos el radio de detección

        % Calcular distancia entre el punto real y el centro del círculo
        distancia = norm(punto_real - centro_detectado);

        % Verificar si el punto real está dentro del círculo ampliado
        if distancia <= radio_detectado
            aciertos = aciertos + 1; % Incrementar aciertos si está dentro
            break;  % Si se encuentra uno, no es necesario buscar más círculos para este punto real
        end
    end
end

% Mostrar en consola
fprintf('Parásitos correctamente detectados: %d de %d reales\n', aciertos, size(solo_parasitos,1));

% Mostrar en la imagen también (opcional)
text(10, 80, ['Aciertos: ' num2str(aciertos)], 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');

%----------------------------------
%%
% Convertir a escala de grises si es necesario
if size(IPF, 3) == 3
    I_gray = rgb2gray(IPF);
else
    I_gray = IPF;
end

% Aplicar K-means con 4 clusters
num_clusters = 4;
[clases_kmeans, centros] = imsegkmeans(I_gray, num_clusters);

% Calcular intensidad media por cluster
intensidades_cluster = zeros(1, num_clusters);
for k = 1:num_clusters
    intensidades_cluster(k) = mean(I_gray(clases_kmeans == k));
end

% Ordenar intensidades y obtener el segundo más oscuro
[~, orden] = sort(intensidades_cluster);
cluster_GB = orden(2); % segundo más oscuro

% Crear máscara binaria de glóbulos blancos
mascara_GB_kmeans = clases_kmeans == cluster_GB;

% Filtrar por tamaño
tam_min = 1500;
tam_max = 100000;
mascara_GB_kmeans_filtrada = bwareafilt(mascara_GB_kmeans, [tam_min tam_max]);

% Visualizar resultado
figure; imshow(IPF); hold on;
props_GB_kmeans = regionprops(mascara_GB_kmeans_filtrada, 'BoundingBox', 'Centroid');
for i = 1:length(props_GB_kmeans)
    rectangle('Position', props_GB_kmeans(i).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 1.5);
    c = props_GB_kmeans(i).Centroid;
    text(c(1), c(2), num2str(i), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
end
title('Glóbulos blancos detectados con K-means (k = 4)');
hold off;
% --- Comparación entre detección automática y anotaciones reales ---

% Cantidad detectada automáticamente con K-means
cantidad_GB_detectados = length(props_GB_kmeans);

% Cantidad real desde el archivo .txt
cantidad_GB_reales = sum(strcmpi(dataTable.Label, 'White_Blood_Cell'));

% Diferencia
diferencia_GB = cantidad_GB_detectados - cantidad_GB_reales;

% Mostrar resultados
fprintf('\n--- Comparación de glóbulos blancos ---\n');
fprintf('Detectados automáticamente (K-means): %d\n', cantidad_GB_detectados);
fprintf('Anotados en el archivo .txt: %d\n', cantidad_GB_reales);
fprintf('Diferencia: %+d\n', diferencia_GB);
%%
% --- Cálculo de métricas por imagen ---
ID_imagen = image_filename;

% 1. Glóbulos blancos

props_GB = regionprops(mascara_GB_kmeans_filtrada, IPF(:,:,1), 'Area', 'MeanIntensity', 'BoundingBox', 'Centroid');
areas_GB = [props_GB.Area];
intensidades_GB = [props_GB.MeanIntensity];
intensidad_media_GB = mean(intensidades_GB);
cantidad_GB = sum(white_cells);

% Filtrar áreas anómalas (muy grandes)
media_area = mean(areas_GB);
desv_area = std(areas_GB);
umbral_max = media_area + 2 * desv_area;

% Filtrar solo las áreas válidas
areas_filtradas = areas_GB(areas_GB < umbral_max);

% Calcular área media solo con las válidas
area_media_GB = mean(areas_filtradas);


% Mostrar áreas detectadas de glóbulos blancos con etiquetas numéricas
figure;
imshow(IPF);
hold on;

% Dibujar rectángulos y etiquetas
for k = 1:length(props_GB)
    % Dibujar rectángulo
    rectangle('Position', props_GB(k).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 1.5);
    
    % Añadir número en el centroide
    c = props_GB(k).Centroid;
    text(c(1), c(2), num2str(k), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
end
title('Glóbulos blancos detectados y numerados con k-means');
hold off;

% 2. Parásitos (solo los reales correctamente detectados)
mascara_parasitos_reales_detectados = false(size(mascara_dilatada));

for i = 1:size(solo_parasitos, 1)
    punto_real = [solo_parasitos{i,6}, solo_parasitos{i,7}];
    for j = 1:size(centros_parasitos, 1)
        centro_detectado = centros_parasitos(j,:);
        radio_detectado = radios_parasitos(j) + radio_extra;
        distancia = norm(punto_real - centro_detectado);
        if distancia <= radio_detectado
            [X, Y] = meshgrid(1:size(mascara_dilatada,2), 1:size(mascara_dilatada,1));
            mascara_parasitos_reales_detectados = mascara_parasitos_reales_detectados | ...
                ((X - punto_real(1)).^2 + (Y - punto_real(2)).^2 <= 5^2); % radio 5 píxeles
            break;
        end
    end
end

% Extraer características solo de los parásitos reales detectados
props_parasitos = regionprops(mascara_parasitos_reales_detectados, IPF(:,:,1), ...
    'Area', 'Perimeter', 'MeanIntensity');

% Calcular métricas si hay detecciones válidas
if ~isempty(props_parasitos)
    areas_parasitos = [props_parasitos.Area];
    intensidades_parasitos = [props_parasitos.MeanIntensity];
    perimetros_parasitos = [props_parasitos.Perimeter];
    circularidades = 4 * pi * areas_parasitos ./ (perimetros_parasitos.^2);
    area_media_parasitos = mean(areas_parasitos);
    intensidad_media_parasitos = mean(intensidades_parasitos);
    circularidad_media_parasitos = mean(circularidades);
else
    area_media_parasitos = NaN;
    intensidad_media_parasitos = NaN;
    circularidad_media_parasitos = NaN;
end

figure;
imshow(mascara_parasitos_reales_detectados);
title('Máscara de parásitos reales correctamente detectados');


fprintf('Área media de parásitos reales detectados: %.2f\n', area_media_parasitos);
fprintf('Intensidad media: %.2f\n', intensidad_media_parasitos);
fprintf('Circularidad media: %.2f\n', circularidad_media_parasitos);

% 3. Ruido (objetos detectados por Otsu que no son parásitos reales)
mascara_ruido = mascara_dilatada & ~mascara_parasitos_reales_detectados;

% Extraer características del ruido
props_ruido = regionprops(mascara_ruido, IPF(:,:,1), 'Area', 'Perimeter', 'MeanIntensity');

% Calcular métricas si hay ruido
if ~isempty(props_ruido)
    intensidades_ruido = [props_ruido.MeanIntensity];
    intensidad_media_ruido = mean(intensidades_ruido);
else
    intensidad_media_ruido = NaN;
end
figure;
imshow(mascara_ruido);
title('Máscara de ruido (segmentos no reales)');

fprintf('Intensidad media del ruido: %.2f\n', intensidad_media_ruido);
% 4. Diagnostico
% Determinar diagnóstico a partir del archivo .txt
labels_presentes = unique(dataTable.Label);

if any(strcmpi(labels_presentes, 'Parasite'))
    diagnostico = "PF";
elseif any(strcmpi(labels_presentes, 'Parasitized'))
    diagnostico = "PV";
else
    diagnostico = "Non-Infected";
end

cantidad_GB = sum(white_cells);
cantidad_parasitos = sum(parasites);

% --- Crear tabla con los resultados de la imagen actual ---
resultados_nuevos = table({ID_imagen}, intensidad_media_GB, area_media_GB, cantidad_GB, ...
    intensidad_media_parasitos, area_media_parasitos, circularidad_media_parasitos, ...
    cantidad_parasitos, intensidad_media_ruido, diagnostico, ...
    'VariableNames', {'ID_Imagen', 'IntensidadMedia_GB', 'AreaMedia_GB', 'Cantidad_GB', ...
    'IntensidadMedia_Parasitos', 'AreaMedia_Parasitos', 'CircularidadMedia_Parasitos', ...
    'Cantidad_Parasitos', 'IntensidadMedia_Ruido', 'Diagnostico'});


% --- Nombre del archivo Excel ---
% Ruta relativa que sube un nivel desde la carpeta actual
filename_excel = fullfile('..', 'resultados_malaria.xlsx');

% --- Verificar si el archivo ya existe ---
if isfile(filename_excel)
    % Leer tabla existente
    tabla_existente = readtable(filename_excel);

    % Verificar si el ID ya está en la tabla
    if any(strcmp(tabla_existente.ID_Imagen, ID_imagen))
        disp(['La imagen "', ID_imagen, '" ya está en la tabla. No se añadirá de nuevo.']);
    else
        % Añadir nueva fila y guardar
        tabla_actualizada = [tabla_existente; resultados_nuevos];
        writetable(tabla_actualizada, filename_excel, 'FileType', 'spreadsheet');
        disp(['Se ha añadido la imagen "', ID_imagen, '" a "', filename_excel, '".']);
    end
else
    % Crear nuevo archivo Excel con la primera fila
    writetable(resultados_nuevos, filename_excel, 'FileType', 'spreadsheet');
    disp(['Archivo "', filename_excel, '" creado con la primera imagen.']);
end
%%
function fov_mask  = segmentCellsStainBased(I)

    if size(I,3) == 2
        stain1 = I(:,:,1);
        stain2 = I(:,:,2);
       
        stain_combined = (stain1 + stain2) / 2;
        stain_gray = mat2gray(stain_combined);
       
        binary_mask = stain_gray > 0.15;

        try

            [centers, radii] = imfindcircles(binary_mask, [floor(min(size(binary_mask))/4) floor(min(size(binary_mask))/2)], ...
                'ObjectPolarity', 'bright', 'Sensitivity', 0.9);
           
            if ~isempty(centers) && ~isempty(radii)
                [y, x] = ndgrid(1:size(binary_mask, 1), 1:size(binary_mask, 2));

                fov_mask = sqrt((x - centers(1,1)).^2 + (y - centers(1,2)).^2) <= (radii(1) * 0.9);
            else      
                throw(MException('FallbackToTraditional', 'No circle detected'));
            end
        catch

            props = regionprops(binary_mask, 'Centroid', 'Area');
            if ~isempty(props)

                center = round(props(1).Centroid);
                estimated_radius = sqrt(props(1).Area / pi);
 
                [y, x] = ndgrid(1:size(binary_mask, 1), 1:size(binary_mask, 2));

                fov_mask = sqrt((x - center(1)).^2 + (y - center(2)).^2) <= (estimated_radius * 0.9);
            else

                fov_mask = imclose(binary_mask, strel('disk', 20));
                fov_mask = imfill(fov_mask, 'holes');

                cc = bwconncomp(fov_mask);
                if cc.NumObjects > 0
                    numPixels = cellfun(@numel, cc.PixelIdxList);
                    [~, indexOfMax] = max(numPixels);
                    fov_mask = false(size(fov_mask));
                    fov_mask(cc.PixelIdxList{indexOfMax}) = true;
                end

                fov_mask = imerode(fov_mask, strel('disk', 10));
               
            end
        end
    end
end

function salida = filtrar_por_tamano(imagen_binaria, umbral_max)
    componentes = bwlabel(imagen_binaria);
    stats = regionprops(componentes, 'Area');
    salida = ismember(componentes, find([stats.Area] <= umbral_max));
end