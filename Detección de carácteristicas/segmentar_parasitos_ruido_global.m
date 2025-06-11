function [mascaraParasitos, mascaraRuido, propsParasitos, propsRuido] = ...
         segmentar_parasitos_ruido_global(imagenSoloGlobulos, fov_mask, anotaciones)
    % Convierte imagen a gris para segmentación
    if size(imagenSoloGlobulos,3) == 3
        gris = rgb2gray(imagenSoloGlobulos);
    else
        gris = imagenSoloGlobulos;
    end

    % Binarizar imagen con Otsu
    nivel = graythresh(gris);
    bin = imbinarize(gris, nivel);

    % Dilatar máscara de glóbulos blancos para evitar zonas cercanas
    elem = strel('disk', 15);
    binDil = imdilate(bin, elem);

    % Crear imagen sin glóbulos blancos dilatados dentro del FOV
    imagenSinGlobulos = double(fov_mask);
    imagenSinGlobulos(binDil) = 0;
    imagenSinGlobulos = uint8(imagenSinGlobulos);

    % Umbral global para segmentar parásitos + ruido
    nivelGlobal = graythresh(imagenSinGlobulos);
    mascaraParasitoyRuido = imbinarize(imagenSinGlobulos, nivelGlobal);

    % Filtrar objetos pequeños (posible ruido)
    areaMinimaParasitos = 400;
    mascaraParasitoyRuido = bwareafilt(mascaraParasitoyRuido, [areaMinimaParasitos Inf]);

    % Detectar parásitos circulares
    [centrosParasitos, radiosParasitos] = imfindcircles(mascaraParasitoyRuido, [2 5], ...
                                                       'Sensitivity', 0.92, 'EdgeThreshold', 0.5);

    % Inicializar máscaras vacías
    mascaraParasitos = false(size(mascaraParasitoyRuido));

    if isempty(anotaciones)
        mascaraRuido = mascaraParasitoyRuido;
        propsParasitos = regionprops(false(size(mascaraParasitoyRuido)));
        propsRuido = regionprops(mascaraRuido, 'Area', 'Perimeter', 'PixelIdxList');
        return;
    end

    radioExtra = 20;

    for i = 1:height(anotaciones)
        puntoReal = [anotaciones.X1(i), anotaciones.Y1(i)];
        for j = 1:size(centrosParasitos,1)
            centroDet = centrosParasitos(j,:);
            radioDet = radiosParasitos(j) + radioExtra;
            distancia = norm(puntoReal - centroDet);
            if distancia <= radioDet
                [X, Y] = meshgrid(1:size(mascaraParasitoyRuido,2), 1:size(mascaraParasitoyRuido,1));
                mascaraCirculo = (X - centroDet(1)).^2 + (Y - centroDet(2)).^2 <= radioDet^2;
                mascaraParasitos = mascaraParasitos | mascaraCirculo;
                break;
            end
        end
    end

    % Ruido = máscara global menos máscara de parásitos
    mascaraRuido = mascaraParasitoyRuido & ~mascaraParasitos;

    % Obtener propiedades
    propsParasitos = regionprops(mascaraParasitos, 'Area', 'Perimeter', 'PixelIdxList');
    propsRuido = regionprops(mascaraRuido, 'Area', 'Perimeter', 'PixelIdxList');
end