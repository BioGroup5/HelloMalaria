function [mascaraGB, mascaraParasitosFinal, mascaraFalsosPositivos] = procesarImagenGlobal(imgPath, txtPath)
    IPF = imread(imgPath);
    Imagen_FILT = im2double(IPF);    
    I_gray = rgb2gray(Imagen_FILT);

    Imagen_mascara_gray = rgb2gray(cat(3, IPF(:,:,1), IPF(:,:,2), zeros(size(IPF(:,:,1)))));
    fov_mask = segmentCellsStainBased(Imagen_mascara_gray);

    OD = -log(Imagen_FILT + 0.01);
    OD = reshape(OD, [], 3);
    W_manual = [
        0.650, 0.072, 0.707;
        0.704, 0.990, 0.707;
        0.286, 0.117, 0.000
    ];
    W = W_manual ./ vecnorm(W_manual, 2, 2);
    W = W';
    C = OD * inv(W);
    deconvolved_manual = reshape(C, size(IPF));
    Imagen_nucleos = deconvolved_manual(:,:,3);

    bw = imbinarize(Imagen_nucleos);
    Imagen_masked = bw .* fov_mask;

    imagenEcualizada = adapthisteq(I_gray);
    mascaraPortaobjetos = imbinarize(imagenEcualizada, graythresh(imagenEcualizada));
    mascaraPortaobjetos = imfill(mascaraPortaobjetos, 'holes');
    mascaraPortaobjetos = bwareafilt(mascaraPortaobjetos, 1);
    imagenEcualizada(~mascaraPortaobjetos) = 0;
    pixelesDentroFOV = imagenEcualizada(fov_mask);
    umbralOscuridad = prctile(pixelesDentroFOV, 5);
    mascaraOscuros = imagenEcualizada < umbralOscuridad;
    mascaraOscuros = mascaraOscuros & logical(fov_mask);
    areaMinimaGlobulos = 1800;
    mascaraGB = bwareafilt(mascaraOscuros, [areaMinimaGlobulos Inf]);

    elementoEstructurante = strel('disk', 15);
    mascaraGBdilatada = imdilate(mascaraGB, elementoEstructurante);

    imagenSinGB = Imagen_masked;
    imagenSinGB(mascaraGBdilatada) = 0;
    imagenSegmentable = uint8(imagenSinGB);

    umbralGlobal = graythresh(imagenSegmentable);
    mascaraParasitosYRuido = imbinarize(imagenSegmentable, umbralGlobal);
    areaMinimaParasitos = 400;
    mascaraParasitosYRuidoFiltrada = filtrar_por_tamano(mascaraParasitosYRuido, areaMinimaParasitos);

    solo_parasitos = leerAnotacionesParasitos(txtPath);
    x_real = solo_parasitos{:,6};
    y_real = solo_parasitos{:,7};
    propsObjetos = regionprops(mascaraParasitosYRuidoFiltrada, 'Centroid', 'Area', 'Perimeter', 'PixelIdxList');

    indicesAsignados = zeros(length(x_real), 1);
    for i = 1:length(x_real)
        centroides = reshape([propsObjetos.Centroid], 2, [])';
        distancias = sqrt((centroides(:,1) - x_real(i)).^2 + (centroides(:,2) - y_real(i)).^2);
        [~, idxMin] = min(distancias);
        indicesAsignados(i) = idxMin;
    end

    mascaraParasitosFinal = false(size(mascaraParasitosYRuidoFiltrada));
    for i = 1:length(indicesAsignados)
        objeto = propsObjetos(indicesAsignados(i));
        mascaraParasitosFinal(objeto.PixelIdxList) = true;
    end

    mascaraParasitos = mascaraParasitosFinal;
    mascaraInvertida = ~mascaraParasitos;
    mascaraFalsosPositivos = mascaraParasitosYRuidoFiltrada & mascaraInvertida;
end