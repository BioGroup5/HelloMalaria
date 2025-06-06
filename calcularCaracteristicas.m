function tabla = calcularCaracteristicas(mascara, imagen, nombreImagen, tipo)
    props = regionprops(mascara, 'Area', 'Perimeter', 'PixelIdxList');
    n = length(props);
    if n == 0
        tabla = table();
        return;
    end

    R = imagen(:,:,1);
    G = imagen(:,:,2);
    B = imagen(:,:,3);

    Area = [props.Area]';
    Perimetro = [props.Perimeter]';
    Circularidad = zeros(n, 1);
    for i = 1:n
        A = props(i).Area;
        P = props(i).Perimeter;
        if P > 0 && A > 0
            Circularidad(i) = 4 * pi * A / (P^2);
        else
            Circularidad(i) = 0;
        end
    end


    IntensidadR = zeros(n, 1);
    IntensidadG = zeros(n, 1);
    IntensidadB = zeros(n, 1);

    for i = 1:n
        pixIdx = props(i).PixelIdxList;
        IntensidadR(i) = mean(R(pixIdx));
        IntensidadG(i) = mean(G(pixIdx));
        IntensidadB(i) = mean(B(pixIdx));
    end

    IntensidadRGBMedia = mean([IntensidadR, IntensidadG, IntensidadB], 2);

    NombreImagen = repmat({nombreImagen}, n, 1);
    Tipo = repmat({tipo}, n, 1);

    tabla = table(NombreImagen, Tipo, Area, Perimetro, Circularidad, ...
                  IntensidadR, IntensidadG, IntensidadB, IntensidadRGBMedia);
end