function datasetFinal = procesarCarpetaGlobal(folderPath)
    imageFiles = [ 
        dir(fullfile(folderPath, '**', '*.jpg')); 
        dir(fullfile(folderPath, '**', '*.tiff'))
    ];

    datasetFinal = table();

    for k = 1:length(imageFiles)
        imgFile = imageFiles(k);
        imgPath = fullfile(imgFile.folder, imgFile.name);
        [~, name, ~] = fileparts(imgFile.name);
        txtPath = fullfile(imgFile.folder, [name, '.txt']);
    
        if exist(txtPath, 'file')
            fprintf('Procesando imagen: %s\n', imgPath);
            try
                [mascaraGB, mascaraParasitos, mascaraRuido] = procesarImagenGlobal(imgPath, txtPath);
                IPF = imread(imgPath);
    
                tablaGB = calcularCaracteristicas(mascaraGB, IPF, imgFile.name, 'GlobuloBlanco');
                tablaParasitos = calcularCaracteristicas(mascaraParasitos, IPF, imgFile.name, 'Parasito');
                tablaRuido = calcularCaracteristicas(mascaraRuido, IPF, imgFile.name, 'Ruido');
    
                caracteristicas = [tablaGB; tablaParasitos; tablaRuido];
                datasetFinal = [datasetFinal; caracteristicas];
            catch ME
                warning('Error procesando %s: %s', imgPath, ME.message);
            end
        else
            warning('No se encontró archivo .txt para %s', imgPath);
        end
    end
end