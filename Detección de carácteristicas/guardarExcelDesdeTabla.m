function guardarExcelDesdeTabla(tablaDatos, carpetaGuardar)
    nombreArchivo = 'datasetCaracteristicasGlobalThresholding.csv';
    rutaCompleta = fullfile(carpetaGuardar, nombreArchivo);

    if exist(rutaCompleta, 'file')
        fprintf('El archivo ya existe en:\n%s\n', rutaCompleta);
        fprintf('Los datos no se volverán a guardar.\n');
    else
        writetable(tablaDatos, rutaCompleta); 
        fprintf('Archivo CSV guardado en:\n%s\n', rutaCompleta);
    end
end