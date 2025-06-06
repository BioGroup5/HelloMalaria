function solo_parasitos = leerAnotacionesParasitos(txtPath)
    fileID = fopen(txtPath, 'r');
    lines = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    lines = lines{1};
    lines(1) = [];
    data = cell(length(lines), 7);
    for i = 1:length(lines)
        tokens = strsplit(lines{i}, ',');
        for j = 1:7
            if j <= length(tokens), data{i,j} = tokens{j}; else, data{i,j} = ''; end
        end
    end
    dataTable = cell2table(data, 'VariableNames', {'ID','Label','Comment','Shape','NumPoints','X1','Y1'});
    dataTable.X1 = str2double(dataTable.X1);
    dataTable.Y1 = str2double(dataTable.Y1);
    solo_parasitos = dataTable(ismember(dataTable.Label, {'Parasite', 'Parasitized'}), :);
end