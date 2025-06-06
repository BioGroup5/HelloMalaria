function fov_mask = segmentCellsStainBased(I)
    if size(I, 3) == 3
        I = rgb2gray(I);
    end
    level = graythresh(I);
    bw = imbinarize(I, level);
    if mean(I(bw)) < mean(I(~bw))
        bw = ~bw;
    end
    bw_clean = imopen(bw, strel('disk', 5));
    bw_clean = imclose(bw_clean, strel('disk', 5));
    bw_clean = imfill(bw_clean, 'holes');
    stats = regionprops(bw_clean, 'Area', 'Centroid', 'EquivDiameter');
    [~, idx] = max([stats.Area]);
    circle_centroid = stats(idx).Centroid;
    circle_radius = stats(idx).EquivDiameter / 2;
    adjusted_radius = max(circle_radius - 15, 0);
    [rows, cols] = size(bw_clean);
    [X, Y] = meshgrid(1:cols, 1:rows);
    fov_mask = (X - circle_centroid(1)).^2 + (Y - circle_centroid(2)).^2 <= adjusted_radius^2;
end