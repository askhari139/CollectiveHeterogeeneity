

mask = im2bw(mask_img);
mask = 1-mask;
cc = bwconncomp(mask,4);
multiply =TMRM_img;
region_prop = regionprops(cc,multiply,'MeanIntensity')

%% show mean image
Mean_image = zeros(size(multiply));
Cells = cc.NumObjects;
Pixel_loc = regionprops(cc,'PixelList');

    for i = 1:Cells
        area = Areas(i).Area;
        intensity = region_prop(i).MeanIntensity;
        for j = 1:area
            Location = [Pixel_loc(i).PixelList(j,1) Pixel_loc(i).PixelList(j,2)];
            Mean_image(Location(2),Location(1)) = intensity;
        end
    end
    imagesc(Mean_image);
%%