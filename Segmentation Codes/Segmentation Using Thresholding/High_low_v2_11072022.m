% Find the mean image prior to thresholding ------ JUNE 2023 ----------
% Author: Basil Thurakkal, TIFR Hyderabad --- Contact: basilt@tifrh.res.in 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% update in v2 - user input files. Keep the files in the same folder as the code. 

bin_file = uigetfile('*.tif','Select bin file');
fprintf(bin_file)
display('\n')
int_file = uigetfile('*.tif','Select int file');
fprintf(int_file)

bin_img = imread(bin_file);
bin_img = im2bw(bin_img); bin_img = 1-bin_img;
int_img = imread(int_file);

%%
bin_img = imcomplement(bin_img);
cc = bwconncomp(bin_img,4);%connected component
%%
cell_prop = regionprops(cc,int_img,'Area','MeanIntensity'); %regionprops
%area filter

%% Displaying image with mean intensity
Mean_image = zeros(size(bin_img,1),size(bin_img,2));
Cells = cc.NumObjects; 
Pixel_loc = regionprops(cc,'PixelList');

    for i = 1:Cells
        area = cell_prop(i).Area;
        intensity = cell_prop(i).MeanIntensity;
        tot_intensity = area*intensity;

        for j = 1:area
            Location = [Pixel_loc(i).PixelList(j,1) Pixel_loc(i).PixelList(j,2)];
            Mean_image(Location(2),Location(1)) = intensity;
        end
    end
    imagesc(Mean_image);
%% Displaying image with total intensity
figure()
Tot_image = zeros(size(bin_img,1),size(bin_img,2));
Cells = cc.NumObjects;
Pixel_loc = regionprops(cc,'PixelList');

    for i = 1:Cells
        area = cell_prop(i).Area;
        intensity = cell_prop(i).MeanIntensity;
        tot_intensity = area*intensity;

        for j = 1:area
            Location = [Pixel_loc(i).PixelList(j,1) Pixel_loc(i).PixelList(j,2)];
            Tot_image(Location(2),Location(1)) = tot_intensity;
        end
    end
    imagesc(Tot_image);