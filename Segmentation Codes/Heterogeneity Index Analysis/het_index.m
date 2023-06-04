% HETEROGENEITY INDEX CALCULATION PROGRAM ------ v1.0; JUNE 2023 ----------
% Author: Basil Thurakkal, TIFR Hyderabad --- Contact: basilt@tifrh.res.in 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% heterogeneity index function
% inputs = int image, mask image (as array)
% output = het index value (normaalized difference between mean of top 25% intensity cells
% and bottom 25% intensity cells

function [het_i,var_het,reg_prop] = het_index(int,mask)
%% segment and make mean image
TMRM = int; % TMRM image
mask = mask; % mask
cc = bwconncomp(mask,4);% connected component
reg_prop = regionprops(cc,TMRM,'MeanIntensity','PixelList','Area','Centroid')% regionprop (Mean intensity );

%% show mean image
    
    Mean_image = zeros(size(TMRM));
    Cells = cc.NumObjects;
    Pixel_loc = regionprops(cc,'PixelList');
    for i = 1:Cells
        area = reg_prop(i).Area;
        intensity = reg_prop(i).MeanIntensity;
        for j = 1:area
            Location = [Pixel_loc(i).PixelList(j,1) Pixel_loc(i).PixelList(j,2)];
            Mean_image(Location(2),Location(1)) = intensity;
        end
    end
    imagesc(Mean_image);
    
%% Top and bottom intensity calculation

int_array = [reg_prop(:).MeanIntensity]; %struct field to array
low_int = prctile(int_array,25);
high_int = prctile(int_array,75);% find the 25% and 75% threshold

T = struct2table(reg_prop); % region prop struct to table
T(T.Area<100,:) = []; % Area filter
idx_low = T.MeanIntensity<low_int; 
idx_high = T.MeanIntensity>high_int; 
low_table = T(idx_low,:); % bottom 25% table
high_table = T(idx_high,:);% top 25% table

mean_low = mean(low_table.MeanIntensity);
mean_high = mean(high_table.MeanIntensity);

var_low = var(low_table.MeanIntensity);
var_high = var(high_table.MeanIntensity);

%% heterogneity index
het_i = (mean_high-mean_low);
var_het = var_high+var_low;
end

