% LOCAL DENSITY,INTENSITY CORRELATTION CALCULATION PROGRAM ------  JUNE 2023 ----------
% Author: Basil Thurakkal, TIFR Hyderabad --- Contact: basilt@tifrh.res.in 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
%Cell_seg_and_intens_calc %Segment the image and find average intensity value per cell

%% !!IMPORTANT!! : Use this section instead of line 3, if you already have segmented mask from cellpose,etc.Otherwise keep it commented out.
% Comment out line 3 if this section is in use.

TMRM_img = imread(''); %import TMRM intensity image
mask_img = imread(''); %import mask
segmented_img_analysis;
%%
centroid_image % Creates cell-centroid image for local density calculation
Loca_Density_TMRM_corr %Finds local density and local intensity fields 
correlation_int_density %finds intensity at high and low density population
plot_fig %output the files in csv format

