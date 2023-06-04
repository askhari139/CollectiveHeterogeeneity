%% Finds local intensity and local density fields. 

Cell_centroid_matrix = Im1;
average_mean_image = Mean_image/mean(mean(Mean_image));
Local_Density = {};
Local_Intensity = {};
init_grid_size = 51
i = init_grid_size;
for i=init_grid_size:50:1001 %grid size is varied here
  
    corr_grid = ones(i,i);
    Local_Density_tmp = imfilter(Cell_centroid_matrix,corr_grid); % Calculate local density by a translating box of size i*i, and counting the number of cells in the box. Linear spatial filtering with unity matrix is equivalent to this.
    Local_Intensity_tmp = imfilter(average_mean_image,corr_grid)/(i*i); % Calculate average intensity by translating a box of size i*i
    Local_Density = [Local_Density;Local_Density_tmp]; %add to local density matrix (array of arrays each corresponding to different size of box)
    Local_Intensity = [Local_Intensity;Local_Intensity_tmp]; % add to local intensity matrix (array of arrays each corresponding to different size of box)
    complete = (i/1001)*100;
    fprintf('| %3.0f%%',complete)
end