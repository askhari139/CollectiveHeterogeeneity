init_grid_size = 51;
pad_arrays = {};
cropped_density_arrays ={};
for i = 1:20
    pad_size = (i*(init_grid_size-1))/2;
    crop_size = length(Mean_image)-2*pad_size;
    pad_array_tmp = ones(crop_size,crop_size);
    pad_array_tmp = padarray(pad_array_tmp,[pad_size,pad_size]);
    pad_arrays = [pad_arrays;pad_array_tmp]; %for cropping out the boundaries of local density fields
    i
end
for i = 1:20
    cropped_array_tmp = cell2mat(Local_Density(i)).*cell2mat(pad_arrays(i));
    cropped_density_arrays = [cropped_density_arrays;cropped_array_tmp]; %density arrays cropped to adjust for the boundary problem
end
%%
cropped_Intensit_arrays = {}; %density arrays cropped to adjust for the boundary problem
for i = 1:20
    cropped_array_tmp = cell2mat(Local_Intensity(i)).*cell2mat(pad_arrays(i));
    cropped_Intensit_arrays = [cropped_Intensit_arrays;cropped_array_tmp];
end

%%
prcntile_cutoffs = [];
for i = 1:20
    crop_beg_ind = (i*(init_grid_size-1))/2;
    crop_end_ind = length(Mean_image)-crop_beg_ind;
    cropped_array_tmp = cell2mat(cropped_density_arrays(i));
    cut_off = prctile(cropped_array_tmp(crop_beg_ind+1:crop_end_ind,crop_beg_ind+1:crop_end_ind),[10 20 25 90 80 75],'all');
    prcntile_cutoffs = [prcntile_cutoffs;cut_off'];
end 
 cutoff_table = array2table(prcntile_cutoffs,'VariableNames',{'10%','20%','25%','90%','80%','75%'}); %top-bottom thresholds
%%
top_bottom_thresh_cell = {};
for i=1:20
    i
    Density_matrix = cropped_density_arrays(i);
    for j=1:3
        j
    Density_matrix_tmp = cell2mat(Density_matrix);
    Density_matrix_tmp(Density_matrix_tmp>prcntile_cutoffs(i,j)) = 0;
    Density_matrix_tmp(Density_matrix_tmp>0) = 1;
    top_bottom_thresh_cell{i,j} = Density_matrix_tmp;
    end
    
    for j=4:6
    Density_matrix_tmp = cell2mat(Density_matrix);
    Density_matrix_tmp(Density_matrix_tmp<prcntile_cutoffs(i,j)) = 0;
    Density_matrix_tmp(Density_matrix_tmp>0) = 1;
    top_bottom_thresh_cell{i,j} = Density_matrix_tmp;
    end
end
%%
top_bottom_int_cell = {};
for i = 1:20
    for j=1:6
    int_top_bottom_tmp = cell2mat(Local_Intensity(i)).*cell2mat(top_bottom_thresh_cell(i,j));
    top_bottom_int_cell{i,j} = int_top_bottom_tmp;
    end
end
%%
top_bottom_int_sorted = {}
for i = 1:20
    for j = 1:6
        tmp_value = cell2mat(top_bottom_int_cell(i,j));
        top_bottom_int_sorted{i,j} = tmp_value(tmp_value~=0);
    end
end 
%%
mean_top_bottom = [];
 std_top_bottom = [];
for i = 1:20
    for j = 1: 6
        mean_tb = mean(cell2mat(top_bottom_int_sorted(i,j)));
        std_tb = std(cell2mat(top_bottom_int_sorted(i,j)));
        mean_top_bottom(i,j) = mean_tb;
        std_top_bottom(i,j) = std_tb;
    end
    mean_top_bottom(i,7) = Cells;
end
%%
%figure()
%for i = 3:8
 %   plot(mean_top_bottom(i,:))
  %  hold on
%end
%hold off