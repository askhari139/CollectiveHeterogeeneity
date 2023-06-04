Int_prcntile_cutoffs = [];
for i = 1:20
    crop_beg_ind = (i*(init_grid_size-1))/2;
    crop_end_ind = length(Mean_image)-crop_beg_ind;
    cropped_array_tmp = cell2mat(cropped_Intensit_arrays(i));
    cut_off = prctile(cropped_array_tmp(crop_beg_ind+1:crop_end_ind,crop_beg_ind+1:crop_end_ind),[5 10 20 95 90 80],'all');
    Int_prcntile_cutoffs = [Int_prcntile_cutoffs;cut_off'];
end 
 Int_cutoff_table = array2table(Int_prcntile_cutoffs,'VariableNames',{'5%','10%','20%','95%','90%','80%'});
%%
Int_top_bottom_thresh_cell = {};
for i=1:20
    i
    Intensity_value_matrix = cropped_Intensit_arrays(i);
    for j=1:3
        j
    Intensity_matrix_tmp = cell2mat(Intensity_value_matrix);
    Intensity_matrix_tmp(Intensity_matrix_tmp>Int_prcntile_cutoffs(i,j)) = 0;
    Intensity_matrix_tmp(Intensity_matrix_tmp>0) = 1;
    Int_top_bottom_thresh_cell{i,j} = Intensity_matrix_tmp;
    end
    
    for j=4:6
    Intensity_matrix_tmp = cell2mat(Intensity_value_matrix);
    Intensity_matrix_tmp(Intensity_matrix_tmp<Int_prcntile_cutoffs(i,j)) = 0;
    Intensity_matrix_tmp(Intensity_matrix_tmp>0) = 1;
    Int_top_bottom_thresh_cell{i,j} = Intensity_matrix_tmp;
    end
end
%%
top_bottom_density_cell = {};
for i = 1:20
    for j=1:6
    dens_top_bottom_tmp = cell2mat(Local_Density(i)).*cell2mat(Int_top_bottom_thresh_cell(i,j));
    top_bottom_density_cell{i,j} = dens_top_bottom_tmp;
    end
end
%%
top_bottom_dens_sorted = {}
for i = 1:20
    for j = 1:6
        tmp_value1 = cell2mat(top_bottom_density_cell(i,j));
        top_bottom_dens_sorted{i,j} = tmp_value1(tmp_value1~=0);
    end
end 
%%
mean_top_bottom_dens = []
for i = 1:20
    for j = 1: 6
        mean_tb_den = mean(cell2mat(top_bottom_dens_sorted(i,j)));
        std_tb_den = std(cell2mat(top_bottom_dens_sorted(i,j)));
        mean_top_bottom_dens(i,j) = mean_tb_den;
    end
end
%%
for i = 1:20
    plot(mean_top_bottom_dens(i,:))
    hold on
end