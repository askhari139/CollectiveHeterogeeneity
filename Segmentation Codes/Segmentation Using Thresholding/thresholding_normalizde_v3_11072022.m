%update from v2 --
% saturated values are imputed with 0. Also the image intensity values are
% normalized  with (max-min). Threshold is 25%

tot_intensity_image = Mean_image;
Max_int = max(tot_intensity_image,[],'all');
img_median = median(tot_intensity_image,'all');

cell_prop_cell = struct2cell(cell_prop);
cell_prop_mat = cell2mat(cell_prop_cell);

img_prc = prctile(cell_prop_mat(2,:),99);

cell_prop_mat(3,:) = cell_prop_mat(2,:)<img_prc;
cell_prop_mat(4,:) = cell_prop_mat(3,:).*cell_prop_mat(2,:);
mean_val = sum(cell_prop_mat(4,:))/sum(cell_prop_mat(3,:));
cell_prop_mat(5,:) = (cell_prop_mat(4,:)-min(cell_prop_mat(2,:)))/(max(cell_prop_mat(4,:))-min(cell_prop_mat(2,:)));

imputed_int = cell_prop_mat(2,:);
imputed_int(imputed_int>img_prc) = mean_val;
imputed_int(2,:) = (imputed_int(1,:)-min(imputed_int(1,:)))/(max(imputed_int(1,:))-min(imputed_int(1,:)));

thresh_img = tot_intensity_image;
thresh_img(thresh_img>img_prc) = mean_val;
thresh_img = (thresh_img-min(cell_prop_mat(2,:)))/(max(cell_prop_mat(4,:))-min(cell_prop_mat(2,:))); %normalization
img_copy = thresh_img;
thresh_img(thresh_img<(0.25)) = 0; %threshold here

%% histogram of cellwise intensities


%hist(cell_prop_mat(5,:),50);
hist(imputed_int(2,:),50);
save_name = strrep(int_file,'.tif','_hist.png')
saveas(gcf,save_name);
%%
thresh_img(thresh_img>0) = 255;
bwimg = im2bw(thresh_img,0.5);
% imagesc(thresh_img)
% figure()
% imagesc(bwimg)

%%
figure()
se = strel('square',4);
er_bwimg = imerode(bwimg,se);
imagesc(er_bwimg);
set(gca,'XTick',[],'YTick',[]);
colormap summer;
%%
bw_cell = im2bw(imcomplement(bin_img),0.5);
bw_cell = imdilate(bw_cell,se);
doub_bwimg = double(bwimg);
doub_bwcell = double(bw_cell);
a = 2.*doub_bwcell + doub_bwimg;
c_map = [1 1 0.4;0 0 0 ;0 0.5 0.4 ];
rgb_img = label2rgb(a,c_map,[0 0.5 0.4]);
imshow(rgb_img);

save_name = strrep(int_file,'.tif','_hilo.png')
saveas(gcf,save_name);