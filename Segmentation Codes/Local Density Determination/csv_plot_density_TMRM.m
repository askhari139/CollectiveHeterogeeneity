BF_Dir = uigetdir;
BF_Files_mean = dir(fullfile(BF_Dir,'*mean_top_bottom.csv'));
mean_top_bottom_cell = {}
BF_Files_std = dir(fullfile(BF_Dir,'*std_top_bottom.csv'));

for x = 1:length(BF_Files_mean)
   FileName = BF_Files_mean(x).name;
   fullFileName = fullfile(BF_Dir, FileName);
   I = csvread(fullFileName);
   mean_top_bottom_cell{x,1} = FileName;
   mean_top_bottom_cell{x,2} = I
end

std_top_bottom_cell = {}
for x = 1:length(BF_Files_std)
   FileName = BF_Files_std(x).name;
   fullFileName = fullfile(BF_Dir, FileName);
   I = csvread(fullFileName);
   std_top_bottom_cell{x,1} = FileName;
   std_top_bottom_cell{x,2} = I 
end
%%
mean_tb_seven = [];
std_tb_seven = [];
for i=1:length(mean_top_bottom_cell)
celldia_pix = 2048/(sqrt(mean_top_bottom_cell{i,2}(1,7)));
box_size = round((7*celldia_pix)/50)%box width size correspoinding to average 7 cell 
mean_tb_seven = [mean_tb_seven; mean_top_bottom_cell{i,2}(box_size,1:6)];
std_tb_seven = [std_tb_seven; std_top_bottom_cell{i,2}(box_size,1:6)];
end



