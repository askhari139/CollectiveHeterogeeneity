% %figure()
corr_size = round(((2048/sqrt(Cells))*7)/50); %size of the grid with around 50 cells (7*7)
% %X = categorical({'Low Density','High Density',});
% %X = reordercats(X,{'Low Density','High Density'});
% %Y = [mean_top_bottom(corr_size,1) mean_top_bottom(corr_size,4) ;mean_top_bottom(corr_size,2)  mean_top_bottom(corr_size,5) ;mean_top_bottom(corr_size,3) mean_top_bottom(corr_size,6) ];
% %std = [std_top_bottom(corr_size,1) std_top_bottom(corr_size,4); std_top_bottom(corr_size,2)  std_top_bottom(corr_size,5) ;std_top_bottom(corr_size,3) std_top_bottom(corr_size,6) ];
% %figure()
% %b = bar(Y, 'grouped');
% %hold on
% % Calculate the number of bars in each group
% nbars = size(Y, 2);
% % Get the x coordinate of the bars
% x = [];
% for i = 1:nbars
%     x = [x ; b(i).XEndPoints];
% end
% % Plot the errorbars
% errorbar(x',Y,std,'k','linestyle','none')'
% hold off
%%
for i=1:6
   filename_csv = [FileName(1:end-4) 'top_bottom_' num2str(i) '.csv'] ;
   csvwrite(filename_csv,cell2mat(top_bottom_int_sorted(corr_size,i)));
end
%%
filename_csv1 = [FileName(1:end-4) '_mean_top_bottom.csv'];
filename_csv2 = [FileName(1:end-4) '_std_top_bottom.csv'];
csvwrite(filename_csv1,mean_top_bottom);
csvwrite(filename_csv2,std_top_bottom);
%%
hold off
xvalues = {'low 10%','low 20%','low 25%','top 10%','top 20%','top 25%'};
yvalues = {};
for i = 1:20
    legend_name = int2str(i*50);
    legend_name = strcat(legend_name, 'px');
    yvalues = [yvalues legend_name]
end
yvalues(1) = {'1 cell'};
%%
figure()
h = heatmap(xvalues,yvalues,mean_top_bottom(:,1:6));
h.ColorScaling = 'scaledcolumns';
h.Colormap = pink;
h.XLabel = 'Circularity';
h.YLabel = 'Grid Size';
h.Title = 'Normalized TMRM Intensity';
filename_hm_sc = [FileName(1:end-4) 'hm_scaled_col.png'];
saveas(gcf,filename_hm_sc)
filename_hm_sc = [FileName(1:end-4) 'hm_scaled_col.fig'];
savefig(gcf,filename_hm_sc)
%%
figure()
h = heatmap(xvalues,yvalues,mean_top_bottom(:,1:6));
%h.ColorScaling = 'scaledrows';
h.Colormap = pink;
h.XLabel = 'Circularity';
h.YLabel = 'Grid Size';
h.Title = 'Normalized TMRM Intensity';
filename_hm = [FileName(1:end-4) '_hm.png'];
saveas(h,filename_hm);