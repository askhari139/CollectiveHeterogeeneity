% HETEROGENEITY INDEX CALCULATION PROGRAM ------ v3.0; JUNE 2023 ----------
% Author: Basil Thurakkal, TIFR Hyderabad --- Contact: basilt@tifrh.res.in 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

TM_Dir = uigetdir;
TM_Files = dir(fullfile(TM_Dir,'*.tif'));
mean_het = {};

Mask_dir = uigetdir;
Mask_files = dir(fullfile(Mask_dir,'*.tif'));

for x = 1:length(TM_Files)
% Read the image
    FileName = TM_Files(x).name;
    fullFileName_TM = fullfile(TM_Dir, FileName);
    fullFileName_mask = fullfile(Mask_dir, FileName);
    FileName

image_int  = imread(fullFileName_TM); %TMRM image
mask_img = imread(fullFileName_mask); %mask image
image_org = image_int;
bin_img = im2bw(mask_img);
bin_img = 1- bin_img;

cc = bwconncomp(bin_img,4)% connected component
reg_prop = regionprops(cc,image_int,'MeanIntensity','PixelList','Area')% regionprop (Mean intensity )
int_array = [reg_prop(:).MeanIntensity]; %struct field to array
image_int = double(image_int)./(mean(int_array));
%% Define the initial grid size and step size
initialGridSize = 256;
gridSizeStep = 256;

% Get the dimensions of the image
[height, width, ~] = size(image_int);
het_index_cell = {};
het_var_cell = {};
gridSizeCount = 1; %
%% Loop over the grid sizes

for gridSize = initialGridSize :gridSizeStep: min(height, width)
    % Calculate the number of rows and columns in the grid    
    numRows = floor(height / gridSize);
    numCols = floor(width / gridSize);

    het_index_cell{gridSizeCount,1} = gridSize;
     het_var_cell{gridSizeCount,1} = gridSize;

    % Loop over each grid cell

    gridCount = 2;

    if gridSize <= height/2 | gridSize == height
    for row = 1 : numRows       
        for col = 1 : numCols
            % Calculate the starting and ending indices for the current grid cell
            startRow = (row - 1) * gridSize + 1;
            endRow = startRow + gridSize - 1;
            startCol = (col - 1) * gridSize + 1;
            endCol = startCol + gridSize - 1;

            % Extract the current grid cell for mask and intensity image
            gridCell_int = image_int(startRow:endRow, startCol:endCol, :);
            gridCell_mask = bin_img(startRow:endRow, startCol:endCol, :);
            % Display the current grid cell
            subplot(numRows, numCols, (row - 1) * numCols + col);
            %imagesc(gridCell_int); axis off;
            
            [het_i,var_het,reg_prop] = het_index(gridCell_int,gridCell_mask); % *HETEROGENEITY INDEX FUNCTION*
            
            het_index_cell{gridSizeCount,gridCount} = het_i;
            het_var_cell{gridSizeCount,gridCount} = var_het;
            gridCount = gridCount + 1;

        end
    
    end

    elseif gridSize > height/2 && gridSize ~= height
    for row = 1 : numRows       
        for col = 1 : numCols
            % Calculate the starting and ending indices for the current grid cell
            startRow = (row - 1) * gridSize + 1;
            endRow = startRow + gridSize - 1;
            startCol = (col - 1) * gridSize + 1;
            endCol = startCol + gridSize - 1;
            
            image_temp = image_int;
            bin_temp = bin_img;
            for i = 1:4             
            % Extract the current grid cell for mask and intensity image
            gridCell_int = image_temp(startRow:endRow, startCol:endCol, :);
            gridCell_mask = bin_temp(startRow:endRow, startCol:endCol, :);
            % Display the current grid cell
            subplot(numRows, numCols, (row - 1) * numCols + col);
            %imagesc(gridCell_int); axis off;
            
            [het_i,var_het,reg_prop] = het_index(gridCell_int,gridCell_mask); % *HETEROGENEITY INDEX FUNCTION*
            het_index_cell{gridSizeCount,gridCount} = het_i;
            het_var_cell{gridSizeCount,gridCount} = var_het;

            gridCount = gridCount + 1;
            image_temp = rot90(image_temp);
            bin_temp = rot90(bin_temp);
            display(i)
            end

        end
    
    end    
    end
    
    % Pause to display each grid size
    gridSizeCount = gridSizeCount+1;
    pause(0.25);

end
%%
mean_het{1,x} = FileName;
SEM_het{1,x} = FileName; 
for i = 1:gridSizeCount-1
    M = [het_index_cell{i,2:end}]; % heterogeneity index values
    var = [het_var_cell{i,2:end}];
    mean_het_i = mean(M);
    sd_het_i = sqrt(sum(var)/length(var)); % SD_of_mean = Sqrt(sumsqr(var)/n)
    mean_het{i+1,x} = mean_het_i;
    SEM_het{i+1,x} = sd_het_i/sqrt(length(var)); %SEM = SD/(sqrt(n))
end
end

writecell(SEM_het,'SEM_het_results.csv');
writecell(mean_het,'results_het.csv');
