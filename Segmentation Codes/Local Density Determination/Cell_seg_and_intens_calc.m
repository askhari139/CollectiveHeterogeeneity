% Use only bright field images. For phase contrast images, invert
% the binary mask and proceed with segmentation.
%Select directory - Put all bright field images and tmrm in separate directories. Select
%the bright field directory initially, followed by TMRM intensity
%directory.
BF_Dir = uigetdir;
BF_Files = dir(fullfile(BF_Dir,'*.tif'));
mkdir(BF_Dir,'intensity');
mkdir(BF_Dir,'overlay');
TM_Dir = uigetdir;
TM_Files = dir(fullfile(TM_Dir,'*.tif'));
%%
for x = 1:length(BF_Files)
    % Read file and adjust contrast
    FileName = BF_Files(x).name;
    fullFileName = fullfile(BF_Dir, FileName);
    fprintf(1, 'Now reading %s\n', FileName);
    I = imread(fullFileName); % reading file
    J = imadjust(I); % contrast adjusting the file
    %imshow(J);
    % Splitting regions into required number of regions
    split = ones(1,8);
    split = split*256;
    Cell = mat2cell(I,split,split);
    rows = numel(split);
    columns = numel(split);
    % Segmenting each region individually:
    for i = 1:rows
        for j = 1:columns
            R = Cell{i, j};
            %imshow(R);
            disp(i);
            % Applying gaussian smoothening to every region
            M = imadjust(R);
            M = imgaussfilt(M,3);
            %M = imguidedfilter(M,'NeighborhoodSize',6);
            %M = wiener2(M,[8 8]);
            %imshow(M);
            % Active contour based edge detection
            level = graythresh(M);
            mask = im2bw(M,level); % For phase contrast images, invert this binary mask and save it as mask
            BW = activecontour(M, mask, 4, 'Chan-Vese');
            maskedImage = M;
            maskedImage(~BW) = 0;
            % Display mask
            %imshow(BW);
            % Filled binary mask
            BWd = imfill(BW, 'holes');
            %imshow(BWd);
            % Dilated binary mask
            se90 = strel('line', 1, 90);
            se0 = strel('line', 1, 0);
            maskd = imdilate(BWd, [se90 se0]);
            %imshow(maskd);
            % Interior noise removed binary mask. Adjust the area to get desired
            % results.
            maska = bwareafilt(maskd,[100 Inf]);
            sero = strel('diamond', 1);
            maska = imerode(maska,se90);
            sedil = strel('diamond', 1);
            maska = imdilate(maska, sedil);
            %imshow(maska);
            
            cmas = imcomplement(maska);
            cmas = bwareafilt(cmas,[1 500]);
            %imshow(cmas);
            maska = maska +cmas;
            %imshow(maska);
            maskb = logical(maska);
            %imshow(maskb);
            % Small areas filtered binary mask
            mask7 = bwareafilt(maskb, [100 1499]);
            %imshow(mask7);
            % Large area filtered binary mask
            mask2 = bwareafilt(maskb, [1500 Inf]);
            %imshow(mask2);
            mask3 = imdilate(mask2, sero);
            %imshow(mask3);
            % Watershed transformed binary mask
            mask4 = mask3;
            mask5 = mask4;
            D = -bwdist(~mask3);
            L = watershed(D);
            mask4(L==0)=0;
            mask6 = imextendedmin(D,4); % Change the numerical argument for the function to get optimum results
            D2= imimposemin(D,mask6);
            Ld2=watershed(D2);
            mask5(Ld2==0) = 0;
            %imshow(mask5);
            
            cmas5 = imcomplement(mask5);
            cmas5 = bwareafilt(cmas5,[1 200]);
            %imshow(cmas5);
            mask5 = mask5 +cmas5;
            %imshow(mask5);
            mask1 = logical(mask5);
            %imshow(mask1);
            % Overlay final mask with input image
            mask8 = mask1+mask7;
            Out_cell{i, j} = mask8;
        end
    end
    
    final_mask_2 = cell2mat(Out_cell);
    %imshow(final_mask_2);
    
    %imwrite(final_mask_2,'mask_256.tif'); %writing mask to file
    
    %
    %X = imread('mask_512.tif');
    %imshow(X)
    %
    %Y = imread('mask_256.tif');
    %imshow(Y);
    %
    %Z = imread('mask_cust.tif');
    %G = X+Y+Z;
    %
    %imwrite(G,'trial_1.tif');
    %imshow(G);
    %
    Y = final_mask_2;
    Tm_file_name = TM_Files(x).name;
    fullname = fullfile(TM_Dir,Tm_file_name);
    tmrm = imread(fullname); % reading tmrm intensity file
    tmrm_dis = imadjust(tmrm); % contrast adjustment of tmrm intensity file
    
    %Skeletonize
    Inv = imcomplement(Y); %inverting binary mask
    %imshow(Inv);
    Skl = bwmorph(Inv,'skel',Inf); %skeletonize binary mask
    Spr = bwmorph(Skl, 'spur',Inf); %remove free edges
    Spr = imcomplement(Spr); %inverting the spur removed skeleton mask
    %imshow(Spr);
    % Calculating region properties
    cc = bwconncomp(Spr,4); % identifying connected components
    Areas = regionprops(cc,'Area','Perimeter','Centroid'); %finding area of every connected component
    Spr_u16 = uint16(Spr); % converting logical skeleton into unsigned integer array
    display = tmrm_dis.*Spr_u16; % Applying mask on contrast adjusted tmrm intensity image for display
    display_1 = J.*Spr_u16;
    multiply = tmrm.*Spr_u16; %Applying mask on original tmrm intensity image
    %imshow(display_1); % Displaying the mask overlayed tmrm intensity image
    pth = '\overlay\';
    in_fil_pth = strcat(BF_Dir,pth,FileName);
    in_fil_pth = strrep(in_fil_pth,'.tif','');
    ext = '.tif';
    in_fil_nme = strcat(in_fil_pth,ext);
    %print(in_fil_pth\FileName,'-png','-dpng','-r200');
    %saveas(gcf,in_fil_pth,'tif');
    imwrite(display_1,in_fil_nme,'tif')
    region = regionprops(cc,multiply,'MeanIntensity'); % calculating average intensity in every cell
    
    
    %% Displaying image with mean intensity
    Mean_image = zeros(2048,2048);
    Cells = cc.NumObjects;
    Pixel_loc = regionprops(cc,'PixelList');
    %area = ones(1,Cells);
    %intensity = ones(1,Cells);
    for i = 1:Cells
        area = Areas(i).Area;
        intensity = region(i).MeanIntensity;
        for j = 1:area
            Location = [Pixel_loc(i).PixelList(j,1) Pixel_loc(i).PixelList(j,2)];
            Mean_image(Location(2),Location(1)) = intensity;
        end
    end
    imagesc(Mean_image);
   
    colorbar;
    colormap jet;
    pth = '\intensity\';
    in_fil_pth = strcat(BF_Dir,pth,FileName);
    in_fil_pth = strrep(in_fil_pth,'.tif','');
    ext2 = '.svg';
    in_fil_pth = strcat(in_fil_pth,ext2);
    saveas(gcf,in_fil_pth,'svg');
     hold on
    figure();
end %Run the code only till this point for routine purposes.
%%
intensity = ones(1,Cells);
for i = 1:Cells
    intensity(1,i) = region(i).MeanIntensity;
end
histogram(intensity);
area = ones(1,Cells);
for i = 1:Cells
    area(1,i) = Areas(i).Area;
end
histogram(area);
hold off