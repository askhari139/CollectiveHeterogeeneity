region3 = regionprops(cc,multiply,'MeanIntensity','Centroid');
at = cat(1,region3.Centroid);
Im1 = zeros(2048,2048);
at = floor(at);
i = 1;

while i<length(at)
   
    Im1(at(i,1),at(i,2)) = 1;
    i = i+1;
end

Im1 = Im1';