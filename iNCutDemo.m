load temp;
imshow(uint8(x)); %¿é¤J¹Ï¤ù
figure;imshow(uint8((y+1)*255/2));  %initial segmentation
z = iNCut1(x,y,[20 100 10 10000]);  %r sigmaX sigmaI iter_max
figure;imshow(uint8((z+1)*255/2));