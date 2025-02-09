function I = recoverImage

% Read as matrix to text file
I = readmatrix("croissant_recovered_matrix.txt");
I = uint8(I);

% Show recovered image
figure(3);
imshow(I);

% Write noisy image
imwrite(I,"croissant_recovered.png");