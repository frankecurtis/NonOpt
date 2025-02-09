function J = addNoiseToImage

% Load image
I = imread("croissant.png");

% Resize image
I = imresize(I,[120 160]);

% Show original image
figure(1);
imshow(I);

% Add noise
J = imnoise(I,'gaussian',0.0,0.005);

% Show noisy image
figure(2);
imshow(J);

% Write noisy image
imwrite(J,"croissant_noisy.png");

% Write as matrix to text file
writematrix(J,'croissant_matrix.txt','Delimiter','tab');