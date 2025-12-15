function addNoiseToImage

% Copyright (C) 2025 Frank E. Curtis
%
% This code is published under the MIT License.
%
% Author(s) : Frank E. Curtis, Lara Zebiane

% Load color image
I = imread("croissant_color.png");

% Convert to grayscale
I = rgb2gray(I);

% Resize grayscale image
I = imresize(I,0.2);

% Show original grayscale image
figure(1);
imshow(I);

% Write original grayscale image
imwrite(I,"croissant.png");

% Write as matrix to text file
writematrix(I,'croissant_original_matrix.txt','Delimiter','tab');

% Set random number generator seed
rng(19970830);

% Add noise
%J = imnoise(I,'gaussian',0.0,0.02);
J = imnoise(I,'salt & pepper');

% Show noisy image
figure(2);
imshow(J);

% Write noisy image
imwrite(J,"croissant_noisy.png");

% Write as matrix to text file
writematrix(J,'croissant_matrix.txt','Delimiter','tab');