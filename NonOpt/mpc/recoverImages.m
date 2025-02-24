function recoverImages

% Copyright (C) 2025 Frank E. Curtis
%
% This code is published under the MIT License.
%
% Author(s) : Frank E. Curtis, Lara Zebiane

% Loop over regularizers
for regularizer = 0:3

  % Read as matrix to text file
  I = readmatrix(sprintf("croissant_recovered_%d.txt",regularizer));
  I = uint8(I);

  % Show recovered image
  figure(3+regularizer);
  imshow(I);

  % Write noisy image
  imwrite(I,sprintf("croissant_recovered_%d.png",regularizer));

end