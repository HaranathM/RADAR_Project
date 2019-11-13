
%% 2-D Transform%

% This example is on Mathworks website too. Please refer that document
% for further details.

% The 2-D Fourier transform is useful for processing 2-D signals and other 2-D data such as images.
% Create and plot 2-D data with repeated blocks.

P = peaks(20);
X = repmat(P,[5 10]);
figure(1);
imagesc(X);

%% Computing the 2-D Fourier transform of the data.  
fft_X = fft2(X);

%% Shift the zero-frequency component to the center of the output
fft_X_shift = fftshift(fft_X);
fft_X_shift = abs(fft_X_shift);  % Complex values are not supported in "imagesc".
imagesc(fft_X_shift);
% The below code gives the image without shifting!!
% figure(2)
% imagesc(abs(fft_X))

