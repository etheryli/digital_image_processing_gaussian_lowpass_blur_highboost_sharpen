% Problem 3 : Frequency Filtering and Image Sharpening 

clear all;
close all;
clc;

% Define D0 for Gaussian Low-pass Filter
D0 = 50;

% Replicate or Symmetric (book says to use symmetric) padding
image = imread('blurry-moon.tif');


% MxN , 
[M, N] = size(image);
% Pad size (PxQ) is [2M 2N]
P = 2*M;
Q = 2*N;

% Find DFT of image with padding F(u,v) using 2-d 
% fast fourier transform (padding and center built-in)
F = fft2(double(image), P, Q);

% Generate H(u,v) filter function of size PxQ with centre (P/2,Q/2)

% Construct D(u,v) with a symmetric PxQ mesh grid centre (P/2,Q/2)
% and apply the Gaussian distribution

% Ranges from zero to size/2 and -size/2 to end for
% symmetric mesh grid 
x = 0:(P-1);
y = 0:(Q-1);

% For any greater than size/2, subtract size to get -size/2 to end
x(x>P/2) = x(x>P/2)-P;
y(y>Q/2) = y(y>Q/2)-Q;

[v, u] = meshgrid(y, x);

% Distance squared equation (D(u,v))^2 will get rid of negative values 
D_squared = (u.^2 + v.^2);

% Gaussian H(u,v) = e^(-(D(u,v)^2)/(2*D_0^2) applied D0 is specified above
H = exp(-(D_squared)/(2*D0^2));

% Form product G(u,v) = H(u,v)F(u,v)
G = H.*F;

% Find padded array output from inverse DFT (2d inverse fast 
% fourier transform) and cast to appropriate type of input image (uint8)
g = uint8(ifft2(G));

% Obtain processed output from M x N Region
gaussian_low_pass_filtered_output_image = g(1:M, 1:N);

% Obtain unsharp mask
unsharp_mask = image - gaussian_low_pass_filtered_output_image;

% Obtain final result sharpened image
output_image = unsharp_mask + image;

% Obtain high-boost filtering with k = 5;
k = 5;
highboost_filtered_output_image = unsharp_mask * k + image;

% Plotting
fig1 = figure('Name', 'Frequency Filtering and Image Sharpening', 'color', [1 1 1]);
subplot(3,3,1);
imshow(image);
title("Input Image");

subplot(3,3,4);
imshow(gaussian_low_pass_filtered_output_image);
title("Blurred Image");

subplot(3,3,2);
imshow(output_image);
title("Sharpened Output Image");

subplot(3,3,5);
imshow(unsharp_mask);
title("Unsharp Mask");

subplot(3,3,3);
imshow(highboost_filtered_output_image);
title("High-Boost Filtered Sharpened Output Image (k = 5)");

fig2 = figure('Name', 'Spectra of Gaussian Low-Pass Filtering Process', 'color', [1 1 1]);
subplot(1,3,1);
F_spectrum = log(1+abs(fftshift(F)));
imshow(F_spectrum, []);
title("F(u,v)");
colormap([0,0,0]);
subplot(1,3,2);
H_spectrum = log(1+abs(fftshift(H)));
imshow(H_spectrum, []);
title("H(u,v)");

subplot(1,3,3);
G_spectrum = log(1+abs(fftshift(G)));
imshow(G_spectrum, []);
title("G(u,v) = H(u,v) F(u,v)");