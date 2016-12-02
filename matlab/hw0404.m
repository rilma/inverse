%
% AUTHOR: RONALD R. ILMA
%
% Program to simulate a CAT inversion of an image
%
function hw0404(graph)

clear global; close all;

try if isempty(graph), graph = 1; end; catch graph = 1; end;

if isunix
    root_dir = '/media/sda5';
    myFontName = 'new century schoolbook';
else
    root_dir = 'D:';
    myFontName = 'NewCenturySchoolBook';
end

myFontSize = 10.0;
gpath = [root_dir '\Users\rilma\work\programs\matlab\cornell\eas5840\hw04\'];
fig_resolution = 300; orientation = 'landscape';

% Creting a image similar to suggested in problem
%
nrad = 32;
nang = 120;

nx = ceil(sqrt(nrad * nang));
ny = nx;

I = encode_image(nx, ny);
%

% Computing the sinogram
theta = 0 : 179;
R = radon(I, theta);

% Recovering the signal applying a backprojection without filtering
IR_nf = iradon(R, theta, 'linear', 'none');

% Recovering the signal applying a filtered backprojection
IR_f = iradon(R, theta, 'linear', 'Ram-Lak');

% Adding a uniform independent normally-distributed noise to 
% original image (zero-mean and 1% of peak value of image
u = 0.0;
omega = 0.1 * max(max(I));
noise = u + omega .* randn(nx, ny);

I_noise = I + noise;
%

% Computing the sinogram of noisy image, and recovering the original using
% backprojection (with or without filtering)
%
R_noise = radon(I_noise, theta, 64);

IR_noise_nf = iradon(R_noise, theta, 'linear', 'none');

IR_noise_f = iradon(R_noise, theta, 'linear', 'Ram-Lak');
%

% Graphical output
%
np_row = 2; np_col = 4;
hf = figure(2); clf(hf, 'reset');
fig_position = [10 40 275*np_col 180*np_row];
set(hf, 'Position', fig_position);

subplot(np_row, np_col, 1);
imagesc(I); title('Original');
colormap(gray); colorbar;

subplot(np_row, np_col, 2);
imagesc(R); title('Sinogram');
xlabel('\theta'); ylabel('x''');
colormap(gray); colorbar;

subplot(np_row, np_col, 3);
imagesc(IR_nf); title('Backprojection (without filtering)');
colormap(gray); colorbar;

subplot(np_row, np_col, 4);
imagesc(IR_f); title('Filtered backprojection');
colormap(gray); colorbar;

subplot(np_row, np_col, 5);
imagesc(I_noise); title('Original (with added noise)');
colormap(gray); colorbar;

subplot(np_row, np_col, 6);
imagesc(R_noise); title('Sinogram');
xlabel('\theta'); ylabel('x''');
colormap(gray); colorbar;

subplot(np_row, np_col, 7);
imagesc(IR_noise_nf); title('Backprojection (without filtering)');
colormap(gray); colorbar; 

subplot(np_row, np_col, 8);
imagesc(IR_noise_f); title('Filtered backprojection');
colormap(gray); colorbar;

% Saving figure
%
figfile = 'hw0404'; orientation = 'landscape';
save_figure(hf, graph, fig_resolution, orientation, [gpath figfile]);
%

function I = encode_image(nx, ny)
%
% Function to encode a image to be used in a CAT scan inversion procedure
%

    I = zeros(nx, ny);

    rx = nx / 2; ry = ny / 2;
    
    if rx >= ry; radii = rx; else radii = ry; end;

    for i = 1 : nx
        
        for j = 1 : ny
            
            x0 = i - rx; y0 = j - ry;
            
            r = sqrt(x0 ^ 2 + y0 ^ 2) / radii;
            
            if r >= 0.2 && r < 0.4; I(i, j) = 1; end;
            
            theta = (pi / 2 + atan(y0 / x0)) * 180 / pi;
            
            if x0 <= 0 && (theta >= 60 && theta < 120) && (r >= 0.4 && r <= 1)
                I(i, j) = 1;
            end

            if x0 >= 0
                theta = theta + 180;
                if ((theta >= 180 && theta < 240) || ...
                        (theta >= 300 && theta < 360)) && (r >= 0.4 && r <= 1)
                    I(i, j) = 1;
                end
            end
            
            if r > 1; I(i, j) = 1; end;
            
        end
    end
