%
% CORNELL UNIVERSITY            EAS 5840: Inverse Methods
% HOMEWORK: 5                   PROBLEM: 3
% AUTHOR: RONALD R. ILMA        DATE: Nov 07, 2008
%
% Application of the Conjugate Gradient Least-Squares (CGLS) technique for
% image deblurring
%
function hw0503(graph, nloop, inloop, nomega)

    clear global; close all;

    try if isempty(graph), graph = 1; end; catch graph = 1; end;
    try if isempty(nloop), nloop = 100; end; catch nloop = 100; end;
    try if isempty(inloop), inloop = 1; end; catch inloop = 1; end;    
    try if isempty(nomega), nomega = 0.01; end; catch nomega = 0.01; end;

    if isunix
        root_dir = '/home/rilma/';
        myFontName = 'new century schoolbook';
    else
        root_dir = 'D:';
        myFontName = 'NewCenturySchoolBook';
    end

    myFontSize = 10.0;
    gpath = [root_dir 'work/programs/matlab/cornell/eas5840/hw05/'];
    fig_resolution = 300; orientation = 'landscape';
    
    % Reading image
    img = double(imread([gpath 'idl_dog_bw_small.jpg']));
    
    % Finding dimensions of image
    [nx, ny] = size(img);
     
    % Adding a uniform independent normally-distributed noise to 
    % original image (zero-mean and 1% of peak value of image
    u = 0.0;
    omega = nomega * max(max(img));
    noise = u + omega .* randn(nx, ny);
    
    img_noisy = img + noise;
    
    % blurred theory
    G = blurred_theory(nx, 3.0, 1.5);
    
    % Generate blurred image
    dimg_blurred = G * reshape(img, nx * ny, 1) + ...
        (u + omega * randn(nx * ny, 1));
    img_blurred = reshape(dimg_blurred, nx, ny);
    
    % Use CGLS to deblurring of an image.
    %
    d0 = zeros(nx*ny, 1); % Initial guess
    [dimg_recovered, l2_error, l2_model] = ...
        cgls_method(G, dimg_blurred, nloop, d0);
    
    img_recovered = reshape(dimg_recovered(:,inloop), nx, ny);
 
    for j = 1 : nloop 
        fprintf('%i %f %f\n', [j l2_error(j) l2_model(j)]);
    end
    
% Graphical routines
%
    
    hf = figure(1);
    np_row = 2; np_col = 2;
    fig_position = [10 40 np_col*300 np_row*300];
    set(hf, 'Position', fig_position);
    
    colormap('gray');

    subplot(np_row, np_col, 1);
    imagesc(img); title('Original');
    set(gca, 'FontName', myFontName, 'FontSize', myFontSize);
        
    subplot(np_row, np_col, 2);
    imagesc(img_blurred); title('Noisy and Blurred');
    set(gca, 'FontName', myFontName, 'FontSize', myFontSize);    
    
    myXLim = [min(l2_error) max(l2_error)];
    subplot(np_row, np_col, 3);
    loglog(l2_error, l2_model, 'ko'); 
    title(['L-Curve (' num2str(nloop, '%i') ' iterations)']);
    xlabel('|| Gm - d ||'); ylabel('|| m ||');
    set(gca, 'FontName', myFontName, 'FontSize', myFontSize, ...
        'XLim', myXLim);
    
    subplot(np_row, np_col, 4);
    imagesc(img_recovered); 
    title(['Recovered (at ' num2str(inloop, '%i') ' iterations)']);
    set(gca, 'FontName', myFontName, 'FontSize', myFontSize);    
    
% Saving figure
%
    figfile = 'hw0503'; orientation = []; fig_resolution = [];
    save_figure(hf, graph, fig_resolution, orientation, [gpath figfile]);
%    
    
function output = blurred_theory(n, u, sigma)
%
% Calculate matrix theory used in digital imaging deblurring
%
    % point-spread function
    
    a = [exp(-((0:u-1).^2)/(2*sigma^2)), zeros(1, n - u)];
    
    g = toeplitz(a); g = sparse(g);
    
    g = ((2*pi*sigma.^2).^(-1)) * kron(g, g);
    
    output = g;
    
    
function [sol, l2_error, l2_model] = cgls_method(G, d, nint, d0) 
%
% Conjugate gradient Least-Squares (CGLS) algorithm. 
%
% The subroutine takes as argument:
%   G: The theory matrix
%   d: data vector
%   nint: number of iterations
%   d0: Initial guess
% 
% It returns:
%   sol: Final solution
%   l2_error: 2-norm prediction error
%   l2_model: 2-norm model complexity length
%

    [~, n] = size(G); sol = zeros(n, nint); 
    l2_model = zeros(nint, 1); l2_error = l2_model; 

% Default value of initial guess
    try if isempty(d0), d0 = zeros(n, 1); end; catch d0 = zeros(n, 1); end;

% Initializing variables to iterate. 
    D = (d' * G)';
    r = d; 
    normr2 = D' * D; 
 
% Loop. 
    for i = 1 : nint         
 
  % Update x and r vectors. 
      Gd = G * D; alpha = normr2 / (Gd' * Gd); 
      d0  = d0 + alpha * D; 
      r  = r - alpha * Gd; 
      s  = (r' * G)';
  
  % Update D vector. 
      normr2_new = s' * s; 
      beta = normr2_new / normr2; 
      normr2 = normr2_new; 
      D = s + beta * D; 
      sol(:, i) = d0; 
   
  % Compute the L2 prediction error.
      l2_error(i) = norm(r); 
  
  % Compute the L2 model complexity
      l2_model(i) = norm(d0);
  
    end
    