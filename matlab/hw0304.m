%
% AUTHOR: RONALD R. ILMA
%
% Program to investigate the vertical seismic profiling (VSP)
% of the bore hole (based on Homework 1 and 2, problems 4 and 2, 
% respectively)
%
%
function hw0304(epsilon, graph, m, n, p)

clear global; close all;

% Damping factor
try if isempty(epsilon), epsilon = 5e-4; end; catch epsilon = 5.e-4; end;

try if isempty(graph), graph = 1; end; catch graph = 1; end;

% # of elements in "Model vector (# of columns in "Thoery matrix")
try if isempty(m), m = 400; end; catch m = 400; end;

% # of rows in "Theory matrix"  or elements in "Data vector"
        % (or # of samples in depth)
try if isempty(n), n = 100; end; catch n = 100; end;

gpath = 'D:\Users\rilma\work\programs\matlab\cornell\eas5840\hw03\';
fig_resolution = 100; orientation = 'landscape';

% Input parameters
%
z0 = 0.0;  % Initial depth [m]
zf = 40.0; % Maximum depth [m]

s0 = 2.0; % Initial Slowness [s/km]
sf = 5.0; % Final Slowness [s/Km]
        
mu = 0.0; % Mean of normal dist.        
sigma = 1e-4; % [s] standard deviation of normally distri. noise

%

% conversion factor from m to km
m2km = 1e-3;

% sampled data in depth [km]
z = (0 : m - 1) ./ (m - 1) .* (zf - z0) .* m2km;
deltaZm = z(2) - z(1); 

Zn = (0 : n - 1) ./ (n - 1) .* (zf - z0) .* m2km;
deltaZn = Zn(2) - Zn(1); 

% Initialize G "operator matrix"
%

G = zeros(n, m);

for i = 1 : n
    for j = 1 : m
        if z(j) <= Zn(i)
            G(i, j) = deltaZm;
        elseif j > 1 && z(j - 1) < Zn(i)
            G(i, j) = Zn(i) - z(j - 1);
        else
        end
    end
end


% Initializes model vector
%

% Getting slowness values from initial model
s = slowness_funct(z, [z0*m2km s0 zf*m2km sf])';

% Imposing an depth interval with constant slowness
%
zrange_ipwf = [3/8 5/8] .* (zf - z0) .* m2km;
ind = find(z >= zrange_ipwf(1) & z < zrange_ipwf(2));
if ~isempty(ind), s(ind) = s0; end;
%

% Applying the forward approach
%
% Travel-Time or delay [s]
T = G * s;
%

% Generating a normally distributed noise
N = mu + sigma .* randn(n,1);

% Adding noise to data
TN = T + N;

% Implementation of SVD approach
%

[m_best, sv] = calc_svd_model(G, T, []);
m_noise = calc_svd_model(G, TN, []);

try if isempty(p), p = []; end; catch p = []; end;
m_noise1 = calc_svd_model(G, TN, p);
m_noise2 = calc_svd_model(G, TN, 10);
%

% Compute data vector from 
TN4_best = G * m_best;
TN4_noise = G * m_noise;
TN4_noise1 = G * m_noise1;
TN4_noise2 = G * m_noise2;

% 1st derivative of travel time (slowness)
ss4_best = diff(TN4_best) ./ deltaZn;
ss4_noise = diff(TN4_noise) ./ deltaZn;
ss4_noise1 = diff(TN4_noise1) ./ deltaZn;
ss4_noise2 = diff(TN4_noise2) ./ deltaZn;
%

% Graphical output commands
%
if graph > 0

    if isunix, myFontName = 'new century schoolbook'; ...
            else myFontName = 'NewCenturySchoolBook'; end    
    
    figh = figure(1); clf(figh, 'reset');    
    
    fig_position = [10 40 1000 450];
    set(figh, 'Position', fig_position);
    
    np_row = 1; np_col = 2;
        
    subplot(np_row,np_col, 1);
    plot(sv, 'Color','k', 'LineWidth',2, 'Marker','o');
    set(gca, 'FontName', myFontName, ...
        'YMinorTick','on');
    title('Singular Values of G');
    xlabel('Index'); ylabel('Value');
    grid on;
        
    subplot(np_row,np_col,2)

    line(Zn(1:(n-1)), ss4_best, 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1, ...
        'Marker', 'o','MarkerSize',2);
    line(Zn(1:(n-1)), ss4_noise, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
    line(Zn(1:(n-1)), ss4_noise1, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(Zn(1:(n-1)), ss4_noise2, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);    
    set(gca, 'FontName', myFontName, 'XLim', [z0 zf].*m2km, ...
        'XMinorTick', 'on', ...
        'YLim', [1.5 5.5], 'YMinorTick', 'on');                    
    title('Slowness vs Depth');
    xlabel('[km]'); ylabel('[s/km]');
    legend('Best','Noisy data (p = 100)',['Truncated (p = ' num2str(p,'%2i') ')'], ...
        'Truncated (p = 10)', 'Location', 'NorthWest');
    grid on;    
    
    % Saving figure
    %
    figfile = ['hw0304_' num2str(n,'%03i') '-' num2str(m,'%03i') '-' ...
        num2str(sigma,'%6.4f')];
    save_figure(figh, graph, fig_resolution, orientation, [gpath figfile]);    
    %
end
%

