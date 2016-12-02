%
% CORNELL UNIVERSITY            EAS 5840: Inverse Methods
% HOMEWORK: 2                   PROBLEM: 3
% AUTHOR: RONALD R. ILMA        DATE: Sept 26, 2008
%
% Program to investigate the vertical seismic profiling (VSP)
% of the bore hole (based on Homework 1 - Problem 4.0, and
% Homework 2 - Problem 2)
%
function hw0203(epsilon, graph, m, n, WmFR)

clear global; close all;

% Default parameters input
%

% Damping factor
try if isempty(epsilon), epsilon = 5e-4; end; catch epsilon = 5e-4; end;

try if isempty(graph), graph = 1; end; catch graph = 1; end;

% # of elements in "Model vector (# of columns in "Thoery matrix")
try if isempty(m), m = 80; end; catch m = 80; end;

% # of rows in "Theory matrix"  or elements in "Data vector"
        % (or # of samples in depth)
try if isempty(n), n = 160; end; catch n = 160; end;

% Assigning a constant weight in the flat region
try if isempty(WmFR), WmFR = 1.0; end; catch WmFR = 1.0; end;
%

if strcmp(computer,'GLNX86'), root_dir = '/media/sda5'; ...
        else root_dir = 'D:'; end;
gpath = [root_dir '/Users/rilma/work/temp/'];       

fig_resolution = 100;

% Some additional input parameters
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

% Defining terms in WDLS equation
%
Mprior = zeros(m, 1);   % Prior of model vector
We = eye(n, n);         % Weighting error matrix
Wm = eye(m, m);         % Weighting model matrix

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

% Assigning some weights in the flat region
%
WmS = ones(m,m); WmS(ind,:) = WmFR;
Wm = Wm .* WmS;
%

% Applying the forward approach

% Travel-Time or delay [s]
T = G * s;

% Generating a normally distributed noise
N = mu + sigma .* randn(n,1);

% Adding noise to data
TN = T + N;

% 1st derivative of travel time (slowness)
ss1 = diff(TN) ./ deltaZn;

% Inverse methods
%

% Least-Squares solution
LS_Mest = inv(G' * G) * G' * TN;

% Imposing the contraint that the model slowness is 
% equal to the value of s0
if ~isempty(ind), LS_Mest(ind) = s0; end;
%

% Compute data vector from 
TN2 = G * LS_Mest;

% 1st derivative of travel time (slowness)
ss2 = diff(TN2) ./ deltaZn;

%

% Weighted Damped Least Square
%
GW = G' * We * G + epsilon^2 .* Wm;
WDLS_Mest = Mprior + inv(GW) * G' * We * (TN - G * Mprior);

% Imposing the contraint that the model slowness is 
% equal to the value of s0
if ~isempty(ind), WDLS_Mest(ind) = s0; end;
%

% Compute data vector from 
TN3 = G * WDLS_Mest;

% 1st derivative of travel time (slowness)
ss3 = diff(TN3) ./ deltaZn;
%

% Graphical output commands
%
if graph > 0
    
    figh = figure(2); clf(figh, 'reset');    
    
    fig_position = [10 40 500 500];
    set(figh, 'Position', fig_position);
    
    line(z, s, 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2, ...
        'Marker', 'o');
    line(Zn(1:(n-1)), ss1, 'Color', 'r', 'LineStyle', 'none', 'LineWidth', 1, ...
        'Marker', 'o', 'MarkerSize',4);
    line(Zn(1:(n-1)), ss2, 'Color', 'b', 'LineStyle', 'none', 'LineWidth', 1, ...
        'Marker', '*', 'MarkerSize', 3);
    line(Zn(1:(n-1)), ss3, 'Color', 'g', 'LineStyle', 'none', 'LineWidth', 1, ...
        'Marker', '*');
    set(gca, 'XLim', [z0 zf].*m2km, 'XMinorTick', 'on', ...
        'YLim', [1.5 5.5], 'YMinorTick', 'on');                    
%        'YLim', [-2.5 10.5], 'YMinorTick', 'on');
    title('Slowness vs Depth');
    xlabel('[km]'); ylabel('[s/km]');
    legend('Model','Derivative','Least-Squares','WD Least-Squares', ...
        'Location', 'NorthWest');
    grid on;     

    myText = ['# constraints: ' num2str(n,'%3i') ...
        ';  # unknows: ' num2str(m,'%3i') ...
        ';  Sigma: ' num2str(sigma,'%6.4f') ' s'];
    text(0.05, 2.55, ...
        myText, 'FontSize', 16, 'HorizontalAlignment','center', ...
        'Units', 'normalized');
    
    
    % Saving figure
    %
    figfile = ['hw0104_' num2str(n,'%03i') '-' num2str(m,'%03i') '-' ...
        num2str(sigma,'%6.4f')];
    save_figure(figh, graph, fig_resolution, [gpath figfile]);    
    %
end
%