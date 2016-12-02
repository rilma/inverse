
function [m, sv] = calc_svd_model(G, T, p)

% Range of G matrix (if no input, default value computed)
    try if isempty(p), p = rank(G); end; catch p = rank(G); end;

% Computing the matrix singular value decomposicion of G
    [U,S,V] = svd(G);

% Singular values (diagonal elements of S matrix)
%
    sv = diag(S);
    
% Compact Form
%

% 1st 'p' columns of U
    Up = U(:,1:p);

% 1st 'p' columns of V
    Vp = V(:,1:p);

% 1st 'p' columns and rows of S
    Sp = S(1:p,1:p);
%

% least-squares
    m = Vp * pinv(Sp) * Up' * T;
    