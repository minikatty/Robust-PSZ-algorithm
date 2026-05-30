function n = pagenorm_tmp(X, type)
%PAGENORM Minimal compatibility for MATLAB R2019b.
% Supports: n = pagenorm(X,'fro') for 3-D arrays.
% Returns: n as 1x1xK so that squeeze(n) gives Kx1 or 1xK.

if nargin < 2
    type = 'fro';
end
type = lower(string(type));
if type ~= "fro"
    error('pagenorm:Type','Only ''fro'' is supported in this minimal implementation.');
end

if ndims(X) ~= 3
    error('pagenorm:Dim','X must be a 3-D array in this minimal implementation.');
end

% Frobenius norm per page: sqrt(sum(|X|^2))
n = sqrt( sum(sum( abs(X).^2, 1), 2 ) );  % size: 1x1xK
end
