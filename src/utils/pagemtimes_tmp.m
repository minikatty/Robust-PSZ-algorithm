function C = pagemtimes_tmp(A, opA, B, opB)
%PAGEMTIMES  Minimal compatibility implementation for MATLAB R2019b.
% Supports calls like:
%   C = pagemtimes(A,'ctranspose',B,'none')
%   C = pagemtimes(A,'transpose', B,'none')
%   C = pagemtimes(A,'none',     B,'none')
%
% A: [m,n,K], B: [n,p,K] (after applying opA/opB per page)

if nargin ~= 4
    error('pagemtimes:InvalidInput','Expected 4 inputs: A, opA, B, opB.');
end

A = apply_op(A, opA);
B = apply_op(B, opB);

if ndims(A) ~= 3 || ndims(B) ~= 3
    error('pagemtimes:Dim','A and B must be 3-D arrays.');
end
if size(A,3) ~= size(B,3)
    error('pagemtimes:Pages','A and B must have the same number of pages.');
end
if size(A,2) ~= size(B,1)
    error('pagemtimes:InnerDim','Inner dimensions must agree after ops.');
end

K = size(A,3);
C = zeros(size(A,1), size(B,2), K, 'like', A);

for k = 1:K
    C(:,:,k) = A(:,:,k) * B(:,:,k);
end
end

function X = apply_op(X, op)
op = lower(string(op));
switch op
    case "none"
        % no-op
    case "transpose"
        X = permute(X, [2 1 3]);
    case "ctranspose"
        X = conj(permute(X, [2 1 3]));
    otherwise
        error('pagemtimes:Op','Unsupported op "%s".', op);
end
end
