function [R, flag] = blockchol(A)
%BLOCKCHOL   Block Cholesky factorization for the GPU.
%   R = BLOCKCHOL(A) calculates the Cholesky factor of full A using
%   the diagonal and upper triangle of A, such that A = R'*R. A must be
%   positive definite, and the lower triangle is assumed to be the (complex
%   conjugate) transpose of the upper triangle.
%
%   This takes a potentially large matrix on the CPU, and internally
%   converts it to blocks of matrices on the GPU that are small enough to
%   fit into the 32-bit memory. The algorithm is applied to these blocks,
%   and then the result is converted back to a matrix on the CPU.

%   Copyright 2025 The MathWorks, Inc.

flag = 0;

n = size(A, 1);
blocksize = 2^14; % 2^14 = 16384 = sqrt(2^32/(2*8))
                  % should allow a blocksize-by-blocksize complex double
                  % matrix to fit in 32-bit memory.
sz = blocksize*ones(1, ceil(n/blocksize));
sz(end) = n - sum(sz(1:end-1));


Rcell = makeTriuGpuarrayCell(A, blocksize);

nb = numel(sz);

for j=1:nb
    for k=1:j-1
        for l=j:nb
            Rcell{j, l} = Rcell{j, l} - Rcell{k, j}'*Rcell{k, l};
        end
    end

    [Rcell{j, j}, flag] = chol(Rcell{j, j});
    if flag ~= 0
        % We don't give a specific column, just a generic flag.
        flag = -1;
        break;
    end

    for k=j+1:nb
        Rcell{j, k} = Rcell{j, j}' \ Rcell{j, k};
    end
end

R = gatherTriuGpuarrayCell(Rcell, n, blocksize);

if nargout == 1 && flag ~= 0
    error(message('MATLAB:posdef'));
end
end

function Rcell = makeTriuGpuarrayCell(R, blocksize)
n = size(R, 1);
nb = ceil(n/blocksize);

Rcell = cell(nb,nb);
for i=1:nb
    for j=i:nb
        Rij = R((i-1)*blocksize+1:min(n,i*blocksize), ...
            (j-1)*blocksize+1:min(n,j*blocksize));
        Rij = gpuArray(Rij);
        Rcell{i, j} = Rij;
    end
    Rcell{i, i} = triu(Rcell{i, i});
end
end

function R = gatherTriuGpuarrayCell(Rcell, n, blocksize)
R = zeros(n, n, "like", gather(Rcell{1, 1}([])));
n = size(R, 1);
nb = ceil(n/blocksize);

for i=1:nb
    for j=i:nb
        Rij = Rcell{i, j};
        Rij = gather(Rij);
        R((i-1)*blocksize+1:min(n,i*blocksize), ...
            (j-1)*blocksize+1:min(n,j*blocksize)) = Rij;
    end
end
end
