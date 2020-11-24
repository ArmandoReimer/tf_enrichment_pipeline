function image = bound2im(b, M, N)

if size(b, 2) ~=2
    error('The boundary must be of size np-by-2');
end

b = round(b);

%Defaults
if nargin==1
    Mmin = min(b(:, 1)) - 1;
    NMin = min(b(:, 2)) - 1;
    H = max(b(:, 1)) - min(b(:,1)) + 1;
    W = max(b(:, 2)) - min(b(:, 2)) + 1;
    M = H + Mmin;
    N = W + NMin;
end

image = false(M,N);
linearIndex = sub2ind([M,N], b(:, 1), b(:, 2));
image(linearIndex) = 1;
  
end