% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function fx = scalar_Gaussian(x, mu, sigma)

if size(x,2)~=1
    error('should be a column vector')
end
if ~isscalar(mu) ||  ~isscalar(sigma)
    error('should be scalars')
end

fx = exp(-(x-mu).^2/(2*sigma^2))/(sigma*sqrt(2*pi));
end