


% Example singular values
sigma = 1./linspace(1, 100).^2;

% Singular values of pseudo-inverse, where SVD is truncated at tol
tol = 0.01;
sigmaTruncate = 1./sigma;
sigmaTruncate(sigma <= tol*sigma(1)) = 0;

% Singular values of inverse computed using Tikhonov regularization
lambda = 0.01;
sigmaTikhonov = sigma ./ (sigma.^2 + lambda.^2);
plot([sigmaTruncate(:), sigmaTikhonov(:)])
title('Singular values of inverse')
legend('Pseudoinverse (truncated SVD)', 'Tikhonov regularization')