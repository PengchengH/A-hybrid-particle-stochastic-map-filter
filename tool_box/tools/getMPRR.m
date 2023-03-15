function TH = getMPRR(Ne)
% constructs a random orthogonal matrix that has ones(Ne,1)
% as an eigenvector, i.e. a Mean Preserving Random Rotation
% Technically it's only a rotation if the determinant is 1,
% but the name is catchy and we don't want to be pedantic
% This implementation is not very efficient. It would be
% faster to do 3 mat/vec multiplications instead of 2 mat/mat
% multiplications followed by a mat/vec.

% persistent U
% if(isempty(U))
  A = [ones(Ne,1) randn(Ne,Ne-1)];
  [U,~] = qr(A);
% end
[Up,~] = qr(randn(Ne-1));
TH = U*[1 zeros(1,Ne-1);zeros(Ne-1,1) Up]*(U');
