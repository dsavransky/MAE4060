function A = skew(x)
%Calculate the skew-symmetric matrix of a vector
%
%INPUT
%   x (1x3 float): vector components
%
%OUTPUT
%   A (3x3 float): skew-symmetric (cross-product equivalent) matrix of the
%                  vector input
%
%NOTE
%   Rotations are CCW positive - these are {}^B{C}^{I} matrices.

% Copyright (c) 2019 Dmitry Savransky (ds264@cornell.edu)

A = [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];

end