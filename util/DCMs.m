function rotMats = DCMs()
%Return array of direction cosine matrices about body principal axes
%
%INPUT
%   None
%
%OUTPUT
%   rotMats - 3x1 cell of lambda functions corresponding to 1st, 2nd and
%   3rd axis rotations
%
%NOTE
%   Rotations are CCW positive - these are {}^B{C}^{I} matrices.

% Copyright (c) 2019 Dmitry Savransky (ds264@cornell.edu)

rotMatE3 = @(ang) [cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
rotMatE2 = @(ang) [cos(ang) 0 -sin(ang);0 1 0;sin(ang) 0 cos(ang)];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];
rotMats = {rotMatE1,rotMatE2,rotMatE3};

end