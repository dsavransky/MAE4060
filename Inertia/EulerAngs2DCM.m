function A = EulerAngs2DCM(angs,rotSet,space)
%Calculate the direction cosine matrix A associated with a specific ordered 
%rotation specified by Euler angles angs with axis order rotSet using
%either body or space axes.
%
%INPUT
%   A         3x1 array or Euler angles (rad)
%   rotSet    3x1 array of rotation axes to use
%   space     logical, true for Space rotation, false for Body (default
%             false)
%
%OUTPUT
%   A         Direction Cosine Matrix
%
%EXAMPLE
% %body-3 1-2-3 Euler angles of 3x3 DCM matrix dcm:
% [th1,th2,th3] = calcEulerAngs(dcm,[1,2,3],false);
%
%NOTE
%   Rotations are CCW positive - these are {}^B{C}^{I} matrices.

%Copyright (c) 2019 Dmitry Savransky (ds264@cornell.edu)

%assume Body rotations unless told otherwise
if ~exist('space','var'),space = false;end

%input check
assert(numel(angs)==3,'EulerAngs2DCM:inputCheck',...
    'angs must have 3 elements')
assert(numel(rotSet)==3,'EulerAngs2DCM:inputCheck',...
    'rotSet must have 3 elements')

rotMats = DCMs(); %get DCM lambdas

if space
    rotSet = flip(rotSet);
end
    
A = rotMats{rotSet(3)}(angs(3))*rotMats{rotSet(2)}(angs(2))*...
    rotMats{rotSet(1)}(angs(1));

end