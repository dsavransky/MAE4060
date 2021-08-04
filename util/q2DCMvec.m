function C = q2DCMvec(q)
%Calculate direction cosine matrices from quaternions
%
%INPUT
%   q (nx4 float): Quaternions, packaged as [v,r] (vector part first)
%
%OUTPUT
%   C (3x3,n float): Direction cosine matrices
%
%NOTE
%   Rotations are CCW positive - these are {}^B{C}^{I} matrices.

% Copyright (c) 2019 Dmitry Savransky (ds264@cornell.edu)

q1 = q(:,1);
q2 = q(:,2);
q3 = q(:,3);
q4 = q(:,4);

C = zeros(3,3,length(q1));

C(1,1,:) = q1.^2 - q2.^2 - q3.^2 + q4.^2;
C(2,2,:) = -q1.^2 + q2.^2 - q3.^2 + q4.^2;
C(3,3,:) = -q1.^2 - q2.^2 + q3.^2 + q4.^2;
C(2,1,:) = 2*(q1.*q2 - q3.*q4);
C(3,1,:) = 2*(q1.*q3 + q2.*q4);
C(1,2,:) = 2*(q1.*q2 + q3.*q4);
C(1,3,:) = 2*(q1.*q3 - q2.*q4);
C(2,3,:) = 2*(q1.*q4 + q2.*q3);
C(3,2,:) = 2*(-q1.*q4 + q2.*q3);

end