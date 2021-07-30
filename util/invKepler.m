function [E,nu] = invKepler(M,e,verb)
%invKepler - Newton-Raphson inversion of the Kepler Equation
%      invKepler(M,e) calculates the eccentric and true anomalies
%      corresponding to mean anomalies (M) and eccentricites (e).
%
%      INPUTS
%      M    n x 1 or 1 x n vector of mean anomalies.
%      e    vector of n eccentricity values, or scalar eccentricity
%      verb print number of iterations taken if true (default false)
%      
%      OUTPUTS
%      E    n x 1 vector of eccentric anomalies
%      nu   n x 1 vector of true anomalies
%
%      Example:
%      %one year of highly eccentric orbit
%      [E, nu] = invKepler(0:pi/100:2*pi,0.8);
%
% See also: keplerSTM

% Written by Dmitry Savransky, 21 Feb 2011, dsavrans@princeton.edu

%condition inputs
if (numel(M) ~= numel(e)) && numel(e) ~= 1
    error('invKepler:inputError',['For M of length n, e must have',...
        '1 or n elements']);
end
M = M(:);
e = e(:);
if ~exist('verb','var') || isempty(verb)
    verb = false;
end

if numel(e) == 1 && e == 0
    E = M;
else
    %initialize
    counter = 0;
    del = 1;
    E = M./(1-e);
    inds = E > sqrt(6*(1-e)./e);
    if length(e) == 1
        einds = 1;
    else
        einds = inds;
    end
    E(inds) = (6*M(inds)./e(einds)).^(1/3);
    
    %iterate
    while ((del > eps(2*pi)) && (counter <1000))
        E = E - (M - E + e.*sin(E))./(e.*cos(E)-1);
        del = max(abs(M - (E - e.*sin(E))));
        counter = counter+1;
    end
    if (counter == 1000)
        error('invKepler:overIteration',...
            'Maximum number of iterations exceeded');
    end
    if verb
        disp([num2str(counter),' iterations taken']);
    end
end
%calculate true anomaly if needed
if (nargout > 1)
    nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
    inds = nu < 0;
    nu(inds) = nu(inds)+2*pi;
end