function [c, ceq] = constraints(parameters)
% sets constraints for fmincon function
%
% INPUT 
% parameters:   vector of parameters for fmincon
%
% OUTPUT
% c:            inequality constraints 
% ceq:          equality constraints

% inequality constraints
c = [-parameters(1); -parameters(2)];

% equality constraints
ceq = [];

end