function Dates_tree = dates_for_tree(initial_date,expiry_date,expiry,dt,month_steps)
% finds the dates for the trinomial tree
%
% INPUT 
% initial_date:  initial date
% expiry_date:   expiry date
% expiry:        expiry in number of months
% dt:            time interval for the tree
% month_steps:   number of monthly steps for the tree
%
% OUTPUT
% dates_tree:    dates in the tree

% mean number of days in a month (between settlement date and expiry date)
mean_number_of_days = (expiry_date - initial_date)/expiry; 

% dates for the tree
fun = @(i) addtodate(initial_date,round(dt*i*mean_number_of_days),'day');  
I = 1:expiry*month_steps;
Dates_tree = arrayfun(fun,I); 

end 

