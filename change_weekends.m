function [dates] = change_weekends(input_dates)
% changes weekend dates in mondays
%
% INPUT
% input_dates:  vector of dates
%
% OUTPUT
% dates:        new vector of shifted dates

% allocation of dates
dates = input_dates;

while sum(isweekend(datetime(dates,'ConvertFrom','datenum'))) > 0  % check if there are some weekend days
    for i = 1:length(dates)
        day = datetime(dates(i),'ConvertFrom','datenum');          % conversion
        if isweekend(day) == 1                                     % if it is a weekend day
            if weekday(day) == 7                                   % saturday
                dates(i) = dates(i) +2;
            elseif weekday(day) == 1                               % sunday
                dates(i) = dates(i) +1;
            end
        end
    end
end

end