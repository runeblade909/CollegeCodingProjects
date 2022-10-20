function [] = displayTime(tH,tM,tS,format,timeOfDay)
if timeOfDay == 0
    timeOfDay = 'AM';
elseif timeOfDay == 1
    timeOfDay = 'PM';
end

if format == 12
    format = 24;
    timeOfDay=[];
elseif format == 24
    format = 12;
end



fprintf('The time is %1.0f:%1.0f:%1.0f %s in a %1.0f hour format.',tH,tM,tS,timeOfDay,format)

end