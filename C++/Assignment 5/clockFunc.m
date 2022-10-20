function [tH,tM,tS,format,timeOfDay] = clockFunc()

clc

tH=input('What is the hour?(0-23)\n');
    if (tH) < 0
        error(' I need a positive number')
    elseif (tH) > 23
        error('I gave you a range fool!')
    end
tM=input('What is the minute?(0-59)\n');
    if (tM) < 0
        error(' I need a positive number')
    elseif (tM) > 59
        error('I gave you a range fool!')
    end
tS=input('What is the second?(0-59)\n');
    if (tS) < 0
        error(' I need a positive number')
    elseif (tS) > 59
        error('I gave you a range fool!')
        
    end

 
format= input('Is that in 24 or 12 hour format? (enter 24 or 12)\n');
%12 hour format and conversion
    if format == 12
    timeOfDay = input(' Is that AM or PM? (0=AM, 1=PM)\n');
        if timeOfDay == 1
            tH24= tH + 12;
        elseif timeOfDay ==0 && tH == 12
                tH24= tH-12;
           
        elseif timeOfDay == 0
            
            tH24= tH;
       
        end
        if tH24>=24
                tt=tH24;
                tH24=tt-12;
        end
    fprintf('The time is %1.0f:%1.0f and %.2f seconds.\n',tH,tM,tS);
    fprintf('In 24 hour format, its %1.0f:%1.0f and %.2f seconds.\n',tH24,tM,tS);
    end
%24 hour format and conversion
    if format == 24
        if tH > 12
            tHour=tH;
            tH= tHour - 12;
            t='PM'
            timeOfDay=string(t);
        fprintf('In 12 hour format, its %1.0f:%1.0f PM and %.2f seconds.\n',tH,tM,tS);     
        elseif tH == 0
            tHour=tH;
            tH= tHour +12;
            t='AM';
            timeOfDay=string(t);
       fprintf('In 12 hour format, its %1.0f:%1.0f AM and %.2f seconds.\n',tH,tM,tS);     
        elseif tH < 12
            tHour=tH;
            fprintf('In 12 hour format, its %1.0f:%1.0f AM and %.2f seconds.\n',tH,tM,tS);
            t='AM';
            timeOfDay=string(t);
        end
     fprintf('The time is %1.0f:%1.0f and %.2f seconds.\n',tHour,tM,tS); 
    end
    
    

end