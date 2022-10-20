function[tH,tM,tS,format,timeOfDay] = convertTime(tH,tM,tS,format,timeOfDay)
%12 hour format and conversion
    if format == 12
        if timeOfDay == 1
            tH= tH + 12;
        elseif timeOfDay ==0 && tH == 12
                tH= tH-12;
           
        elseif timeOfDay == 0
            
            tH= tH;
       
        end
        if tH>=24
                tt=tH;
                tH=tt-12;
        end
    end
%24 hour format and conversion
    if format == 24
        if tH > 12
            tHour=tH;
            tH= tHour - 12;
            timeOfDay = 1;  
        elseif tH == 0
            tHour=tH;
            tH= tHour +12;
            timeOfDay = 0;     
        elseif tH < 12
            tHour=tH;
            timeOfDay = 0;
        end
    end




end