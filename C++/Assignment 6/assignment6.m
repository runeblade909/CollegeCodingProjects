%Assignment 6

clc
clear

choice = 1;

compressedExists = false;

            audiovector = 0;

while choice ~= 0
    choice = menu('What would you like to do?',...
        'Choose File',...
        'Play audio',...
        'Audio echo',...
        'Compress file',...
        'fftBar',...
        'Play Compressed Audio');
    switch choice
        case 1
            file = input('What audio file are we working with today?\n','s');
            [audioVec,Fs] = audioread(strcat(file,'.wav'));
        case 2
            player = audioplayer(audioVec,Fs);
            if audiovector == 0
            play(player);
            audiovector = 1;
            else
            stop(player);
            end
        case 3
            delay = input('Enter delay in seconds: ');
            echoGain = input('Enter gain between 0-1.5: ');
            echo = audioEcho(audioVec,Fs,delay,echoGain);
            player = audioplayer(echo,Fs);
            play(player);
        case 4
            [compressedAudioVec,newFs] = compressMono(audioVec,Fs);
            compressedExists = true;
        case 5 
            fftBar(audioVec,Fs)  
        case 6
            if compressedExists
                
                player = audioplayer(compressedAudioVec,newFs);
                play(player);
            else
                disp('No compressed audio to play.');
            end
    end
end