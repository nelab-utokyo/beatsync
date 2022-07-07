%%

music = struct();

[y,Fs] = audioread('./material/sonata_s75.wav');
disp(mean(y.*y))
ymean = mean(y.*y);
music.y0 = y.*sqrt(ymean./mean(y.*y));

[y,Fs] = audioread('./material/Sonata_2r_org.wav');
disp(mean(y.*y));
music.y1 = y.*sqrt(ymean./mean(y.*y));

[y,Fs] = audioread('./material/sonata_s200.wav');
disp(mean(y.*y))
music.y2 = y.*sqrt(ymean./mean(y.*y));
[y,Fs] = audioread('./material/sonata_s300.wav');
disp(mean(y.*y))
music.y3 = y.*sqrt(ymean./mean(y.*y));
[y,Fs] = audioread('./material/sonata_s400.wav');
disp(mean(y.*y))
music.y4 = y.*sqrt(ymean./mean(y.*y));

%}
%[y,Fs] = audioread('./material/sonata_2r_p100.wav');
%music.y6 = y.*sqrt(ymean./mean(y.*y));
%[y,Fs] = audioread('./material/sonata_2r_r100.wav');
%music.y7 = y.*sqrt(ymean./mean(y.*y));


%%
s = daq.createSession('ni');% Create a session of NI device
addAnalogOutputChannel(s,'dev2',0:1,'Voltage');% Add an analog output channel(sound)
% Set session and channel properties
s.Rate = Fs;% set the sampling rate to 500kHz
s.IsContinuous = false;% you can stop the operation in the middle
%below=10000; s.NotifyWhenScansQueuedBelow = below;
rest = 5;
tones = cell(6*2,1); TRG = cell(6*2,1);

sponlen = 10*Fs;

tones{1} = zeros(sponlen,1);
temp = ones(sponlen,1);
temp(0.01*Fs:sponlen-0.01*Fs,1) = 0;
TRG{1} = temp;

tones{2} = zeros(round(rest*Fs),1);
TRG{2} = zeros(round(rest*Fs),1);

tones_ = cell2mat(tones);
TRG_ = cell2mat(TRG);

%%

%lh = addlistener(s,'DataRequired',@(src,event) src.queueOutputData([tones, TRG]));
disp("fixed_music_ex2")

% start serial communication
serialobj=serial('com3');
set(serialobj, 'BaudRate', 115200);
fopen(serialobj);
% % start recording
% d=hex2dec('FD');
% fwrite(serialobj, d, 'uint8');


for n = 1:10
% start recording
d=hex2dec('FD');
fwrite(serialobj, d, 'uint8');

xx = randperm(length(fieldnames(music)));
    for i = 2:(length(fieldnames(music))+1)
        tones{i*2-1} = music.("y"+string(xx(i-1)-1))(:,1);
        TRG{i*2-1} = music.("y"+string(xx(i-1)-1))(:,1);
        tones{i*2} = zeros(round(rest*Fs),1);
        TRG{i*2} = zeros(round(rest*Fs),1);
    end
    tones_ = cell2mat(tones);
    TRG_ = cell2mat(TRG);
    
    s.queueOutputData([tones_,TRG_]); % sound set: require some time
    s.startBackground
    disp([n,0,xx])
    while ~s.IsDone
        pause(5);
    end
%     % end recording
%     d=hex2dec('FE');
%     fwrite(serialobj, d, 'uint8');

end

% end serial communication
d=hex2dec('FE');
fwrite(serialobj, d, 'uint8');
fclose(serialobj);
delete(serialobj);

%%

% delete(s);
s.stop

% end serial communication
d=hex2dec('FE');
fwrite(serialobj, d, 'uint8');
fclose(serialobj);
delete(serialobj);


