%% Device Setting (NIdaqUSB6361:sound device)
s = daq.createSession('ni');% Create a session with NI device
addAnalogOutputChannel(s,'dev2',0:1,'Voltage');% ch0_sound stimuli, ch1_trigger
SamplingRate = 10e4; s.Rate = SamplingRate;  % set the sampling rate (Max 2MHz)
% s.NotifyWhenScansQueuedBelow = 1000;
Fs = s.Rate;
% s.IsContinuous = false;
%% click sound
duration_on = 0.05e-3;  %[sec]  %0.2e-3
duration_off= 0.01e-3;  %[sec]
dur = 1e-2;
amp = 0.2;  %[volt]

mode = 5;  %1:periodic, 2:4 chunk, 3:8 chunk, 4:16 chunk, 5:64 chunk

N1=duration_off*SamplingRate; 
N2=duration_on*SamplingRate;
trg=dur*SamplingRate; 
click = amp*[zeros(N1,1); ones(N2,1);-ones(N2,1);zeros(N1-N2*2+trg,1)];
trigger =1*[zeros(N1,1); ones(trg,1); zeros(N1,1)];


%% random configuration

N = 64;
chunk = struct();
chunk2 = struct();
dur = 1e-2;

freqs = [32,16,8,4,2,1];
times = [15,5];
itis = freqs.^(-1) - dur;

count = 1;
for i = 1:length(freqs)
    for i2 = 1:length(times)
        rep = freqs(i)*times(i2);
        tmp = [];
        for j = 1:round(rep)
            tmp = [tmp, itis(i)];
        end
        chunk.("m"+count) = tmp;
        chunk2.("m"+count) = strrep("p"+freqs(i)+"_"+times(i2),".","_");
        count = count + 1;
    end
end

freqs2 = [32,16,8,4];
itis2 = (freqs2.^(-1) - dur).';
itis2_ = [itis2,itis2,itis2+freqs2.^(-1).'];
for i = 1:length(freqs2)
    for i2 = 1:length(times)
        rep = freqs2(i)*times(i2)/4;
        tmp = [];
        for j = 1:round(rep)
            tmp = [tmp, itis2_(i,:)];
        end
        chunk.("m"+count) = tmp;
        chunk2.("m"+count) = strrep("f"+freqs2(i)+"_"+times(i2),".","_");
        count = count + 1;
    end
end

itis3_ = [itis2,itis2,itis2,itis2,itis2+freqs2.^(-1).',itis2+freqs2.^(-1).'];
for i = 1:length(freqs2)
    for i2 = 1:length(times)
        rep = round(freqs2(i)*times(i2)/8);
        tmp = [];
        for j = 1:max([2,round(rep)])
            tmp = [tmp, itis3_(i,:)];
        end
        chunk.("m"+count) = tmp;
        chunk2.("m"+count) = strrep("e"+freqs2(i)+"_"+times(i2),".","_");
        count = count + 1;
    end
end


%%
groupN = 10;

mode = size(struct2table(chunk),2);
rest = 5;
spon = 10;

% start serial communication
serialobj=serial('com3');
set(serialobj, 'BaudRate', 115200);
fopen(serialobj);
% start recording
d=hex2dec('FD');
fwrite(serialobj, d, 'uint8');

disp("periodic_click2")
for k=1:groupN
    xx = randperm(mode);
    disp([k,0,xx])
    s.queueOutputData([zeros(spon*Fs,1),ones(spon*Fs,1)*0.5]); % sound set: require some time
    s.startBackground
    pause(spon+2);
    for i = 1:mode
        ITIs = chunk.("m"+string(xx(i)));
        tones = cell(length(ITIs)*2+2,1); 
        TRG = cell(length(ITIs)*2+2,1);
        count = 2;
        tones{1} = zeros(round(SamplingRate),1);
        TRG{1}   = zeros(round(SamplingRate),1);
        for j=1:length(ITIs)
            tones{count} = click; 
            TRG{count} = trigger;
            count = count + 1;
            tones{count} = zeros(round(ITIs(j)*SamplingRate),1);
            TRG{count}   = zeros(round(ITIs(j)*SamplingRate),1);
            count = count + 1;
        end
        tones{count} = zeros(round(rest*SamplingRate),1);
        TRG{count}   = zeros(round(rest*SamplingRate),1);
        count = count + 1;
        tones_ = cell2mat(tones);
        TRG_ = cell2mat(TRG);
        %disp([i,chunk2.("m"+string(xx(i))),mode])
        s.queueOutputData([tones_,TRG_]); % sound set: require some time
        s.startBackground
        while ~s.IsDone
            pause(1);
        end
    end
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

