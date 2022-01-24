
function [complete] = CF_latency(patha)
    clearvars FirstSpikeRaster FirstSpikeLatency
    files = dir(patha + "/TC_*/*RasterAndFSL.mat");
    load(files(1).folder + "/" +files(1).name)
    disp(files(1).name)
    elenum = 96;
    dbnum = 7;
    frqnum = 18;

    CF = zeros(elenum,1);
    latency = zeros(elenum,1);


    for e = 1:elenum
        firstspikenum = zeros(frqnum,2);
        for f = 1:frqnum
            tmp = 0;
            tmp2 = 0;
            for d = 1:dbnum
                tmp = tmp + sum(sum(FirstSpikeRaster{d,f,e}));
                tmp2 = tmp2 + sum(sum(Raster{d,f,e}));
            end
            firstspikenum(f,:) = [tmp,tmp2];
        end
        while 1
            [M,I] = max(firstspikenum(:,1));
            I2 = find(firstspikenum(:,1) == M);
            if M < 10
                [M,I] = max(firstspikenum(:,2));
                I2 = find(firstspikenum(:,2) == M);
            end
            if max(FirstSpikeLatency(1,I2,e)) > 2
                if length(I2) > 1
                    cand = FirstSpikeLatency(1,I2,e);
                    cand(cand == 0) = NaN;
                    [Mn,I3] = min(cand);
                    I4 = find(cand == Mn);
                    CF(e) = mean(I2(I4));
                    latency(e) = Mn;
                else
                    CF(e) = I;
                    latency(e) = FirstSpikeLatency(1,I,e);
                end
                break
            else
                firstspikenum(I,:) = NaN;
                if sum(isnan(firstspikenum(:,1))) == frqnum
                    CF(e) = NaN;
                    latency(e) = NaN;
                    break
                end
            end
        end
    end

    save(files(1).folder+"/CF_latency.mat","CF","latency")
    complete = true;
end
