%%
%
paths = ["20200908","20200910","20200915","20200917","20201006","20210106","20210107"] + "/";
dpath = "/Users/TP/research/Tlab/data/";
spath = "/Users/TP/research/Tlab/sound/material/";
scand = [75,100,200,300,400];
wsec = 10;
numker = 499;
numband = 31;
numhist = 40;

trainy = [7,9,3,1,5,4,2,6,10,8];  %randam value for training dataset
div = 5;

%%

for d = 1:length(paths)
    load(dpath + paths(d) + "shuf_param.mat")
    load(dpath + paths(d) + "CF_latency.mat")
    
    load(dpath + paths(d) + "processed/exdata_ex1.mat")
    muas = struct();
    %
    rmv = find(latency(TC) > 20);
    TC(rmv) = [];
    init = 1500;
    
    for usenum = 1:div
        trainind = trainy;
        trainind(usenum*2-1:usenum*2) = [];
        for c = 1:length(scand)
            data = df.("c"+scand(c));
            [resp,~,~] = data2matrix3(data,1:length(df.("c"+scand(c)).trg),TC,1);
            muas.("ctr"+usenum+scand(c)) = squeeze(sum(mean(resp(:,trainind,:),2)));
            bin = [init-wsec, init];
            tmua = [];
            while bin(2) < length(muas.("ctr"+usenum+scand(c)))
                tmp = mean(muas.("ctr"+usenum+scand(c))(bin(1):bin(2)-1));
                tmua = [tmua, tmp];
                
                bin = bin + wsec/2;
            end
            muas.("ttr"+usenum+scand(c)) = tmua';
            muas.("cts"+usenum+scand(c)) = squeeze(sum(mean(resp(:,trainy(usenum*2-1:usenum*2),:),2)));
            
            bin = [init-wsec, init];
            tmua = [];
            while bin(2) < length(muas.("cts"+usenum+scand(c)))
                tmp = mean(muas.("cts"+usenum+scand(c))(bin(1):bin(2)-1));
                tmua = [tmua, tmp];
                bin = bin + wsec/2;
            end
            muas.("tts"+usenum+scand(c)) = tmua';
        end
        
        %
        for c = 1:length(scand)
            if c == 2
                [y, Fs] = audioread(spath + "sonata_2r_org.wav");
            else
                [y, Fs] = audioread(spath + "sonata_s" + scand(c) + ".wav");
            end    
            window = Fs*wsec/1000;
            noverlap = Fs*wsec/2000;
    
            S = melSpectrogram(mean(y,2),Fs, ...
                               'Window',hann(window,'periodic'), ...
                               'OverlapLength',noverlap, ...
                               'NumBands',numband, ...
                               'FrequencyRange',[1000,32e3]);
            
            S(S < mean(S(:,end))) = mean(S(:,end));
            S = log10(S/mean(S(:,end)));
            %imagesc(S);
            
            %make adaptation kernel
            kernel = [];
            for n = 1:numband
                kernel = [kernel;(1:numker)/(500 - 105*log10(1000*((2^((n-1)/6)))))];
            end
            
            kernel = exp(kernel);
            kernel = kernel./sum(kernel,2);
            
            S2 = zeros(size(S));
            %for q = 2:length(S)
            %    tmp = max([1,q-numker+1]);
            %    S2(:,q) = dot(kernel(:, numker-(q-1-tmp):numker)',S(:,tmp:q-1)');
            %end
            
            S3 = S - S2;
            S3(S3 < 0) = 0;
            
            %make tensor
            Stens = zeros(length(S3), numband, numhist);
            for n = 1:length(S)
                tmp = max([n - numhist + 1, 1]);
                Stens(n,:,end - n + tmp:end) = reshape(S3(:,tmp:n),[1,numband,n - tmp + 1]);
            end
            
            %fitting
            
            Lambda = [0.001,0.01,0.1,1,10,100,1000];
            tmp = min([length(muas.("ttr"+usenum+scand(c))), length(Stens)]);
            ydata = muas.("ttr"+usenum+scand(c))(1:tmp);
            xdata = reshape(Stens(1:tmp,:,:),tmp,[]);
            
            CVMdl = fitrlinear(xdata,ydata,'ObservationsIn','rows','KFold',5,'Lambda',Lambda,...
                'Learner','leastsquares','Solver','sgd','Regularization','ridge');
            
            mse = kfoldLoss(CVMdl);
            
            %plot result
            
            [M,I] = min(mse);
            disp([mse,Lambda(I)])
            %{
            figure()
            imagesc(reshape(CVMdl.Trained{1}.Beta(:,I),[numband, numhist]))
            colorbar("Fontsize",15)
            xticks((1:10:numhist+1)*(numhist-1)/numhist)
            xticklabels(-wsec/2*numhist:wsec/2*10:0)
            
            yticks(1:6:numband+1)
            yticklabels(1*((2.^((0:6:numband+1)/6))))
            xlabel("time [ms]","Fontsize",25)
            ylabel("Frequency [kHz]","Fontsize",25)
            title(scand(c) + "%", "Fontsize", 25)
            ax = gca;
            ax.YAxis.FontSize = 20;
            ax.XAxis.FontSize = 20;
            mkdir(dpath + paths(d) + "summary")
            saveas(gca, dpath + paths(d) + "summary/STRF_"+scand(c) + ".png")
            %}
            muas.("b"+usenum+scand(c)) = reshape(CVMdl.Trained{1}.Beta(:,I),[numband, numhist]);
        end
    end
    mkdir(dpath + paths(d) + "buildmat")
    save(dpath + paths(d) + "buildmat/mua_summary3.mat", "muas")
end


%% average of all animals

strfs = zeros(numband, numhist,length(scand),length(paths),div);

for d = 1:length(paths)
    load(dpath + paths(d) + "buildmat/mua_summary3.mat", "muas")
    for c = 1:length(scand)
        for b = 1:div
            strfs(:,:,c,d,b) = zscore(muas.("b"+b+scand(c)),0,"all");
        end
    end
end

%% plot summary

fpath = "/Users/TP/research/Tlab/data/STRF/";

for c = 1:length(scand)
    figure()
    imagesc(mean(mean(strfs(:,:,c,:,:),5),4))
    colorbar("Fontsize",15,"Ticks",-5:1:5)
    xticks((1:10:numhist+1)*(numhist-1)/numhist)
    xticklabels(-wsec/2*numhist:wsec/2*10:0)
    
    yticks(1:6:numband+1)
    yticklabels(1*((2.^((0:6:numband+1)/6))))
    xlabel("time [ms]","Fontsize",25)
    ylabel("Frequency [kHz]","Fontsize",25)
    title(scand(c) + "%", "Fontsize", 25)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    saveas(gca, fpath + "/zSTRF_"+scand(c) + ".png")
end


%% make music input
spect = struct();

for c = 1:length(scand)
    if c == 2
        [y, Fs] = audioread(spath + "sonata_2r_org.wav");
    else
        [y, Fs] = audioread(spath + "sonata_s" + scand(c) + ".wav");
    end
    window = Fs*wsec/1000;
    noverlap = Fs*wsec/2000;
    numband = 31;
    numhist = 40;
    
    S = melSpectrogram(mean(y,2),Fs, ...
                       'Window',hann(window,'periodic'), ...
                       'OverlapLength',noverlap, ...
                       'NumBands',numband, ...
                       'FrequencyRange',[1000,32e3]);
    
    S(S < mean(S(:,end))) = mean(S(:,end));
    S = log10(S/mean(S(:,end)));
    spect.("s"+scand(c)) = S;
end

%% make mua output and strf 

ydata_ = struct();
strfs_ = struct();
for d = 1:length(paths)
    load(dpath + paths(d) + "buildmat/mua_summary3.mat", "muas")
    load(dpath + paths(d) + "shuf_param.mat", "TC")
    for c = 1:length(scand)
        tmp = length(spect.("s"+scand(c)))-160;
        if d == 1
            ydata_.("d"+scand(c)) = zeros(length(paths),tmp,div);
        end
        for b = 1:div
            ydata_.("d"+scand(c))(d,:,b) = muas.("tts"+b+scand(c))(1:tmp)/length(TC);
        end
    end
    for c = 1:length(scand)
        if d == 1
            strfs_.("d"+scand(c)) = nan(numband, numhist, length(paths), length(scand), div);
        end
        for c2 = 1:length(scand)
            for b = 1:div
                strfs_.("d"+scand(c))(:,:,d,c2,b) = muas.("b"+b+scand(c2));
            end
        end
        strfs_.("d"+scand(c)) = squeeze(nanmean(nanmean(strfs_.("d"+scand(c)),4),3));
    end
end

%%

fun = @(x,xdata)(x(1) + (x(2)./exp(-(xdata-x(3))/x(4))));
opt = optimoptions('lsqcurvefit','MaxIterations',10000);
xdata_ = struct();

fitres = struct();
for c = 1:length(scand)
    for b = 1:div
        ydata = squeeze(mean(ydata_.("d"+scand(c))(:,:,b)));
        mstrf = squeeze(strfs_.("d"+scand(c))(:,:,b));
        xdata = zeros(1, length(ydata));
        numhist_ = 40;
        for q = 2:length(ydata)
            tmp = max([1,q-numhist_+1]);
            xdata(q) = sum(sum(mstrf(:, numhist_-(q-1-tmp):numhist_).*spect.("s"+scand(c))(:,tmp:q-1)));
        end
        xdata_.("x" + b + scand(c)) = xdata;
        [x,resnorm,residual,exitflag,output] = lsqcurvefit(fun,[-76,110,2000,10500],xdata,ydata,-[1,1,1,1]*1e+5,[1,1,1,1]*1e+5,opt);
        times = 1:length(xdata);
        figure()
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.7, 0.55]);
        plot(times,ydata,'k-',times,fun(x,xdata),'c-')
        fitres.("f"+b+scand(c)).x = x;
        fitres.("f"+b+scand(c)).err = resnorm/length(residual);
        fitres.("f"+b+scand(c)).R2 = 1-resnorm/sum((ydata - mean(ydata)).^2);
        if c < 4
            xticks(0:2000:length(ydata))
            xticklabels(0:10:length(ydata)/20)
        else
            xticks(0:1000:length(ydata))
            xticklabels(0:5:length(ydata)/10)
        end
        xlim([-0.01,1.01]*length(ydata))
        ax = gca;
        ax.YAxis.FontSize = 22;
        ax.XAxis.FontSize = 22;
        xlabel("time [s]","fontsize",25)
        ylabel("MUA [spike/s]","fontsize",25)
        saveas(gca, dpath + "STRF/strf_fit_"+scand(c)+"_"+b+".png")
    end
end

%% plot STRF kernel

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.45, 0.6]);

mstrf = squeeze(mean(strfs_.d200,3));
%mstrf = squeeze(mean(zscore(strfs_.d200,0,[1,2]),3));
imagesc(mstrf)
colorbar("Fontsize",20,"Ticks",-4:2:10)
xticks((1:10:numhist+1)*(numhist-1)/numhist)
xticklabels(-wsec/2*numhist:wsec/2*10:0)

ax = gca;
ax.YAxis.FontSize = 25;
ax.XAxis.FontSize = 25;

yticks(1:6:numband+1)
yticklabels(1*((2.^((0:6:numband+1)/6))))
xlabel("time [ms]","Fontsize",28)
ylabel("Frequency [kHz]","Fontsize",28)


saveas(gca, dpath + "STRF/strf_kernel.png")


%% plot R^2 


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.45, 0.6]);
tmp = [];
for i = 1:length(scand)
    tmp2 = [];
    for b = 1:div
        tmp2 = [tmp2; fitres.("f"+b+scand(i)).R2];
    end
    tmp = [tmp, tmp2];
end


h = boxplot(tmp);
hold on
set(h,'LineWidth',2)
for i = 1:length(scand)
    scatter(ones(div,1)*i+rand(div,1)*0.2-0.1, tmp(:,i),200, "g.")
end
xticks(1:length(scand))
xticklabels(scand + "%")

xlim([0.6,5.4])
ylim([0.1,0.62])
yticks(0.1:0.1:1)
ax = gca;
ax.YAxis.FontSize = 25;
ax.XAxis.FontSize = 25;
ylabel("R^2","Fontsize",28)

saveas(gca, dpath + "STRF/strf_R2.png")

%%

load('/Users/TP/research/Tlab/code/sound/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];

beatcont = zeros(length(scand), div, 5);
for c = 1:length(scand)
    trg = round(clspeed.("s"+scand(c))*1000/5) + 4;
    for b = 1:div
        pred = fun(fitres.("f"+b+scand(c)).x, xdata_.("x"+b+scand(c)));
        soundres = [];
        for t = 1:length(trg)
            soundres = [soundres, mean(pred(trg(t)+1:trg(t)+5)) - pred(trg(t))];
        end
        soundres(soundres < -1) = -1;
        for p = 1:5
            if p < 5
                beat = onbeat(:,p);
                beat(isnan(beat)) = [];
            else
                beat = allbeats;
            end
            offbeat = 1:489;
            offbeat(allbeats) = [];
            tmp = mean(soundres(offbeat))/mean(soundres(beat));
            disp([tmp,(1-tmp)/(1+tmp)])
            beatcont(c, b, p) = (1 - tmp)/(1 + tmp);
        end
    end
end


%% plot musical on-beat response of all position

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);

h = boxplot(beatcont(:,:,5)');
hold on
set(h,'LineWidth',2)
for i = 1:length(scand)
    scatter(ones(div,1)*i+rand(div,1)*0.2-0.1, beatcont(i,:,5),200, "g.")
end

xticks(1:length(scand))
xticklabels(scand + "%")
yticks(0:0.2:1)
ylim([0.1,0.9])

xlim([0.5,length(scand)+0.5])
ax = gca;
ax.YAxis.FontSize = 25;
ax.XAxis.FontSize = 25;

ylabel("Beat contrast", "FontSize",28)
%xlabel("Play speed (BPM)", "FontSize",28)
saveas(gca, dpath + "STRF/strf_beatcontrast.png")


%% plot musical beat contrast of each position (with color)

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);

clors = [0.8500,0.3250,0.0980;0,0.4470,0.7410;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560];
mrk = ["o-","^-","d-","v-"];
h1s = [];
for i = 1:4
    h1 = plot((1:length(scand))-0.25+0.1*i,mean(beatcont(:,:,i),2),mrk(i),"Color",clors(i,:),'LineWidth',2,"Markersize",20);
    set(h1, 'markerfacecolor', get(h1, 'color'))
    h1s = [h1s,h1];
    hold on
end

xticks(1:length(scand))
xticklabels(scand + "%")
yticks(0:0.2:1)

xlim([0.5,length(scand)+0.5])
ylim([0.1,0.9])
%legend(h1s,["pos.1","pos.2","pos.3","pos.4"], "location","south west","fontsize", 23);
[h,icons] = legend(h1s,["pos.1","pos.2","pos.3","pos.4"], "location","north west","fontsize", 23);
%icons = findobj(icons,'Type','line');
%icons = findobj(icons,'Marker','none','-xor');
%set(icons,'MarkerSize',15,'LineWidth',2);
ax = gca;
ax.YAxis.FontSize = 25;
ax.XAxis.FontSize = 25;

ylabel("Beat contrast", "fontsize", 28)
%xlabel("Play speed (BPM)", "FontSize",28)
saveas(gca, dpath + "STRF/beat_contrast_4.png")



