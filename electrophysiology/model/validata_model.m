%%
ipath = "/Users/yoshiki/research/Tlab/data/IPC/IPC/";
idr = dir(ipath + "21*.mat");


%% data processing
maxtime = 300;
meandata = struct();
times = struct();
isis = [32,56,100,178,316];

for y = 1:length(idr)
    load(ipath + idr(y).name)
    TC = ch_cortex;
    counter = [1,0,0,0,1];
    for x = 1:length(trials)
        cnt = find((isis - trials{x}.isi)==0);
        counter(cnt) = counter(cnt) + 1;
        if counter(cnt) ~= 2
            continue
        end
        click = find(trials{x}.click_or_silence == 1);
        time = trials{x}.trig_timings(click) / trials{x}.record_sF;
        time = time - time(1);
        click(time > maxtime) = [];
        time(time > maxtime) = [];
        times.("s"+trials{x}.isi) = time;
        data = trials{x}.spike_rates(TC, click);
        if isfield(meandata,"s"+trials{x}.isi)
            %meandata.("s"+trials{x}.isi) = [meandata.("s"+trials{x}.isi);mean(data)/mean(data(:,1))];
            meandata.("s"+trials{x}.isi) = [meandata.("s"+trials{x}.isi);mean(data)];
        else
            meandata.("s"+trials{x}.isi) = mean(data);
        end
    end
end


%% model calcuration of suggested model

ppath = "/Users/yoshiki/research/Tlab/data/yuta/";
load(ppath + "experiment_data.mat")
dtpkup = 1;
regstr = "";

r0s = [120:1:138];           %134,152,166,158,92
alphas = (5:1:10) + floor(pdata(1,1));         
betas = [0.042:0.001:0.047];   %0.037,0.042,0.03,0.035,0.051
inpx = 1000./[0.2,1,2,4,8,16,32];

z1 = 5; %5
z2 = 6; %6
z3 = 6; %6

criteria = [2,2,2,2,2];
Fs = 1000;
for r = 1:length(r0s)
    ratio = r0s(r);
    rsbs0 = struct();
    mse = [];
    for i = 1:length(isis)
        tt = times.("s"+isis(i))*Fs;
        rsb = ones(length(tt),1);
        lbound = 1 - finalv2(dtpkup)/ratio;
        X = alphas(z2);
        ys = X./pdata(:,1) - 1;
        st = betas(z3);
        inpy = calcurate(ys, st);
        for t = 1:(length(tt)-1)
            eks = 1:t;
            eks(tt(t+1) - tt(eks) > 5000) = [];
            kern = interp1(inpx, [0;inpy], tt(t+1) - tt(eks));
            I = dot(rsb(eks), kern);
            I = min([lbound,I]);
            rsb(t+1) = max([0,1 - I]);
        end
        rsbs0.("s"+isis(i)) = rsb*ratio;
        answer = mean(meandata.("s"+isis(i)))';
        tmp = mean((rsb*ratio - answer).^2)/var(answer);
        mse = [mse, tmp]; 
    end
    if mean(mse) < mean(criteria)
        rsbs = rsbs0;
        criteria = mse;
        ratios = ratio;
        disp([mean(criteria),ratio])
    end
end



%% model calcuration of Drew & Abbott (2006)

r0s_ = [140:1:170];   
alphas_ = [0.9:0.01:1.0];         
betas_ = 0.04:0.001:0.06;  

y1 = 7;
y2 = 3;
y3 = 15;

Fs = 1000;
criteria_ = [2,2,2,2,2]*0.6;
for r = 1:length(r0s_)
    rsb_s0 = struct();
    ratio_ = r0s_(r);
    mse = [];
    for i = 1:length(isis)
        tt = times.("s"+isis(i))*Fs;
        rsb_ = ones(length(tt),1);
        lbound = 1 - finalv2(dtpkup)/ratio_;
        alpha = alphas_(y2);
        beta = betas_(y3);
        for t = 1:(length(tt)-1)
            mn = 1;
            kernel = Integrator(tt(t+1)-tt(mn:t), beta);
            I = alpha*dot(rsb_(mn:t), kernel);
            I = min([lbound,I]);
            rsb_(t+1) = max([0,1 - I]);
        end
        rsb_s0.("s"+isis(i)) = rsb_*ratio_;
        answer = mean(meandata.("s"+isis(i)))';
        tmp = mean((rsb_*ratio_ - answer).^2)/var(answer);
        mse = [mse, tmp]; 
    end
    if mean(mse) < mean(criteria_)
        rsb_s = rsb_s0;
        ratios0 = ratio_;
        criteria_ = mse;
        disp([mean(criteria_),ratio_])
    end
end

%% model calcuration of Zuk

ppath = "/Users/yoshiki/research/Tlab/data/yuta/";
load(ppath + "experiment_data.mat")
dtpkup = 1;
regstr = "";

r0s = [105:1:115];           %134,152,166,158,92
alphaz = 0.48;

zkernelx = 0:200:800;
zkernely_ = [180, 134, 51, 11, 0]/180;

inpx = 1000./[0.2,1,2,4,8,16,32];

criteriaz = [2,2,2,2,2]*5;
Fs = 1000;
for r = 1:length(r0s)
    zkernely = zkernely_*alphaz;
    ratio = r0s(r);
    rsbz0 = struct();
    mse = [];
    for i = 1:length(isis)
        tt = times.("s"+isis(i))*Fs;
        rsb = ones(length(tt),1);
        lbound = 1 - finalv2(dtpkup)/ratio;
        for t = 1:(length(tt)-1)
            eks = 1:t;
            eks(tt(t+1) - tt(eks) > 800) = [];
            kern = interp1(zkernelx, zkernely, tt(t+1) - tt(eks));
            I = dot(rsb(eks), kern);
            I = min([lbound,I]);
            rsb(t+1) = max([0,1 - I]);
        end
        rsbz0.("s"+isis(i)) = rsb*ratio;
        answer = mean(meandata.("s"+isis(i)))';
        tmp = mean((rsb*ratio - answer).^2)/var(answer);
        mse = [mse, tmp]; 
    end
    if mean(mse) < mean(criteriaz)
        rsbz = rsbz0;
        criteriaz = mse;
        ratioz = ratio;
        disp([mean(criteriaz),ratio])
    end
end


%% plot suggested model

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 1]);

for i = 1:length(isis)
    subplot(length(isis),1,i)
    plot(times.("s"+isis(i)), mean(meandata.("s"+isis(i))))
    hold on
    plot(times.("s"+isis(i)), rsbs.("s"+isis(i)),"color",[1,0.2,0.2])
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    if i == 3
        ylabel("mean MUA","FontSize",23)
    end
    ylim([0,150])
    xlim([0,maxtime])
    xticks(0:60:maxtime)
    title(isis(i)*2 + " [ms]", "Fontsize",20)
end
xlabel("time [s]","FontSize",23)
legend(["Experiment","Fitting"],"FontSize",15,"position",[0.355,0.14,0,0])
saveas(gca, ipath + "figure/fitted.png")



%% plot Drew & Abbott 

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 1]);

for i = 1:length(isis)
    subplot(length(isis),1,i)
    plot(times.("s"+isis(i)), mean(meandata.("s"+isis(i))))
    hold on
    plot(times.("s"+isis(i)), rsb_s.("s"+isis(i)),"color",[0.6,0.3,0])
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    if i == 3
        ylabel("mean MUA","FontSize",23)
    end
    ylim([0,150])
    xlim([0,maxtime])
    xticks(0:60:maxtime)
    title(isis(i)*2 + " [ms]", "Fontsize",20)
end
xlabel("time [s]","FontSize",23)
legend(["Experiment","Fitting"],"FontSize",15,"position",[0.355,0.14,0,0])
%saveas(gca, ipath + "figure/fitted_abbott.png")


%% plot Zuk

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 1]);

for i = 1:length(isis)
    subplot(length(isis),1,i)
    plot(times.("s"+isis(i)), mean(meandata.("s"+isis(i))))
    hold on
    plot(times.("s"+isis(i)), rsbz.("s"+isis(i)),"color",[0.5,0.2,0.4])
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    if i == 3
        ylabel("mean MUA","FontSize",23)
    end
    ylim([0,150])
    xlim([0,maxtime])
    xticks(0:60:maxtime)
    title(isis(i)*2 + " [ms]", "Fontsize",20)
end
xlabel("time [s]","FontSize",23)
legend(["Experiment","Fitting"],"FontSize",15,"position",[0.355,0.14,0,0])
saveas(gca, ipath + "figure/fitted_zuk.png")

%%

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.55]);
plot(1-criteria,"o-","color",[1,0.2,0.2],"linewidth",2)
hold on
plot(1-criteria_,"o-","color",[0.6,0.3,0],"linewidth",2)
plot(1-criteriaz,"o-","color",[0.5,0.2,0.4],"linewidth",2)

plot([1,length(isis)],[0,0],"k--")
xticks(1:length(isis))
xticklabels(isis*2)
yticks(-1:0.5:0.5)
ylim([-1.1,0.5])
xlim([0.8, length(isis)+0.2])
ax = gca;
ax.YAxis.FontSize = 22;
ax.XAxis.FontSize = 22;

legend(["Present","Drew & Abbott(2006)", "Zuk et. al. (2018)"],"FontSize",20,"location","southwest")
%xlabel("dt [ms]","FontSize",25)
xlabel("mean ISI [ms]","FontSize",25)
ylabel("R^2","FontSize",25)
%saveas(gca, ipath + "figure/R2_error.png")

%% kernel plot (abbot vs model)

tts = 32:1030;

X = alphas(z2);
ys = X./pdata(:,1) - 1;
st = betas(z3);
inpx = 1000./[0.2,1,2,4,8,16,32];
inpy = calcurate(ys, st);

beta = betas_(y3);
kernel = Integrator(tts, beta);

zkernelx = [0,32,200,400,600,800,1200];
zkernely = [180, 172.64, 134, 51, 11, 0,0]/180 * 0.48;

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.55]);

plot(inpx, [0;inpy], "-","color",[1,0.2,0.2],"linewidth",4)
hold on
plot(tts, kernel, "-","color",[0.6,0.3,0],"linewidth",4)
plot(zkernelx(2:end), zkernely(2:end), "-","color",[0.5,0.2,0.4],"linewidth",4)

ylim([0,0.7])
xlim([0,1020])
yticks(0:0.2:1)
xlim([0,1020])
legend("Present model","Drew & Abbott (2006)", "Zuk et. al. (2018)","Fontsize",22)
ax = gca;
ax.YAxis.FontSize = 22;
ax.XAxis.FontSize = 22;

xlabel("t [ms]","fontsize",30)
%ylabel("Kernel function","fontsize",25)
ylabel("$K (t)$","fontsize",32, 'Interpreter','latex','fontweight','bold')

saveas(gca, ipath + "figure/kernel_compare3.png")


%%
function inpy = calcurate(ys_, st)
    inpx = 1000./[0.2,1,2,4,8,16,32];
    ys = ys_ * st / ys_(1);
    inpy = zeros(6,1);
    rps = zeros(4);
    inpy(1) = st;
    inpy(2) = ys(2) - ys(1);
    rps(1) = mean(inpy(1:2));
    inpy(3) = ys(3) - 3*rps(1);
    rps(2) = mean(inpy(2:3));
    inpy(4) = ys(4) - 5*rps(1) - 3*rps(2) + inpy(2);
    rps(3) = mean(inpy(3:4));
    inpy(5) = ys(5) - 9*rps(1) - 5*rps(2) - 3*rps(3) + inpy(2) + inpy(3);
    rps(4) = mean(inpy(4:5));
    inpy(6) = ys(6) - 17*rps(1) - 9*rps(2) - 5*rps(3) - 3*rps(4)+ inpy(2) + inpy(3)+inpy(4);
    [M,I] = max(inpy);
    inpy((I+1):6) = M;
    inpy(6) = inpy(5)/2;
    %[M,I] = min(inpy(1:4));
    %inpy(1:(I-1)) = M;
end

function kernel= Integrator(tt, beta)
     kernel = beta./(tt/1000+beta);
end

