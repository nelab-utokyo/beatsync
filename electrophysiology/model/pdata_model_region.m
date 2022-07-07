
%%%%%%%%%%%%%%
% new kernel
%%%%%%%%%%%%%%

load("experiment_data.mat")

dtpkup = 4;
if dtpkup == 1
    regstr = "";
elseif dtpkup == 2
    regstr = "a1";
elseif dtpkup == 3
    regstr = "aaf";
elseif dtpkup == 4
    regstr = "core";
elseif dtpkup == 5
    regstr = "belt";
end

if dtpkup > 1
    pdata = squeeze(rpdata(dtpkup-1,:,:));
    fdata = squeeze(rfdata(dtpkup-1,:,:));
    fratio = squeeze(rfratio(dtpkup-1,:,:));
    fcontrast = squeeze(rfcontrast(dtpkup-1,:,:));
end

if dtpkup == 1
    r0s = [128:2:138];           %134,152,166,158,92
    alphas = (5:1:10) + floor(pdata(1,1));         
    betas = [0.042:0.001:0.047];   %0.037,0.042,0.03,0.035,0.051
elseif dtpkup == 5
    r0s = [90:2:100];           %134,152,166,158,92
    alphas = (5:1:15) + floor(pdata(1,1));         
    betas = [0.055:0.001:0.065];   %0.037,0.042,0.03,0.035,0.051
else
    r0s = [148:2:153];           %134,152,166,158,92
    alphas = (7:1:12) + floor(pdata(1,1));         
    betas = [0.035:0.001:0.045];   %0.037,0.042,0.03,0.035,0.051
end

inpx = 1000./[0.2,1,2,4,8,16,32];

frqs = [1000, 500, 250, 125, 125/2, 125/4, 250, 125, 125/2, 125/4];
pfs = [1,1,1,1,1,1,2,2,2,2];
abbotp = zeros(sum(pfs==1),length(r0s),length(alphas),length(betas));
abbotf = zeros(sum(pfs==2),2,length(r0s),length(alphas),length(betas));

for i = 1:length(frqs)
    frq = frqs(i);
    Fs = 1000;
    time = 15;
    inp = zeros(time*Fs,1);
    if pfs(i) == 1
        inp(frq:frq:length(inp)) = 1;
    elseif pfs(i) == 2
        inp(frq:frq:length(inp)) = 1;
        temp = find(inp);
        inp(temp(4:4:end)) = 0;
    end
    tt = find(inp);
    for r0 = 1:length(r0s)
        rsb = ones(length(tt),length(alphas),length(betas));
        ratio = r0s(r0);
        lbound = 1 - finalv2(dtpkup)/ratio;
        for a = 1:length(alphas)
            X = alphas(a);
            ys = X./pdata(:,1) - 1;
            for b = 1:length(betas)
                st = betas(b);
                inpy = calcurate(ys, st);
                for t = 1:(length(tt)-1)
                    eks = 1:t;
                    eks(tt(t+1) - tt(eks) > 5000) = [];
                    %kern = interp1(inpx, alpha*[0;inpy], tt(t+1) - tt(eks), "pchip");
                    kern = interp1(inpx, [0;inpy], tt(t+1) - tt(eks));
                    I = dot(rsb(eks,a,b), kern);
                    I = min([lbound,I]);
                    %I = min([0.8,I]);
                    rsb(t+1,a,b) = max([0,1 - I]);
                end
            end
        end
        if pfs(i) == 1
            abbotp(i,r0,:,:) = mean(rsb)*ratio;   %pick up the mean value for idealres2
        elseif pfs(i) == 2
            abbotf(i-6,1,r0,:,:) = mean(rsb(1:3:end,:,:))*ratio;
            abbotf(i-6,2,r0,:,:) = (mean(rsb(2:3:end,:,:)) + mean(rsb(3:3:end,:,:)))*ratio/2;
        end
    end
end

tmp = squeeze(mean((abbotp - pdata(:,1)).^2));
tmpf = squeeze(mean(mean((abbotf -fdata(:,[1,3])).^2)));

tmpf_ = abbotf(:,1,:,:,:)./abbotf(:,2,:,:,:);
tmpf_ = squeeze(tmpf_);
tmpf_ = squeeze(mean((tmpf_ - fratio(:,1)).^2));


[M, I] = min(tmp(:));
[p1, p2, p3] = ind2sub(size(tmp),I);
%disp([r0s(p1),alphas(p2),betas(p3), M/(norm(pdata(:,1)).^2/6)])

[M, I] = min((tmp(:)/(norm(pdata(:,1)).^2/6) + tmpf(:)/(norm(fdata(:,[1,3])).^2/8))/2);
[f1, f2, f3] = ind2sub(size(tmpf),I);
%disp([r0s(f1),alphas(f2),betas(f3), M*1000])

%[M, I] = min((tmp(:)/(norm(pdata(:,1)).^2/6) + tmpf(:)/(norm(fdata(:,[1,3])).^2/8) ...
%           + tmpf_(:)/(norm(fratio(:,1)).^2/4))/3);
[M, I] = min((tmpf(:)/(norm(fdata(:,[1,3])).^2/8)));
[m1, m2, m3] = ind2sub(size(tmpf),I);
disp([r0s(m1),alphas(m2),betas(m3),M*1000])

[M, I] = min(((tmpf(:)/(norm(fdata(:,[1,3])).^2/8))+tmpf_(:)/(norm(fratio(:,1)).^2/4))/2);
[m1_, m2_, m3_] = ind2sub(size(tmpf),I);
%disp([r0s(m1_),alphas(m2_),betas(m3_),M*1000])

z1 = m1; 
z2 = m2; 
z3 = m3; 


%% check the model accuracy

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.32, 1]);

subplot(2,1,1)
labels = ["1/1","1/2","1/4","1/8","1/16","1/32"];
frqs = [1000, 500, 250, 125, 125/2, 125/4];
h1 = errorbar(log2(1000./frqs), pdata(:,1), pdata(:,2), "color",[0.2,0.6,0.2],'LineWidth',4);
hold on
h2 = plot(log2(1000./frqs), abbotp(:,z1,z2,z3),"-.k" ,'LineWidth',3);
title("Periodic stim.","Fontsize",27)
xlim([-0.2,5.2])
ylim([20,130])
if dtpkup == 4
    ylim([20,155])
elseif dtpkup == 5
    ylim([15,90])
end
yticks(0:30:150)
xticks(0:length(frqs)-1)
xticklabels(labels)
ax = gca;
ax.YAxis.FontSize = 22;
ax.XAxis.FontSize = 22;
ylabel("mean MUA","FontSize",25)
xlabel("ISI [s]","Fontsize",25)
legend([h1,h2],["Data","Model"],"FontSize",20)

subplot(2,1,2)

labels = [60,120,240,480];

h1 = errorbar(1:size(fratio,1), fdata(:,1), fdata(:,2), "-","color",[0.7,0.8,0.3],'LineWidth',4);
hold on
h2 = errorbar(1:size(fratio,1), fdata(:,3), fdata(:,4),"-","color",[0.3,0.4,0], 'LineWidth',4);

h3 = plot(1:size(fratio,1), abbotf(:,1,z1,z2,z3),"k-.",'LineWidth',3);
plot(1:size(fratio,1), abbotf(:,2,z1,z2,z3),"k-.",'LineWidth',3)
legend([h1,h2,h3], ["On beat","Off beat","Model"],"FontSize",20)
title("Rhythmic stim.","FontSize",27)
xlim([0.8,4.2])
xticks(1:size(fratio,1))
xticklabels(labels)
yticks(30:30:150)
ylim([20,130])
if dtpkup == 4
    ylim([20,155])
elseif dtpkup == 5
    ylim([15,90])
end
ax = gca;
ax.YAxis.FontSize = 22;
ax.XAxis.FontSize = 22;
ylabel("mean MUA","FontSize",25)
xlabel("BPM","Fontsize",25)

saveas(gca,"./summary/compare_data_pdata_"+regstr+".fig")
saveas(gca,"./summary/compare_data_pdata_"+regstr+".png")

%% plot MUA beat contrast and kernel function

labels = [60,120,240,480];
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.48]);
if dtpkup > 3
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.32, 1]);
end
subplot(2,1,2)
inpx = 1000./[0.2,1,2,4,8,16,32];
clrs = [1,0.5,0;1,0.3,0.3;1,0.3,0.3;1,0.3,0.3;0.3,0.3,1];
tts = 32:1030;

isis = 1./[1,2,4,8,16];
for i = 1:length(isis)
    plot([1000,1000]*isis(i),[0,1],"k--","Linewidth",1.2)
    hold on
end

X = alphas(z2);
st = betas(z3);
inpy = calcurate(X./pdata(:,1)-1,st);
inpy(end) = 0;
inpx(end) = 0;
h3 = plot(inpx, interp1(inpx, [0;inpy], inpx), "o-","color",clrs(dtpkup,:), 'LineWidth',4);

ylim([0,0.6])
xlim([0,1020])
yticks(0:0.2:1)
xticks(0:250:1000)
xlim([0,1020])
%legend([h1,h2],"Drew & Abbott (2006)","Suggested","Fontsize",22)
ax = gca;
ax.YAxis.FontSize = 22;
ax.XAxis.FontSize = 22;

%xlabel("t [ms]","fontsize",40,'Interpreter','latex','fontweight','bold')
ylabel("$K (t)$","fontsize",30,'Interpreter','latex','fontweight','bold')
xlabel("t [ms]","fontsize",25)
%ylabel("Kernel function","fontsize",25)


subplot(2,1,1)
predicted = pdata(2:5,1)./pdata(3:6,1);

h1 = errorbar(1:length(labels),fcontrast(:,1),fcontrast(:,2),"-g", 'LineWidth',4);
hold on
model = (abbotf(:,1,z1,z2,z3)-abbotf(:,2,z1,z2,z3))./(abbotf(:,1,z1,z2,z3)+abbotf(:,2,z1,z2,z3));
h2 = plot(model, "-.k", 'LineWidth',3);

xlim([0.8,4.2])
xticks(1:length(labels))
xticklabels(labels)
yticks(0:0.1:0.5)
ylim([0,0.23])
legend([h1,h2],["Data","Model"],"FontSize",20)
if dtpkup == 4
    %ylim([0.93,1.65])
    ylim([0,0.3])
    legend([h1,h2],["Data","Model"],"FontSize",20,"location","northeast")
elseif dtpkup == 5
    %ylim([0.85,1.5])
    ylim([0,0.3])
end
ax = gca;
ax.YAxis.FontSize = 22;
ax.XAxis.FontSize = 22;
ylabel("Neural beat contrast", "fontsize",25)
xlabel("BPM","Fontsize",25)

saveas(gca,"./summary/model_kernel_"+regstr+"2.fig")
saveas(gca,"./summary/model_kernel_"+regstr+"2.png")


%% check kernel function

figure()
inpx = 1000./[0.2,1,2,4,8,16,32];
clrs = [1,0.5,0;1,0.3,0.3;1,0.3,0.3;1,0.3,0.3;0.3,0.3,1];
tts = 0:1030;

%figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.5]);

isis = 1./[1,2,4,8,16,32];
for i = 1:length(isis)
    plot([1000,1000]*isis(i),[0,1],"k--","Linewidth",1.2)
    hold on
end

X = alphas(z2);
st = betas(z3);
inpy = calcurate(X./pdata(:,1)-1,st);
h3 = plot(tts, interp1([inpx,0], [0;inpy;0], tts), "color",clrs(dtpkup,:), 'LineWidth',4);

ylim([0,0.6])
xlim([0,1020])
yticks(0:0.2:1)
xticks(0:250:1000)

xlim([0,1020])
%legend([h1,h2],"Drew & Abbott (2006)","Suggested","Fontsize",22)
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;

xlabel("t [ms]","fontsize",40,'Interpreter','latex','fontweight','bold')
%ylabel("Kernel function","fontsize",25)
ylabel("$K (t)$","fontsize",40,'Interpreter','latex','fontweight','bold')



%% plot the time course of suggested model (all auditory neurons only)

%load("experiment_data.mat")

frqs = [1000, 500, 250, 125, 125/2, 125/4, 250, 125, 125/2, 125/4];
pfs = [1,1,1,1,1,1,2,2,2,2];
pfnames = ["p1","p2","p4","p8","p16","p32","f4","f8","f16","f32"];

ratio = r0s(z1);
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 1]);
labels = ["60","120","240","480"] + " BPM ";

for i = 7:10
    if pfs(i) == 1
        target = evpsum.(pfnames(i));
    else
        target = evfsum.(pfnames(i));
    end
    frq = frqs(i);
    Fs = 1000;
    time = 15;
    inp = zeros(time*Fs,1);

    if pfs(i) == 1
        inp(round(frq:frq:length(inp))) = 1;
    elseif pfs(i) == 2
        inp(round(frq:frq:length(inp))) = 1;
        temp = find(inp);
        inp(temp(4:4:end)) = 0;
    end
    tt = find(inp);
    stms = ones(length(tt),1);
    rsb = ones(length(tt),1);
    X = alphas(z2);
    st = betas(z3);
    ys = X./pdata(:,1) - 1;
    inpy = calcurate(ys,st);
    for t = 1:(length(tt)-1)
        eks = 1:t;
        eks(tt(t+1) - tt(eks) > 5000) = [];
        kern = interp1(inpx, [0;inpy], tt(t+1) - tt(eks));
        I = dot(rsb(eks,1), kern);
        I = min([1-finalv2(dtpkup)/ratio,I]);
        %I = min([0.8,I]);
        rsb(t+1,1) = max([0,1 - I]);
    end
    rsb = rsb * ratio;
    
    subplot(4,1,i-6)
    %plot(target,"r", "Linewidth",2)
    plot(target,"g", "Linewidth",2)
    hold on
    plot(rsb,"k", "Linewidth",1.2)
    ylim([0,160])
    xlim([0,size(rsb,1)])
    xticks(0:(length(rsb)/3):length(rsb))
    xticklabels([0,5,10,15])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel("mean MUA", "Fontsize",25)
    if i-6 == 1
        %title("Suggested model", "Fontsize",25)
        legend("Data","Model","Fontsize",23, "Location","Southwest")
    end
    text(size(rsb,1),140,labels(i-6),'FontSize',25,"HorizontalAlignment","right")
end

xlabel("Time [s]", "Fontsize",25)

%saveas(gca,"./summary/pdatamodel_fit_f.fig")
%saveas(gca,"./summary/pdatamodel_fit_f.png")

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


function inpy = calcurate2(ys_, st)
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
    %[M,I] = max(inpy);
    %inpy((I+1):6) = M;
    %inpy(6) = inpy(5)/2;
    %[M,I] = min(inpy(1:4));
    %inpy(1:(I-1)) = M;
end

function kernel= Integrator(tt, beta)
     kernel = beta./(tt/1000+beta);
end
