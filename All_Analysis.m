clear all
clc


% ---- Load Stimulation outputs
% Here


params = [];
N = 400;

t_burn = 1500;
cort_phi_e = cort_phi_e(t_burn:end,:,:); % Remove burn-in



% --- Haemodynamic Modeling

TR = 0.586;
fNY=1/(2*dt);
[Bfilt,Afilt] = butter(3,0.1/fNY,'low');

[hrf,~] = spm_hrf(TR);
bold_ts = [];

for pp =1: size(cort_phi_e,3)
    for rr = 1: size(cort_phi_e,2)
        x = zscore(cort_phi_e(:,rr,pp));
        conv_ts = conv(hrf,x);
        xx = conv_ts(1:end-length(hrf)+1,1);
        bold_ts(:,rr,pp) = xx;
    end
end


%% HRF Correlations
[hrf,~] = spm_hrf(0.586);

J1 = corr(bold_ts(:,:,1));
J2 = corr(bold_ts(:,:,2));
J3 = corr(bold_ts(:,:,3));
J4 = corr(bold_ts(:,:,4));

J1(logical(eye(400))) = 0;
J2(logical(eye(400))) = 0;
J3(logical(eye(400))) = 0;
J4(logical(eye(400))) = 0;

max_corr = max([max(J1,[],'all'),max(J2,[],'all'),max(J3,[],'all'),max(J4,[],'all')]);
min_corr = min([min(J1,[],'all'),min(J2,[],'all'),min(J3,[],'all'),min(J4,[],'all')]);

figure
subplot(141)
imagesc(J1)
caxis([min_corr max_corr])
xticks = [];
yticks = [];

subplot(142)
figure
imagesc(J2)
caxis([min_corr max_corr])
xticks = [];
yticks = [];

subplot(143)
figure
imagesc(J3)
caxis([min_corr max_corr])
xticks = [];
yticks = [];

subplot(144)
figure
imagesc(J4)
caxis([min_corr max_corr])
xticks = [];
yticks = [];

%% MSD Landscape

ndt = 100; % Number of lags
ds = 0:1:20; % Number of MSD divisions

nrgSig = nan(ndt,numel(ds));
nrgBase = nan(ndt,numel(ds));

for mm = 1:size(bold_ts,3)

    for tt = 1:ndt
    
        % Normalize data
        cortSig = squeeze(zscore(bold_ts(:,:,mm)));

        % MSD calculation
        MSD = mean((cortSig(1+tt:end,:) - cortSig(1:end-tt,:)).^2,2);
            
        % Calculate probability distribution and energy for each dt
        [nrgSigdt] = PdistGaussKern(MSD,ds);
            
        % Pool results across time
        nrgSig(tt,:,mm) = nrgSigdt;

    end
end

%
%  Plot estimated MSD energy landscape 

x = [1:ndt].*dt;
y = ds;
[X,Y] = meshgrid(x,y);

xmin = x(1);
xmax = x(end);

figure
subplot(2,4,1)
mesh(X,Y,squeeze(nrgSig(:,:,2))','EdgeColor', 'k')
xlabel('Time [s] ')
ylabel('MSD')
zlabel('MSD  energy')
xlim([xmin xmax])
%zlim([2 78])
ylim([1 ds(end)])
view(40,26)   % XZ
title('Propofol')

subplot(2,4,2)
mesh(X,Y,squeeze(nrgSig(:,:,4))','EdgeColor', 'm')
xlabel('Time [s]')
ylabel('MSD')
zlabel('MSD  energy')
xlim([xmin xmax])
%zlim([2 78])
ylim([1 ds(end)])
view(40,26)   % XZ
title('Low Matrix')

subplot(2,4,3)
mesh(X,Y,squeeze(nrgSig(:,:,3))','EdgeColor', 'r')
xlabel('Time [s]')
ylabel('MSD')
zlabel('MSD  energy')
xlim([xmin xmax])
%zlim([2 78])
ylim([1 ds(end)])
view(40,26)   % XZ
title('High Matrix')

subplot(2,4,4)
mesh(X,Y,squeeze(nrgSig(:,:,1))','EdgeColor', [0.4660 0.6740 0.1880])
xlabel('Time [s]')
ylabel('MSD')
zlabel('MSD  energy')
xlim([xmin xmax])
%zlim([2 78])
ylim([1 ds(end)])
view(40,26)   % XZ
title('Wake')


x_slice = [find(x >= 0.125,1),find(x >= 0.155,1)];

ymax = 100;

subplot(245)
hold on
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,2),1)),'LineWidth',1.5,'Color','k')
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,1),1)),'--','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880])
ylim([0 ymax])
xlabel('MSD')
ylabel('MSD energy')

subplot(246)
hold on
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,4),1)),'LineWidth',1.5,'Color','m')
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,1),1)),'--','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880])
ylim([0 ymax])
box on
xlabel('MSD')
ylabel('MSD energy')

subplot(247)
hold on
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,3),1)),'-b','LineWidth',1.5)
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,2),1)),'--','LineWidth',1.5,'Color','k')

ylim([0 ymax])
xlabel('MSD')
ylabel('MSD energy')

subplot(248)
hold on
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,1),1)),'-','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880])
plot(y, squeeze(mean(nrgSig(x_slice(1):x_slice(2),:,2),1)),'--','LineWidth',1.5,'Color','k')
ylim([0 ymax])
box on
xlabel('MSD')
ylabel('MSD energy')


%% Coherence LIP-FEF / Wake vs Propofol

lip_ids = [124:126,136:138];
fef_ids = [127,130:132];

ts1 = cort_phi_e(:,:,1);
ts2 = cort_phi_e(:,:,2);

coh_all1 = [];
coh_all2 = [];

for xx = 1:length(lip_ids)
    for yy = 1:length(fef_ids)


        x = ts1(:,lip_ids(xx));
        y = ts1(:,fef_ids(yy));

        n_windows = 50;
        frac_overlap = 0.5;
        window_idx = nf.partition(size(x,1),n_windows,[],frac_overlap,1,1);
        K = length(window_idx);

        Npts = window_idx(1,2) - window_idx(1,1);

        time = (1:size(x,1))/fs;
        T = time(window_idx(1,2))-time(1);

        Sxx = zeros(K,Npts);
        Syy = zeros(K,Npts);
        Sxy = zeros(K,Npts);


        for k=1:K

            x_w = x(window_idx(k,1):window_idx(k,2)-1);
            y_w = y(window_idx(k,1):window_idx(k,2)-1);


            Sxx(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(x_w));
            Syy(k,:) = 2*dt^2/T*fft(y_w).*conj(fft(y_w));
            Sxy(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(y_w));

        end

        Sxx = Sxx(:,1:floor(Npts/2)+1);
        Syy = Syy(:,1:floor(Npts/2)+1);
        Sxy = Sxy(:,1:floor(Npts/2)+1);


        Sxx = mean(Sxx,1);
        Syy = mean(Syy,1);
        Sxy = mean(Sxy,1);

        coh1 = abs(Sxy).^2./(Sxx.*Syy);

        df = 1/T;

        fNQ = 1/dt/2;
        coh_f = [0:df:fNQ];

        coh1 = real(coh1);


        % Propofol
        x = ts2(:,lip_ids(xx));
        y = ts2(:,fef_ids(yy));

        %n_windows = 20;
        frac_overlap = 0.5;
        window_idx = nf.partition(size(x,1),n_windows,[],frac_overlap,1,1);
        K = length(window_idx);

        Npts = window_idx(1,2) - window_idx(1,1);

        time = (1:size(x,1))/fs;
        T = time(window_idx(1,2))-time(1);

        Sxx = zeros(K,Npts);
        Syy = zeros(K,Npts);
        Sxy = zeros(K,Npts);

        for k=1:K

            x_w = x(window_idx(k,1):window_idx(k,2)-1);
            y_w = y(window_idx(k,1):window_idx(k,2)-1);


            Sxx(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(x_w));
            Syy(k,:) = 2*dt^2/T*fft(y_w).*conj(fft(y_w));
            Sxy(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(y_w));

        end

        Sxx = Sxx(:,1:floor(Npts/2)+1);
        Syy = Syy(:,1:floor(Npts/2)+1);
        Sxy = Sxy(:,1:floor(Npts/2)+1);


        Sxx = mean(Sxx,1);
        Syy = mean(Syy,1);
        Sxy = mean(Sxy,1);

        coh2 = abs(Sxy).^2./(Sxx.*Syy);

        df = 1/T;

        fNQ = 1/dt/2;
        coh_f = [0:df:fNQ];

        coh2 = real(coh2);

        %coh_diff1 = [coh_diff1; coh1 - coh2];
        coh_all1 = [coh_all1; coh1];
        coh_all2 = [coh_all2; coh2];



    end
    disp(xx)
end


f_min = find(coh_f > 1, 1,'first');
f_min = 1;

figure
hold on

shadedErrorBar(coh_f(f_min:end),mean(coh_all2(:,f_min:end),1),std(coh_all2(:,f_min:end),[],1)./sqrt(size(coh_all2,1)),'lineprops',{'-r'})
shadedErrorBar(coh_f(f_min:end),mean(coh_all1(:,f_min:end),1),std(coh_all1(:,f_min:end),[],1)./sqrt(size(coh_all1,1)),'lineprops',{'-k'})
xlim([0 90])
ylim([0 0.25])
box on

legend('Propofol','Wake')
set(gca,'Fontsize',15)

set(0,'defaultTextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');


%% Coherence LIP-FEF / Matrix Strong vs Weak

lip_ids = [124:126,136:138];
fef_ids = [127,130:132];

ts1 = cort_phi_e(:,:,3);
ts2 = cort_phi_e(:,:,4);

d_50 = designfilt('bandstopiir','FilterOrder',2, ...
       'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
       'DesignMethod','butter','SampleRate',fs);

%ts1 = filtfilt(d_50,ts1);
%ts2 = filtfilt(d_50,ts2);

%coh_diff1 = [];
coh_all3 = [];
coh_all4 = [];

for xx = 1:length(lip_ids)
    for yy = 1:length(fef_ids)


        x = ts1(:,lip_ids(xx));
        y = ts1(:,fef_ids(yy));

        n_windows = 50;
        frac_overlap = 0.5;
        window_idx = nf.partition(size(x,1),n_windows,[],frac_overlap,1,1);
        K = length(window_idx);

        Npts = window_idx(1,2) - window_idx(1,1);

        time = (1:size(x,1))/fs;
        T = time(window_idx(1,2))-time(1);

        Sxx = zeros(K,Npts);
        Syy = zeros(K,Npts);
        Sxy = zeros(K,Npts);


        for k=1:K

            x_w = x(window_idx(k,1):window_idx(k,2)-1);
            y_w = y(window_idx(k,1):window_idx(k,2)-1);


            Sxx(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(x_w));
            Syy(k,:) = 2*dt^2/T*fft(y_w).*conj(fft(y_w));
            Sxy(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(y_w));

        end

        Sxx = Sxx(:,1:floor(Npts/2)+1);
        Syy = Syy(:,1:floor(Npts/2)+1);
        Sxy = Sxy(:,1:floor(Npts/2)+1);


        Sxx = mean(Sxx,1);
        Syy = mean(Syy,1);
        Sxy = mean(Sxy,1);

        coh1 = abs(Sxy).^2./(Sxx.*Syy);

        df = 1/T;

        fNQ = 1/dt/2;
        coh_f = [0:df:fNQ];

        coh1 = real(coh1);


        % Propofol
        x = ts2(:,lip_ids(xx));
        y = ts2(:,fef_ids(yy));

        %n_windows = 20;
        frac_overlap = 0.5;
        window_idx = nf.partition(size(x,1),n_windows,[],frac_overlap,1,1);
        K = length(window_idx);

        Npts = window_idx(1,2) - window_idx(1,1);

        time = (1:size(x,1))/fs;
        T = time(window_idx(1,2))-time(1);

        Sxx = zeros(K,Npts);
        Syy = zeros(K,Npts);
        Sxy = zeros(K,Npts);

        for k=1:K

            x_w = x(window_idx(k,1):window_idx(k,2)-1);
            y_w = y(window_idx(k,1):window_idx(k,2)-1);


            Sxx(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(x_w));
            Syy(k,:) = 2*dt^2/T*fft(y_w).*conj(fft(y_w));
            Sxy(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(y_w));

        end

        Sxx = Sxx(:,1:floor(Npts/2)+1);
        Syy = Syy(:,1:floor(Npts/2)+1);
        Sxy = Sxy(:,1:floor(Npts/2)+1);


        Sxx = mean(Sxx,1);
        Syy = mean(Syy,1);
        Sxy = mean(Sxy,1);

        coh2 = abs(Sxy).^2./(Sxx.*Syy);

        df = 1/T;

        fNQ = 1/dt/2;
        coh_f = [0:df:fNQ];

        coh2 = real(coh2);

        %coh_diff1 = [coh_diff1; coh1 - coh2];
        coh_all3 = [coh_all3; coh1];
        coh_all4 = [coh_all4; coh2];



    end
    disp(xx)
end


f_min = find(coh_f > 1, 1,'first');
f_min = 1;

figure
hold on

shadedErrorBar(coh_f(f_min:end),mean(coh_all4(:,f_min:end),1),std(coh_all4(:,f_min:end),[],1)./sqrt(size(coh_all4,1)),'lineprops',{'-m'})
shadedErrorBar(coh_f(f_min:end),mean(coh_all3(:,f_min:end),1),std(coh_all3(:,f_min:end),[],1)./sqrt(size(coh_all3,1)),'lineprops',{'-b'})
xlim([0 40])
%ylim([0 0.2])
box on

legend('Low Matrix','High Matrix')
set(gca,'Fontsize',15)

set(0,'defaultTextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');


%% CAD measures


bsv = [];
sus = [];
exp_var = [];
v_stab = [];
m_stab = [];
time_sc = [];
alpha_lip_fef_coh = [];
part_coeff_7 = [];
part_coeff_17 = [];
kc = [];

lip_ids = [124:126,136:138];
fef_ids = [127,130:132];

for pp = 1: size(cort_phi_e,3)

    % -------------------- Functional Connectivity  --------------
    fC = corr(bold_ts(:,:,pp));
    fC(logical(eye(N))) = 0;
    

    % -------------------- Regional diversity  --------------
    % Grab the upper triangle of the regional FC
    up_idx  = ones(size(fC));
    up_idx = triu(up_idx,1);
    state_region = fC(logical(up_idx));
    
    % Standard deviation of the regional FC
    reg_div(pp) = std(state_region(:));

    % -------------------- Susceptibility  --------------
    data = cort_phi_e(:,:,pp);
    [Z,~,~] = zscore(data,0,1);
    density = sum(Z > 0,2)/N; % Number of nodes above mean
    
    sus(pp) = (mean(density.^2) - mean(density)^2)/mean(density);

    % -------------------- Participation --------------------
    load('schaef_id.mat')
    part_coeff_7(pp) = mean(participation_coef_sign(fC,schaef_id));


    % -------------------- Timescale  --------------
    acf = [];
    for rr = 1: N
        [acf(:,rr),lags] = autocorr(cort_phi_e(:,rr,pp));
    end
    avg_acf(:,pp) = mean(acf,2);
    grad_acf = gradient(avg_acf(:,pp));
    x_max = find(grad_acf >= 0,1,'first');

    x = lags(1:x_max);
    y = avg_acf(1:x_max,pp);
    g = fittype('a-b*exp(-c*x)');
    f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
    fitvalues = coeffvalues(f0);
    time_sc(pp) = fitvalues(3); % save c parameter

    % -------------------- BOLD Timescale  --------------
    acf = [];
    for rr = 1: N
        [acf(:,rr),lags] = autocorr(zscore(bold_ts(:,rr,pp)));
    end
    avg_acf(:,pp) = mean(acf,2);
    grad_acf = gradient(avg_acf(:,pp));
    x_max = find(grad_acf >= 0,1,'first');

    x = lags(1:x_max);
    y = avg_acf(1:x_max,pp);
    g = fittype('a-b*exp(-c*x)');
    f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
    fitvalues = coeffvalues(f0);
    bold_time_sc(pp) = fitvalues(3); % save c parameter


    % -------------------- Metastability  --------------
    phi = [];
    for rr = 1: N
        x = bandpass(downsample(bold_ts(:,rr,pp),0.01/dt),[0.01, 0.1]);
        X = hilbert(x);
        phi(:,rr) = angle(X);
    end
    m_stab(pp) = var(abs(mean(exp(1i.*phi),2))); % regional average

    % -------------------- Kolmogorov Complexity -------------------
    temp_kc = [];
    for rr = 1: N
        temp_kc(rr) = kolmogorov(zscore(cort_phi_e(:,rr,pp)) > 0);
    end
    kc(pp) = mean(temp_kc);

    % -------------------- Coherence -------------------

    coh_all1 = [];
    for xx = 1:length(lip_ids)
        for yy = 1:length(fef_ids)
    
    
            x = cort_phi_e(:,lip_ids(xx),pp);
            y = cort_phi_e(:,fef_ids(yy),pp);
    
            n_windows = 50;
            frac_overlap = 0.5;
            window_idx = nf.partition(size(x,1),n_windows,[],frac_overlap,1,1);
            K = length(window_idx);
    
            Npts = window_idx(1,2) - window_idx(1,1);
    
            time = (1:size(x,1))/fs;
            T = time(window_idx(1,2))-time(1);
    
            Sxx = zeros(K,Npts);
            Syy = zeros(K,Npts);
            Sxy = zeros(K,Npts);
    
    
            for k=1:K
    
                x_w = x(window_idx(k,1):window_idx(k,2)-1);
                y_w = y(window_idx(k,1):window_idx(k,2)-1);
    
    
                Sxx(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(x_w));
                Syy(k,:) = 2*dt^2/T*fft(y_w).*conj(fft(y_w));
                Sxy(k,:) = 2*dt^2/T*fft(x_w).*conj(fft(y_w));
    
            end
    
            Sxx = Sxx(:,1:floor(Npts/2)+1);
            Syy = Syy(:,1:floor(Npts/2)+1);
            Sxy = Sxy(:,1:floor(Npts/2)+1);
    
    
            Sxx = mean(Sxx,1);
            Syy = mean(Syy,1);
            Sxy = mean(Sxy,1);
    
            coh1 = abs(Sxy).^2./(Sxx.*Syy);
    
            df = 1/T;
    
            fNQ = 1/dt/2;
            coh_f = [0:df:fNQ];
    
            coh1 = real(coh1);
            coh_all1 = [coh_all1; coh1];
            
        end
    end
    avg_coh = mean(coh_all1,1);

    f_min = find(coh_f > 8,1,'first');
    f_max = find(coh_f > 13,1,'first');

    alpha_lip_fef_coh(pp) = sum(avg_coh(1,f_min:f_max),2);

    disp(pp)
end



% Spider-Plot
X = [reg_div.*1e3; sus.*1e3; part_coeff_17; time_sc; 1./bold_time_sc; m_stab.*1e3; kc; alpha_lip_fef_coh];
X1 = X(:,1:2);
X2 = X(:,3:4);
tick_data = {'','','',''};
figure
spider_plot(X.','AxesLabels', {'FC Variance [$10^{-3}$]', 'Susceptibility [$10^{-3}$]', 'Network Participation', 'Time Scale','BOLD Time Scale','Synchronization Variance [$10^{-3}$]','Kolmogorov Complexity','8-13Hz LIP-FEF Coherence'},'AxesInterpreter','latex','FillOption', 'on')
legend('Wake','Propofol','High Matrix','Low Matrix')


lims = [50, 5, 0.91, 0.6, 0.085, 10, 0.94, 0.4; 70, 11, 0.93, 1.4, 0.12, 18, 0.99, 1.3];
ax_prec = [2,2,2,2,2,2,2,2];

figure
spider_plot(X.','AxesLabels', {'FC Variance [$10^{-3}$]', 'Susceptibility [$10^{-3}$]', 'Network Participation', 'Time Scale','BOLD Time Scale','Synchronization Variance [$10^{-3}$]','Kolmogorov Complexity','8-13Hz LIP-FEF Coherence'},'AxesLimits', lims,'AxesPrecision', ax_prec,...
    'AxesInterpreter','latex','FillOption', 'on','AxesTickLabels',tick_data);

legend('Wake','Propofol','High Matrix','Low Matrix')
%legend('High Matrix','Low Matrix')


%% Transfer Entropy - Gaussian Estimator
clc
% Change location of jar to match yours:
javaaddpath('infoDynamics.jar');
clear teCalc


% ---- BOLD-transformed Time series
data = [];
for pp =1: 4
    for rr = 1: 400
        data(:,rr,pp) = downsample(bold_ts(:,rr,pp),0.01/dt);
    end
end


xx_ids = [1:400];
yy_ids = [1:400];

% ----------- Set ACF Timescale
acf = [];
for rr = 1: N
    [acf(:,rr),lags] = autocorr(data(:,rr,pp));
end
avg_acf(:,pp) = mean(acf,2);
grad_acf = gradient(avg_acf(:,pp));
x_max = find(grad_acf >= 0,1,'first');


% ------------------------------------------

acf_length = x_max;% Set ACF Timescale

% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian');
teCalc.setProperty('NORMALISE', 'true'); % Normalise the individual variables
teCalc.setProperty('DYN_CORR_EXCL', string(acf_length)); % 

TE_mtx = zeros(length(yy_ids), length(xx_ids), 4);
for mm = 1: 4
   for yy = 1:length(xx_ids)
        teCalc.setProperty('AUTO_EMBED_METHOD', 'MAX_CORR_AIS_DEST_ONLY'); % AutoEmbedding target
        teCalc.setProperty('AUTO_EMBED_K_SEARCH_MAX', '10'); % AutoEmbedding target
        teCalc.setProperty('AUTO_EMBED_TAU_SEARCH_MAX', '2'); % AutoEmbedding target

        k = 1;
        k_tau = 1;

        for xx = 1: length(yy_ids)
            if xx  > 1
               teCalc.setProperty('AUTO_EMBED_METHOD', 'NONE'); % AutoEmbedding target
               teCalc.setProperty('k_HISTORY', string(k)); % AutoEmbedding target
               teCalc.setProperty('k_TAU', string(k_tau)); % AutoEmbedding target
            end

            teCalc.initialise(); % Use history length 1 (Schreiber k=1)
            teCalc.setObservations(squeeze(data(:,xx_ids(xx),mm)), squeeze(data(:,yy_ids(yy),mm)));
            % For copied source, should give something close to 1 bit:
            TE_mtx(xx,yy,mm) = teCalc.computeAverageLocalOfObservations();

            if xx == 1
                k = str2num(teCalc.getProperty('k_HISTORY')); % AutoEmbedding target - TO DO: Store
                k_tau = str2num(teCalc.getProperty('k_TAU')); % AutoEmbedding target
            end
        
        end
        disp(yy)
    end
end

%% Active Information Storage - Gaussian Estimator
clc
javaaddpath('infoDynamics.jar');

xx_ids = [1:400];
yy_ids = [1:400];

% Create a TE calculator and run it:
infoStor=javaObject('infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorGaussian');
infoStor.setProperty('NORMALISE', 'true'); % Normalise the individual variables
infoStor.setProperty('DYN_CORR_EXCL', string(acf_length)); % Set AR Timescale - TO DO: Check AR Timescale

IA_mtx = zeros(length(xx_ids), 4);
for mm = 1: 4
    for xx = 1: length(xx_ids)

            infoStor.setProperty('AUTO_EMBED_METHOD', 'MAX_CORR_AIS'); % AutoEmbedding target
            infoStor.setProperty('AUTO_EMBED_K_SEARCH_MAX', '10'); % AutoEmbedding target
            infoStor.setProperty('AUTO_EMBED_TAU_SEARCH_MAX', '2'); % AutoEmbedding target

            infoStor.initialise(); % Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised unit
        
            infoStor.setObservations((squeeze(cort_phi_e(:,xx_ids(xx),mm))));
            % For copied source, should give something close to 1 bit:
            IA_mtx(xx,mm) = infoStor.computeAverageLocalOfObservations();
        
        disp(xx)
    end
end

%% Plot TE and AIS
nbins = 20;
edges = [0:0.3/nbins: 0.3];

for pp =1:4
    figure
    histogram(sum(TE_mtx(:,:,pp),2),edges,'Normalization','probability')
    ylim([0 0.35])
end


nbins = 20;
edges = [0:0.4/nbins: 0.4];

for pp =1:4
    figure
    histogram(IA_mtx(:,pp),edges,'Normalization','probability')
    ylim([0 0.2])
end
