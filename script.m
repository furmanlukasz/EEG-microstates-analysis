clearvars;
close all;

% 64 channels labels
channels = {'FC5','FC3','FC1','FCZ','FC2','FC4','FC6','C5','C3','C1','CZ',...
    'C2' ,'C4' ,'C6' ,'CP5' ,'CP3' ,'CP1','CPZ','CP2','CP4','CP6','FP1','FPZ',...
    'FP2','AF7','AF3','AFZ' ,'AF4' ,'AF8','F7','F5','F3','F1','FZ' ,'F2','F4',...
    'F6','F8','FT7','FT8','T7','T8','T9','T10','TP7','TP8','P7' ,'P5',...
    'P3' ,'P1','PZ' ,'P2' ,'P4' ,'P6' ,'P8' ,'PO7','PO3','POZ','PO4','PO8',...
    'O1' ,'OZ' ,'O2','IZ'};  

% Color map red-blu
map = zeros(21,3);
map(21,:) = [1,0,0];
map(20,:) = [1,1/10,1/10];
map(19,:) = [1,2/10,2/10];
map(18,:) = [1,3/10,3/10];
map(17,:) = [1,4/10,4/10];
map(16,:) = [1,5/10,5/10];
map(15,:) = [1,6/10,6/10];
map(14,:) = [1,7/10,7/10];
map(13,:) = [1,8/10,8/10];
map(12,:) = [1,9/10,9/10];
map(11,:) = [1,1,1];
map(10,:) = [9/10,9/10,1];
map(9,:) = [8/10,8/10,1];
map(8,:) = [7/10,7/10,1];
map(7,:) = [6/10,6/10,1];
map(6,:) = [5/10,5/10,1];
map(5,:) = [4/10,4/10,1];
map(4,:) = [3/10,3/10,1];
map(3,:) = [2/10,2/10,1];
map(2,:) = [1/10,1/10,1];
map(1,:) = [0,0,1];

col = ['r','b','g','c','y'];

n_channels = 64;
T = 9632;

% "head plot" utility function
addpath('plot_topography');

% Selected subject data
[hdr, record] = edfread('S092R02.edf');
record = record(1:n_channels,1:T);


% ------- Alternative filtering solutions --------- %
% fs = 160;
% fpass = [2 25];
% record = bandpass(record, fpass, fs);
% 
% [b, a] = butter(4, fpass/80, 'bandpass');
% y = filtfilt(b,a,record);
% record = y;
% ---------------------------------------------------------- %

record_t = record'; 


%% Finding instants of max variation from rest condition: GFP

% Compute Global Field Potential (GFP)
V_mean = mean(record_t,2);   % instantaneous potential averaged between electrodes
GFP = zeros(1,T);
for i = 1 : T
    GFP(i) = sqrt((sum((record_t(i,:)-V_mean(i)).^2))/n_channels);
end

% Compute GFP peaks
[pks, locs] = findpeaks(GFP);

% Segnali eeg dei diversi canali nei picchi di GFP
eeg_at_peaks = record_t(locs,:);

%% Example: eeg signal from certain channels and GFP

% to_plot = eeg signals together with GFP
to_plot = record_t;               
to_plot(:,65) = GFP';
figure; 
% 80 samples plot: 0.5s time interval (fs=160/s)
s = stackedplot(to_plot(100:180,[1,12,27,34,45,51,63,65]));         
s.DisplayLabels = {'FC5','C2','AFZ','FZ','TP7','PZ','O2','GFP'};
s.LineProperties(8).Color = 'red';

%% Principal Component Analysis
[U,X,S]=pca(zscore(eeg_at_peaks));

% Pareto
figure
bar(100*S./sum(S))
xlabel('# components')
ylabel('variance explained [%]')
box off

% Lorenz
figure
plot(100*cumsum(S)./sum(S))
xlabel('# components')
ylabel('variance explained [%]')
box off

figure
plot3(X(:,1),X(:,2),X(:,3),'.b')
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
axis equal
hold on
% Biplot
biplot(100*U(:,1:2),'varlabels',cellstr(channels))

%% Cluster Analysis 

% Optimal number of cluster evaluation 
E1 = evalclusters(eeg_at_peaks, 'kmeans', 'silhouette','klist', 1:7);
figure; plot(E1); 
% k = E1.OptimalK; 

% [idx, C] = kmeans(eeg_at_peaks, k);
% figure; silhouette(eeg_at_peaks, idx);

k = 4; 
% By looking at data in PC space we can assume K=4
idx = kmeans(eeg_at_peaks,k);  
cols = {'b','r','g','m','c'};
figure
for ii = 1:k
    is_ii = find(idx == ii);
    mu_ii = nanmean(X(is_ii,1:3));
    hold on
    plot3(X(is_ii,1),X(is_ii,2),X(is_ii,3),'.','MarkerFaceColor',cols{ii},'MarkerEdgeColor',cols{ii})
    plot3(mu_ii(1),mu_ii(2),mu_ii(3),'bo','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    view(3)
    xlabel('PC_1')
    ylabel('PC_2')
    zlabel('PC_3')
end

k = 8;
for ii = 1:k
    [idx,C,sumd] = kmeans(eeg_at_peaks,k);
    hold on
    figure
    scatter(ii, sum(sumd),'red','filled')
    xlabel('# clusters')
    ylabel('avg distance')
end

k = 8;
figure
for ii = 1:k
    [idx,C,sumd] = kmeans(eeg_at_peaks,k);
    hold on
    scatter(ii,sum(sumd),'red','filled')
    xlabel('# clusters')
    ylabel('avg distance')
end

% --------------------- Gaussian Mixture ------------------------------
k = 4;
GMModel = fitgmdist(eeg_at_peaks,k);
P = posterior(GMModel,eeg_at_peaks);
idxi = nan(1,size(eeg_at_peaks,1));
for ii = 1:size(eeg_at_peaks,1)
    [m,idx] = max(P(ii,:));
    idxi(ii) = idx;
end

figure
for ii = 1:k
    is_ii = find(idxi == ii);
    mu_ii = nanmean(X(is_ii,1:3));
    hold on
    plot3(X(is_ii,1),X(is_ii,2),X(is_ii,3),'o','MarkerFaceColor',col(ii),'MarkerEdgeColor',col(ii))
    xlabel('PC_1')
    ylabel('PC_2')
    zlabel('PC_3')
    view(3)
end

k = 5;
figure
for ii = 1:k
    GMModel = fitgmdist(eeg_at_peaks,ii);
    [P,nlogL] = posterior(GMModel,eeg_at_peaks);
    hold on
    scatter(ii,nlogL,'red','filled')
    xlabel('# clusters')
    ylabel('-log(likelihood)')
end

% -------------------------------------------------------------- %

k = 4;
[idx,C] = kmeans(eeg_at_peaks,k);

% Identified topographic maps (associated to a microstate) plot 
for i = 1 : k
   figure; plot_topography(channels, C(i,:));
   colormap(map);
end

%% Association of a general "template" map to each instant of the original eeg signals
% Each instant of the original eeg signal is uniquely associated
% with a microstate among those identified.

ms_at_peaks = zeros(length(eeg_at_peaks),1);
for i = 1 : length(eeg_at_peaks)
    corrcoefs = zeros(1,k);
    for j = 1 : k
        M = corrcoef(eeg_at_peaks(i,:),C(j,:));
        % Compute correlation between instant 'i' of eeg signal and map 'j'
        corrcoefs(1,j) = M(1,2);
    end
    index = find(corrcoefs==max(corrcoefs));
    ms_at_peaks(i,1) = index;
    % Microstates index with max correlation with eeg signal at instant 'i' 
    % is assigned to the same instant
end

% We assume that the point of change between two consecutive microstates/maps 
% is at half way between the two peaks

% il punto di cambio di microstato/mappa è a metà tra due picchi i cui 
% microstati/mappe corrispondenti sono diversi.
ms_sequence = init_ms_seq(ms_at_peaks,locs);

%% Statistical Analysis

% Average Lifespan of MS
ms_avg_lifespan = zeros(1,k);

% Frequency of Appearance of MS
f_appear = zeros(1,k);

% Fraction Total Covered Time of MS
total_covered_t_ms = zeros(1,k);

% Compute parameters for each microstates 
seq = cell(1,k); % vettori contenenti le lunghezze di ogni sequenza per ogni microstato
for i = 1 : k
   idx = 1;
   l = 0;
   % Lunghezze delle sequenze di microstati   
   for j = 1 : length(ms_sequence)
        if ms_sequence(j,1) == i
            l = l+1;
        elseif ms_sequence(j,1) ~= i && l~=0
             seq_lengths(idx) = l; 
             l = 0;
             idx = idx+1;
        end
   end
   
   seq(1,i) = {seq_lengths};
   
   ms_avg_lifespan(1,i) = (mean(seq_lengths)/160*10^3); 
                                                        
   n = find(ms_sequence==i);
   n_peaks = find(ms_at_peaks==i);
   f_appear(1,i) = length(n_peaks)/(9632/160);    
   
   total_covered_t_ms(1,i) = length(n)/9632*100; 
   
end

T = table(ms_avg_lifespan',f_appear',total_covered_t_ms', 'RowNames',...
    {'1';'2';'3';'4'}, 'VariableNames', {'Avg_Lifespan'; 'Freq_Appear'; 'Tot_Covered_time'});
figure;
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1], 'ColumnWidth', {60 60 60});
set(gcf,'position',[500,500,230,130]);

%% Plot GFP con sequenze di microstati 

count = zeros(1,k);
figure;
width=900;
height=150;
for i = 1000 : 1160%length(ms_sequence)/80
    k = ms_sequence(i,1);
    count(1,k) = count(1,k)+1;
    seq_l = cell2mat(seq(1,k));
    if count(1,k) < 601
        ar = area(i:i+seq_l(1,count(1,k)),GFP(1,i:i+seq_l(1,count(1,k))));
        ar.FaceColor = col(k);
        ar.EdgeColor = col(k);
        hold on;
        i = i+seq_l(1,count(1,k));
    end
    xlabel('t');
    ylabel('GFP');
    set(gcf,'position',[100,100,width,height])
end
