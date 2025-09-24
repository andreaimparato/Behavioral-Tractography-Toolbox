%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PURPOSE:        Preprocess EMA and SIPS data, perform dimensionality
%                 reduction (PCA) on symptom ratings, and classify subjects
%                 into two clusters using k-means on the first/second PCA dimension.
%
% INPUTS:         
%   - EMA_DEF.mat      : Ecological Momentary Assessment dataset
%   - SIPS_DEF.mat     : Structured Interview for Prodromal Syndromes dataset
%   - IQ_DEF           : Path to IQ dataset (not directly loaded in this script)
%   - sips_names_2.mat : Symptom variable names for PCA interpretation
%
% PROCESSING STEPS:
%   1. EMA preprocessing:
%        - Keep only completed sessions
%        - Retain subjects with ≥1/3 valid beeps
%        - Identify valid time points (adjacent notifications)
%   2. SIPS preprocessing:
%        - Remove ID and non-symptom variables
%        - Replace missing values
%        - Remove "P3 Grandiosity" item
%   3. Dimensionality reduction:
%        - Principal Component Analysis (PCA) on SIPS items
%   4. Clustering:
%        - k-means clustering (k=2) performed on the first PCA dimension
%        - Evaluation of cluster validity with silhouette criterion
%   5. Visualization:
%        - Scatterplot of subjects in PCA space
%        - Clusters labeled as Cluster 1 and Cluster 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Initialization
clc                % Clear command window
clear all          % Remove all variables from workspace
close all          % Close all open figure windows

%% Load datasets
% Load ecological momentary assessment (EMA) and SIPS symptom data
load C:\Users\ACER\Documents\Clinic_and_networks\EMA_DEF
load C:\Users\ACER\Documents\Clinic_and_networks\SIPS_DEF

% Define path for IQ data (not yet loaded here)
iq_file = 'C:\Users\ACER\Documents\Clinic_and_networks\IQ_DEF';

%% EMA preprocessing: completion, filtering, valid time points
% Keep only completed EMA sessions
EMA = EMA(EMA.CompletedSession == 1,:);

% Retain participants with at least 1/3 of the maximum expected beeps
max_beeps = 48;
threshold = floor(max_beeps/3);
u = unique(EMA.subj_id); % Unique subject IDs
s = arrayfun(@(id) sum(EMA.CompletedSession(EMA.subj_id == id)), u);
to_keep = ~ismember(EMA.subj_id, u(s < threshold));
EMA = EMA(to_keep,:);

% Identify valid time points (adjacent notifications)
vtp = false(numel(EMA.NotificationNo),1);
vtp(1) = abs(EMA.NotificationNo(1) - EMA.NotificationNo(2)) == 1;
for i = 2:numel(EMA.NotificationNo)-1
    prev_adj = abs(EMA.NotificationNo(i) - EMA.NotificationNo(i-1)) == 1;
    next_adj = abs(EMA.NotificationNo(i) - EMA.NotificationNo(i+1)) == 1;
    vtp(i) = prev_adj || next_adj;
end
vtp(end) = abs(EMA.NotificationNo(end) - EMA.NotificationNo(end-1)) == 1;

% Keep only valid EMA entries
EMA = EMA(vtp,:);

% Extract unique subject IDs and corresponding age values
[unique_ema_subj_id, idx] = unique(EMA.subj_id,'stable');
unique_age_y_ema = EMA.age_years(idx);

%% Symptom table (SIPS) preprocessing
ts_names = SIPS.Properties.VariableNames;

% Prepare SIPS data for PCA and clustering
input_pca = SIPS;
input_pca(:,[1,21]) = [];       % Remove ID and non-symptom variable
input_pca = table2array(input_pca);
input_pca(isnan(input_pca)) = 0; % Replace NaNs with 0

% Extract age variable
ts_age = SIPS.VCFSDatabase20_01_20052__Age;

%% PCA on symptom data
% Load symptom names for plotting/interpretation
load C:\Users\ACER\Documents\Clinic_and_networks\sips_names_2.mat
ts_names = strrep(sips_names_2,'_',' '); % Clean underscores in labels

% Frequency of non-zero values per symptom
freq_not_zero_sips = sum(input_pca>0)/size(input_pca,1);

% Remove "P3 Grandiosity" item
to_delete = strcmp(ts_names, 'P3 Grandiosity');
input_pca(:,to_delete) = [];
ts_names(to_delete) = [];

% Perform principal component analysis
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(input_pca);
XYZ_pca = SCORE(:,1:3); % First 3 PCs for visualization

%% Clustering analysis
rng(0); % For reproducibility

% Evaluate clustering quality across 1–20 clusters (silhouette method)
eva = evalclusters(input_pca,'kmeans','silhouette','KList',1:20);

% Define colormap for 2-cluster solution
colormap_1  = [0 0.4470 0.7410; 1 0 0]; % Blue, Red

% Compute PCA scores from original data
score_from_original = (input_pca - mean(input_pca))*COEFF;

% Select PCA dimension(s) for clustering (currently 1st dimension only)
ik = score_from_original(:,1);
cluster_names = {'LSI', 'HSI'}; % Labels: Low vs. High Symptom Intensity

% Evaluate silhouette across 1–20 clusters on the selected dimension
eva = evalclusters(ik,'kmeans','silhouette','KList',1:20);

% Keep first 2 PCA dimensions for visualization
XYZ_pca = score_from_original(:,1:2);

% Perform k-means clustering (k=2)
[clust, centroids] = kmeans(ik, 2, 'OnlinePhase','on', 'MaxIter', 10000);

%% Visualization of PCA and cluster results
figure
scatter_handles = gobjects(2,1); % Preallocate for legend handles
for i = 1:numel(unique(clust))
    % Extract coordinates of points belonging to cluster i
    x = XYZ_pca(clust == i,1);
    y = XYZ_pca(clust == i,2);
    n = SIPS.Sujet_ID(clust == i); % Subject IDs
    c = colormap_1(i,:);           % Cluster color
    
    % Scatter plot of subjects by cluster
    scatter_handles(i) = scatter(x, y, 35, c, 'filled', ...
        'Marker', '^', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1);
    xlabel('1 dim: Symptoms Severity')
    ylabel('2 dim: Positive - Negative')
    zlabel('3 dim')
    hold on
    
    % Plot cluster centroids
    plot(centroids(i,1), 1, 'o', 'Color',c, 'MarkerSize',15, ...
        'MarkerFaceColor',c, 'MarkerEdgeColor', 'k');
    
    % Add subject IDs as labels
    text(x, y, string(n), 'FontSize', 8, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
legend(scatter_handles, cluster_names, 'Location','northeast');
