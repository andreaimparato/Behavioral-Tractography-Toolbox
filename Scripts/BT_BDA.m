clc
clear all
close all


%% load and set variables
for iter = 1:1
    
% load matrix (MLN matrix computed in MLN_network_1)
load matrix_1.mat
cross_net(cross_net<0) = 0;
NET20_CONTROL = cross_net;
longi_net(longi_net<0) = 0;
NET40_CONTROL = longi_net;

% load names
load names_3_v4

% load coordinates (PCA first and second dimesion as computed in
% MLN_Network_1 and stored in "XY_CROSS" variable
load coord_1
XY_CROSS_x2 = [XY_CROSS;XY_CROSS];

% set variables
dim_cross = max(size(NET40_CONTROL))/2;
dim_longi = max(size(NET40_CONTROL));

SOURCE = 1:dim_cross;
TARGET = (dim_cross + 1) : dim_longi;

% color quadrants
color_quadrants = "#000000";


end


%% choose population
for iter = 1:1

net1 = NET40_CONTROL;
net1_cross = NET20_CONTROL;

end


%% creation of the final adj matrix
for iter = 1:1

net1(1:dim_cross, 1:dim_cross) = net1_cross;
net1(dim_cross + 1:dim_longi, dim_cross + 1:dim_longi) = net1_cross;

end


%% delete connections from layer2 to layer1
for iter = 1:1

net1(TARGET,SOURCE)=0;

end


%% paths
for iter = 1:1
% each path is characterized by 4 coordinates: start, exit node layer1, entrance node layer2, end.

% paths net1
D=1./posweights(net1);
[SPL,hops,Pmat] = distance_wei_floyd(D);
c = 1;
B_net1 = zeros(1,10);
for iter1=1:max(size(SOURCE))
    for iter2=1:max(size(TARGET))
        B=retrieve_shortest_path(SOURCE(iter1),TARGET(iter2),hops,Pmat);
        B_net1(c,:) = [B',zeros(1,10-max(size(B)))];
        for iter = 1:max(size(B))
            if B(iter) > dim_cross
                break
            end
        end
        coord_entry_t1_net1(c,:) = XY_CROSS_x2(B(1),:);
        coord_exit_t1_net1(c,:) = XY_CROSS_x2(B(iter-1),:);
        coord_entry_t2_net1(c,:) = XY_CROSS_x2(B(iter),:);
        coord_exit_t2_net1(c,:) = XY_CROSS_x2(B(max(size(B))),:);
        c = c + 1;
    end
end


end


%% kmeans option A
% no manual inputs
for iter = 1:1


% coord_k = [coord_entry_t1_net1, coord_exit_t1_net1, coord_entry_t2_net1, coord_exit_t2_net1];
% 
% % seed
% rng(1);
% 
% % zscore input variables
% diff_coord2 = zscore(coord_k);
% 
% % kmeans
% eva = evalclusters(coord_k,'kmeans','Gap','KList',1:10);
% 
% % number of cluster
% n_cluster = eva.OptimalK;
% 
% % cluster
% [clust, centroids] = kmeans(diff_coord2, n_cluster, 'OnlinePhase','on', 'MaxIter', 10000);
% 
% % count of sp per cluster
% for iter = 1:max(clust)
%     n_of_sp_per_clust(iter) = max(size(find(clust==iter)));
% end


end


%% kmeans option B
%number of cluster as manual input
for iter = 1:1


% coord_k = [coord_entry_t1_net1, coord_exit_t1_net1, coord_entry_t2_net1, coord_exit_t2_net1];
% 
% % seed
% rng(1);
% 
% % zscore input variables
% diff_coord2 = zscore(coord_k);
% 
% % kmeans
% %eva = evalclusters(coord_k,'kmeans','Gap','KList',1:10);
% 
% % number of cluster
% n_cluster = 9;
% 
% % cluster
% [clust, centroids] = kmeans(diff_coord2, n_cluster, 'OnlinePhase','on', 'MaxIter', 10000);
% 
% % count of sp per cluster
% for iter = 1:max(clust)
%     n_of_sp_per_clust(iter) = max(size(find(clust==iter)));
% end


end


%% kmeans option C
%cluster as manual input
for iter = 1:1

% load clusterization (if it already exists)
load cluster_hc.mat

% count of sp per cluster
for iter = 1:max(clust)
    n_of_sp_per_clust(iter) = max(size(find(clust==iter)));
end


end


%% optional: cluster reorder / matching
for iter = 1:1


for iter = 1:max(size(clust))
    if clust(iter)==1
        clust2(iter) = 2;
    end
    if clust(iter)==2
        clust2(iter) = 5;
    end
    if clust(iter)==3
        clust2(iter) = 9;
    end
    if clust(iter)==4
        clust2(iter) = 7;
    end
    if clust(iter)==5
        clust2(iter) = 4;
    end
    if clust(iter)==6
        clust2(iter) = 6;
    end
    if clust(iter)==7
        clust2(iter) = 3;
    end
    if clust(iter)==8
        clust2(iter) = 8;
    end
    if clust(iter)==9
        clust2(iter) = 1;
    end
end
clear clust
clust = clust2;

% count of sp per cluster
for iter = 1:max(clust)
    n_of_sp_per_clust(iter) = max(size(find(clust==iter)));
end


end


%% BDA
net1_bda = [];
for cluster_to_test = 1:max(clust)

    %% net 1
    %% vector
    for iter = 1:1
    
    % Compute shortest paths and distances for net1
    D = 1 ./ posweights(net1);
    [SPL, hops, Pmat] = distance_wei_floyd(D);
    
    % Initialize variables for layer 1
    % Iterate over nodes
    for HUB = 1:dim_cross
        cnd = 1;
        cc = 1;
        clear coord_entr_t1_nd
        clear coord_exit_t1_nd
        clear coord_entr_t2_nd
        clear coord_exit_t2_nd
        
        % Iterate over source and target nodes
        for iter1 = 1:max(size(SOURCE))
            for iter2 = 1:max(size(TARGET))
                % Check if the paths belongs to the cluster of interest
                if clust(cc) == cluster_to_test
                    % Retrieve the shortest path between source and target
                    B = retrieve_shortest_path(SOURCE(iter1), TARGET(iter2), hops, Pmat);
                    % Check if the node is part of the shortest path
                    if min(size(intersect(B, HUB))) > 0
                        % Extract relevant coordinates for entry and exit nodes
                        p = find(B == HUB);
                        pl = 0;
                        for iter = 1:max(size(B))
                            if B(iter) > dim_cross
                                pl = iter;
                                break
                            end
                        end
                        coord_entr_t1_nd(cnd,:) = XY_CROSS_x2(B(1),:);
                        coord_exit_t1_nd(cnd,:) = XY_CROSS_x2(B(pl-1),:);
                        coord_entr_t2_nd(cnd,:) = XY_CROSS_x2(B(pl),:);
                        coord_exit_t2_nd(cnd,:) = XY_CROSS_x2(B(max(size(B))),:);
                        cnd = cnd + 1;
                    end
                end
                cc = cc + 1;
            end
        end
        
        % Compute mean coordinates for each symptom
        if exist('coord_entr_t1_nd','var') == 1
            mL1_coord_entr_t1_nd(HUB,:) = mean(coord_entr_t1_nd(:,:),1);
            mL1_coord_exit_t1_nd(HUB,:) = mean(coord_exit_t1_nd(:,:),1);
            mL1_coord_entr_t2_nd(HUB,:) = mean(coord_entr_t2_nd(:,:),1);
            mL1_coord_exit_t2_nd(HUB,:) = mean(coord_exit_t2_nd(:,:),1);
        else
            mL1_coord_entr_t1_nd(HUB,:) = [nan, nan];
            mL1_coord_exit_t1_nd(HUB,:) = [nan, nan];
            mL1_coord_entr_t2_nd(HUB,:) = [nan, nan];
            mL1_coord_exit_t2_nd(HUB,:) = [nan, nan];
        end
    end
    
    % Similar process for layer 2
    for HUB = dim_cross+1:dim_longi
        cnd = 1;
        cc = 1;
        clear coord_entr_t1_nd
        clear coord_exit_t1_nd
        clear coord_entr_t2_nd
        clear coord_exit_t2_nd
    
        for iter1 = 1:max(size(SOURCE))
            for iter2 = 1:max(size(TARGET))
                if clust(cc) == cluster_to_test
                    B = retrieve_shortest_path(SOURCE(iter1), TARGET(iter2), hops, Pmat);
                    if min(size(intersect(B, HUB))) > 0
                        p = find(B == HUB);
                        pl = 0;
                        for iter = 1:max(size(B))
                            if B(iter) > dim_cross
                                pl = iter;
                                break
                            end
                        end
                        coord_entr_t1_nd(cnd,:) = XY_CROSS_x2(B(1),:);
                        coord_exit_t1_nd(cnd,:) = XY_CROSS_x2(B(pl-1),:);
                        coord_entr_t2_nd(cnd,:) = XY_CROSS_x2(B(pl),:);
                        coord_exit_t2_nd(cnd,:) = XY_CROSS_x2(B(max(size(B))),:);
                        cnd = cnd + 1;
                    end
                end
                cc = cc + 1;
            end
        end
        
        if exist('coord_entr_t1_nd','var') == 1
            mL2_coord_entr_t1_nd(HUB-dim_cross,:) = mean(coord_entr_t1_nd(:,:),1);
            mL2_coord_exit_t1_nd(HUB-dim_cross,:) = mean(coord_exit_t1_nd(:,:),1);
            mL2_coord_entr_t2_nd(HUB-dim_cross,:) = mean(coord_entr_t2_nd(:,:),1);
            mL2_coord_exit_t2_nd(HUB-dim_cross,:) = mean(coord_exit_t2_nd(:,:),1);
        else
            mL2_coord_entr_t1_nd(HUB-dim_cross,:) = [nan, nan];
            mL2_coord_exit_t1_nd(HUB-dim_cross,:) = [nan, nan];
            mL2_coord_entr_t2_nd(HUB-dim_cross,:) = [nan, nan];
            mL2_coord_exit_t2_nd(HUB-dim_cross,:) = [nan, nan];
        end
    end
    
    end
    
    %% Save variables
    for iter = 1:1
    
    % Save computed values for net1
    net1_mL1_coord_entr_t1_nd = mL1_coord_entr_t1_nd;
    net1_mL1_coord_exit_t1_nd = mL1_coord_exit_t1_nd;
    net1_mL1_coord_entr_t2_nd = mL1_coord_entr_t2_nd;
    net1_mL1_coord_exit_t2_nd = mL1_coord_exit_t2_nd;
    
    net1_mL2_coord_entr_t1_nd = mL2_coord_entr_t1_nd;
    net1_mL2_coord_exit_t1_nd = mL2_coord_exit_t1_nd;
    net1_mL2_coord_entr_t2_nd = mL2_coord_entr_t2_nd;
    net1_mL2_coord_exit_t2_nd = mL2_coord_exit_t2_nd;

    a = [net1_mL1_coord_entr_t1_nd; net1_mL1_coord_exit_t1_nd; net1_mL1_coord_entr_t2_nd; net1_mL1_coord_exit_t2_nd; 
        net1_mL2_coord_entr_t1_nd; net1_mL2_coord_exit_t1_nd; net1_mL2_coord_entr_t2_nd; net1_mL2_coord_exit_t2_nd];

    % net1_bda is the variables where BDA metrics are stored.
    % In each pair of columns, there is the BDA (Bidimensional Data Array) related to a cluster.
    % For example, the first two columns represent cluster 1.
    % Within each pair, the first column represents the x coordinate while the second column represents the y coordinate.
    % Each row represents a symptom, and every twenty rows, the coordinate in question changes
    % (e.g., net1_mL1_coord_entr_t1_nd, net1_mL1_coord_exit_t1_nd, ...).
    % For instance, element (22,3) represents the average exit x-coordinate from layer 1, related to the second symptom (loneliness).

    net1_bda = [net1_bda, a];
    
    end

end


%% Define colors 
for iter = 1:1

% Define colors for different networks
color_t1 = "#008000"; % Color for T1 network
color_longi = "#7CFC00"; % Color for longitudinal network
color_t2 = "#FDDA0D"; % Color for T2 network
color_allnet = "#A9A9A9"; % Color for all networks

end


%% choose cluster to plot
for iter = 1:1

    iter_cluster = 5;
    hub_to_plot = [1:1:dim_longi]; % Nodes to plot
    cl_to_plot = iter_cluster; % Cluster to plot
    plot_avg_vector = 1; % Flag to plot average vector

end


%% coordinates avg vector
for iter=1:1

    % Compute average entry and exit coordinates for a specific cluster
    coord_avgv_entry_t1 = mean(coord_entry_t1_net1(find(clust==cl_to_plot),:));
    coord_avgv_exit_t1 = mean(coord_exit_t1_net1(find(clust==cl_to_plot),:));
    coord_avgv_entry_t2 = mean(coord_entry_t2_net1(find(clust==cl_to_plot),:));
    coord_avgv_exit_t2 = mean(coord_exit_t2_net1(find(clust==cl_to_plot),:));

end


%% Create structure to plot CLUSTER vectors and shortpaths (s1)
% This section creates the structure that will then be exported as a json. 
for iter = 1:1
    


%% input data for the export
for iter=1:1


%% adj matrix    
% Zero Diagonal: Set diagonal elements of net1 to zero
for iter = 1:max(size(net1))
    net1(iter,iter) = 0;
end

% Remove Connections: Remove connections from TARGET nodes to SOURCE nodes
net1(TARGET,SOURCE) = 0;

% Compute Shortest Paths:
% Invert the weights of net1 to get distance matrix D
% Compute shortest path lengths SPL, number of hops hops, and shortest path matrices Pmat using Floyd's algorithm
D = 1./posweights(net1);
[SPL, hops, Pmat] = distance_wei_floyd(D);

% Initialize Variables:
% Initialize matrices ACROSS_BETW, BTW, and nc to store network properties
ACROSS_BETW = zeros(max(size(net1)), max(size(net1)));
BTW = zeros(max(size(net1)), 1);
nc = zeros(max(size(net1)), 1);

cc = 1; % Counter for cluster index
aa = 1; % Counter for hub index
for iter1 = 1:max(size(SOURCE))
    for iter2 = 1:max(size(TARGET))
        B = retrieve_shortest_path(SOURCE(iter1), TARGET(iter2), hops, Pmat);
        % Check if the current cluster matches cl_to_plot and there's an intersection between the path B and the hub nodes hub_to_plot
        if clust(cc) == cl_to_plot
            if min(size(intersect(B, hub_to_plot))) > 0
                % Adjacency Matrix: Construct a binary adjacency matrix ACROSS_BETW2 for the current path B
                ACROSS_BETW2 = zeros(max(size(net1)), max(size(net1)));
                for iter3 = 1:(max(size(B)) - 1)
                    ACROSS_BETW2(B(iter3), B(iter3 + 1)) = 1;
                end
                % Update ACROSS_BETW by adding ACROSS_BETW2
                ACROSS_BETW = ACROSS_BETW + ACROSS_BETW2;

                % Betweenness Centrality: Update BTW by setting nodes in B to 1
                BTW2 = zeros(max(size(net1)), 1);
                BTW2(B) = 1;
                BTW = BTW + BTW2;

                % Nodes Color: Update node colors nc based on node membership
                for iter4 = 1:(max(size(B)))
                    if nc(B(iter4)) ~= 1 && nc(B(iter4)) ~= 4
                        if B(iter4) <= dim_cross
                            nc(B(iter4)) = 2;
                        else
                            nc(B(iter4)) = 3;
                        end
                    end
                end
                % Set the first and last nodes of B to colors 1 and 4 respectively
                nc(B(1)) = 1;
                nc(B(max(size(B)))) = 4;
            end
        end
        cc = cc + 1;
    end
end



%% node
% Calculate node sizes based on the logarithm of betweenness centrality
for iter = 1:dim_cross
    node_size1(iter) = log(BTW(iter));
    node_size2(iter) = log(BTW(iter+dim_cross));
end

% Apply threshold to node sizes
for iter = 1:dim_cross
    if node_size1(iter) <= 1
        node_size1(iter) = 1;
    end
    if node_size2(iter) <= 1
        node_size2(iter) = 1;
    end
end

% Create unique node names for layers 1 and 2
count = 1;
for iter=1:20
    nodes_unic_names_layer_1(count) = strcat(names_3(iter),'+Layer1');
    nodes_unic_names_layer_2(count) = strcat(names_3(iter),'+Layer2');
    count = count + 1;
end
nodes_unic_names_layer_1 = string(nodes_unic_names_layer_1)';
nodes_unic_names_layer_2 = string(nodes_unic_names_layer_2)';

% Assign node colors
for iter = 1:max(size(nc))
    if nc(iter) == 1
        nodecolors(iter,:) = color_t1;
    elseif nc(iter) == 2 || nc(iter) == 3
        nodecolors(iter,:) = color_longi;
    elseif nc(iter) == 4
        nodecolors(iter,:) = color_t2;
    else
        nodecolors(iter,:) = "#E0E0E0";
    end
end

% Adjust node coordinates to center them around (0,0) and scale them
max_x = max(XY_CROSS(:,1));
max_y = max(XY_CROSS(:,2));
min_x = min(XY_CROSS(:,1));
min_y = min(XY_CROSS(:,2));
centre_x = (max_x + min_x)/2;
centre_y = (max_y + min_y)/2;

XY_CROSS(:,1) = XY_CROSS(:,1) - centre_x;
XY_CROSS(:,2) = XY_CROSS(:,2) - centre_y;
XY_CROSS(:,1) = XY_CROSS(:,1)*1.2;
XY_CROSS(:,2) = XY_CROSS(:,2)*1.2;

% Adjust average entry and exit coordinates similarly
coord_avgv_entry_t1(:,1) = coord_avgv_entry_t1(:,1) - centre_x;
coord_avgv_entry_t1(:,2) = coord_avgv_entry_t1(:,2) - centre_y;
coord_avgv_exit_t1(:,1) = coord_avgv_exit_t1(:,1) - centre_x;
coord_avgv_exit_t1(:,2) = coord_avgv_exit_t1(:,2) - centre_y;
coord_avgv_entry_t2(:,1) = coord_avgv_entry_t2(:,1) - centre_x;
coord_avgv_entry_t2(:,2) = coord_avgv_entry_t2(:,2) - centre_y;
coord_avgv_exit_t2(:,1) = coord_avgv_exit_t2(:,1) - centre_x;
coord_avgv_exit_t2(:,2) = coord_avgv_exit_t2(:,2) - centre_y;

coord_avgv_entry_t1(:,1) = coord_avgv_entry_t1(:,1)*1.2;
coord_avgv_entry_t1(:,2) = coord_avgv_entry_t1(:,2)*1.2;
coord_avgv_exit_t1(:,1) = coord_avgv_exit_t1(:,1)*1.2;
coord_avgv_exit_t1(:,2) = coord_avgv_exit_t1(:,2)*1.2;
coord_avgv_entry_t2(:,1) = coord_avgv_entry_t2(:,1)*1.2;
coord_avgv_entry_t2(:,2) = coord_avgv_entry_t2(:,2)*1.2;
coord_avgv_exit_t2(:,1) = coord_avgv_exit_t2(:,1)*1.2;
coord_avgv_exit_t2(:,2) = coord_avgv_exit_t2(:,2)*1.2;

% Determine label positions (top or bottom) based on node coordinates
for iter = 1:max(size(XY_CROSS(:,2)))
    if XY_CROSS(iter,2) < 0
        labelPosition(iter) = "bottom";
    else
        labelPosition(iter) = "top";
    end
end

% Calculate label lengths and adjust X positions accordingly
f = 130/11; % Scaling factor

for iter = 1:20
    label_lenght(iter) = strlength(names_3(iter));
end

for iter = 1:max(size(XY_CROSS(:,1)))
    if XY_CROSS(iter,1) < 0
        labelAdjustX(iter) = -label_lenght(iter)*f;
    else
        labelAdjustX(iter) = "0";
    end
end

% Adjust specific label positions manually
labelAdjustX(11) = "0"; % feeling unsafe
labelAdjustX(12) = "0"; % confusion
labelAdjustX(4) = -(label_lenght(4)*f+100);
labelAdjustX = [labelAdjustX,labelAdjustX];
labelAdjustX(24) = 0;



%% edges

% Define NAMES_NEW2 as names_3
NAMES_NEW2 = names_3;

% Generate source and target names for edges within Layer 1
count = 1;
for iter = 1:20
    for iter2 = 1:20
        src_names_edges_cross_1(count) = strcat(NAMES_NEW2(iter), '+Layer1');
        trg_names_edges_cross_1(count) = strcat(NAMES_NEW2(iter2), '+Layer1');
        count = count + 1;
    end
end

% Convert source and target names to strings
src_names_edges_cross_1 = src_names_edges_cross_1';
src_names_edges_cross_1 = string(src_names_edges_cross_1);
trg_names_edges_cross_1 = trg_names_edges_cross_1';
trg_names_edges_cross_1 = string(trg_names_edges_cross_1);

% Generate source and target names for edges from Layer 1 to Layer 2
count = 1;
for iter = 1:20
    for iter2 = 21:40
        src_names_edges_longi_1_2(count) = strcat(NAMES_NEW2(iter), '+Layer1');
        trg_names_edges_longi_1_2(count) = strcat(NAMES_NEW2(iter2 - 20), '+Layer2');
        count = count + 1;
    end
end

% Convert source and target names to strings
src_names_edges_longi_1_2 = src_names_edges_longi_1_2';
src_names_edges_longi_1_2 = string(src_names_edges_longi_1_2);
trg_names_edges_longi_1_2 = trg_names_edges_longi_1_2';
trg_names_edges_longi_1_2 = string(trg_names_edges_longi_1_2);

% Generate source and target names for edges within Layer 2
count = 1;
for iter = 1:20
    for iter2 = 1:20
        src_names_edges_cross_2(count) = strcat(NAMES_NEW2(iter), '+Layer2');
        trg_names_edges_cross_2(count) = strcat(NAMES_NEW2(iter2), '+Layer2');
        count = count + 1;
    end
end

% Convert source and target names to strings
src_names_edges_cross_2 = src_names_edges_cross_2';
src_names_edges_cross_2 = string(src_names_edges_cross_2);
trg_names_edges_cross_2 = trg_names_edges_cross_2';
trg_names_edges_cross_2 = string(trg_names_edges_cross_2);

% Calculate edge weights
W = ACROSS_BETW;

% Generate edge weights and colors for edges within Layer 1
count = 1;
for iter = 1:20
    for iter2 = 1:20
        WEIGHT_ARENA1(count) = W(iter, iter2);
        COLORS1(count) = color_t1;
        count = count + 1;
    end
end
WEIGHT_ARENA1 = WEIGHT_ARENA1';

% Generate edge weights and colors for edges from Layer 1 to Layer 2
count = 1;
for iter = 1:20
    for iter2 = 21:40
        WEIGHT_ARENA1_2(count) = W(iter, iter2);
        COLORS1_2(count) = color_longi;
        count = count + 1;
    end
end
WEIGHT_ARENA1_2 = WEIGHT_ARENA1_2';

% Generate edge weights and colors for edges within Layer 2
count = 1;
for iter = 21:40
    for iter2 = 21:40
        WEIGHT_ARENA2(count) = W(iter, iter2);
        COLORS2(count) = color_t2;
        count = count + 1;
    end
end
WEIGHT_ARENA2 = WEIGHT_ARENA2';

% Clear unnecessary variables
clear a b c d;

% Remove edges with zero or negative weights
src_names_edges_cross_1(WEIGHT_ARENA1 == 0) = [];
trg_names_edges_cross_1(WEIGHT_ARENA1 == 0) = [];
COLORS1(WEIGHT_ARENA1 == 0) = [];
WEIGHT_ARENA1(WEIGHT_ARENA1 == 0) = [];

src_names_edges_cross_1(WEIGHT_ARENA1 < 0) = [];
trg_names_edges_cross_1(WEIGHT_ARENA1 < 0) = [];
COLORS1(WEIGHT_ARENA1 < 0) = [];
WEIGHT_ARENA1(WEIGHT_ARENA1 < 0) = [];

src_names_edges_longi_1_2(WEIGHT_ARENA1_2 == 0) = [];
trg_names_edges_longi_1_2(WEIGHT_ARENA1_2 == 0) = [];
COLORS1_2(WEIGHT_ARENA1_2 == 0) = [];
WEIGHT_ARENA1_2(WEIGHT_ARENA1_2 == 0) = [];

src_names_edges_longi_1_2(WEIGHT_ARENA1_2 < 0) = [];
trg_names_edges_longi_1_2(WEIGHT_ARENA1_2 < 0) = [];
COLORS1_2(WEIGHT_ARENA1_2 < 0) = [];
WEIGHT_ARENA1_2(WEIGHT_ARENA1_2 < 0) = [];

src_names_edges_cross_2(WEIGHT_ARENA2 == 0) = [];
trg_names_edges_cross_2(WEIGHT_ARENA2 == 0) = [];
COLORS2(WEIGHT_ARENA2 == 0) = [];
WEIGHT_ARENA2(WEIGHT_ARENA2 == 0) = [];

src_names_edges_cross_2(WEIGHT_ARENA2 < 0) = [];
trg_names_edges_cross_2(WEIGHT_ARENA2 < 0) = [];
COLORS2(WEIGHT_ARENA2 < 0) = [];
WEIGHT_ARENA2(WEIGHT_ARENA2 < 0) = [];



end




%% field 3 layers
for iter = 1:1

field3 = 'layers';

c1.name = 'Layer1';
c1.position_x = '-300';
c1.position_y = '0';
c1.position_z = '0';
c1.last_layer_scale = '1';
c1.rotation_x = '0';
c1.rotation_y = '0';
c1.rotation_z = '0';
c1.floor_current_color = '#777777';
c1.geometry_parameters_width = '1001.90476190476';


c2.name = 'Layer2';
c2.position_x = '300';
c2.position_y = '0';
c2.position_z = '0';
c2.last_layer_scale = '1';
c2.rotation_x = '0';
c2.rotation_y = '0';
c2.rotation_z = '0';
c2.floor_current_color = '#777777';
c2.geometry_parameters_width = '1001.90476190476';



value3 = struct([c1,c2]);

end


%% field 4 nodes
for iter = 1:1

field4 = 'nodes';


% nodes 1 %%%%%%%%%
nodes1 = [];
for iter = 1:20
    nodes1(iter).name = string(names_3(iter));
    nodes1(iter).labelPosition = "Top";
    nodes1(iter).labelAdjustX = '0';
    nodes1(iter).labelAdjustY = '0';
    nodes1(iter).isBold = "true";
    %nodes1(iter).isBold = bold_no_corr_str(iter);
    nodes1(iter).layer = 'Layer1';
    nodes1(iter).position_x = '0';
    nodes1(iter).position_y = string(XY_CROSS(iter,2)*1000);
    nodes1(iter).position_z = string(XY_CROSS(iter,1)*1000);
    %nodes1(iter).scale_x = string(BTW_size(iter)/25);
    nodes1(iter).scale_x = node_size1(iter);
    nodes1(iter).color = nodecolors(iter);
    nodes1(iter).url = '';
    nodes1(iter).descr = '';
    nodes1(iter).unic_name = nodes_unic_names_layer_1(iter);
    nodes1(iter).labelSize = '8';
end


% nodes 2 %%%%%%%%%
nodes2 = [];
for iter = 1:20
    nodes2(iter).name = string(names_3(iter));
    nodes2(iter).labelPosition = "Top";
    nodes2(iter).labelAdjustX = '0';
    nodes2(iter).labelAdjustY = '0';
    nodes2(iter).isBold = "true";
    %nodes2(iter).isBold = bold_no_corr_str(iter+20);
    nodes2(iter).layer = 'Layer2';
    nodes2(iter).position_x = '0';
    nodes2(iter).position_y = string(XY_CROSS(iter,2)*1000);
    nodes2(iter).position_z = string(XY_CROSS(iter,1)*1000);
    %nodes2(iter).scale_x = string(BTW_size(iter+20)/25);
    nodes2(iter).scale_x = node_size2(iter);
    nodes2(iter).color = nodecolors(iter+dim_cross);
    nodes2(iter).url = '';
    nodes2(iter).descr = '';
    nodes2(iter).unic_name = nodes_unic_names_layer_2(iter);
    nodes2(iter).labelSize = '8';
end




%intelaiatura L1
for iter = 1:1


m = 0;
m = max(size(nodes_unic_names_layer_1));

%vertex
for iter = 1:1

nodes1(m+1).name = "4";
nodes1(m+1).labelPosition = "Top";
nodes1(m+1).labelAdjustX = '0';
nodes1(m+1).labelAdjustY = '0';
nodes1(m+1).isBold = "true";
nodes1(m+1).layer = 'Layer1';
nodes1(m+1).position_x = '0';
nodes1(m+1).position_y = "500";
nodes1(m+1).position_z = "500";
nodes1(m+1).scale_x = "0";
nodes1(m+1).color = "#000000";
nodes1(m+1).url = '';
nodes1(m+1).descr = '';
nodes1(m+1).unic_name = "4+Layer1";
nodes1(m+1).labelSize = '';

nodes1(m+2).name = "3";
nodes1(m+2).labelPosition = "Top";
nodes1(m+2).labelAdjustX = '0';
nodes1(m+2).labelAdjustY = '0';
nodes1(m+2).isBold = "true";
nodes1(m+2).layer = 'Layer1';
nodes1(m+2).position_x = '0';
nodes1(m+2).position_y = "500";
nodes1(m+2).position_z = "-500";
nodes1(m+2).scale_x = "0";
nodes1(m+2).color = "#000000";
nodes1(m+2).url = '';
nodes1(m+2).descr = '';
nodes1(m+2).unic_name = "3+Layer1";
nodes1(m+2).labelSize = '0';

nodes1(m+3).name = "1";
nodes1(m+3).labelPosition = "Top";
nodes1(m+3).labelAdjustX = '0';
nodes1(m+3).labelAdjustY = '0';
nodes1(m+3).isBold = "true";
nodes1(m+3).layer = 'Layer1';
nodes1(m+3).position_x = '0';
nodes1(m+3).position_y = "-500";
nodes1(m+3).position_z = "-500";
nodes1(m+3).scale_x = "0";
nodes1(m+3).color = "#000000";
nodes1(m+3).url = '';
nodes1(m+3).descr = '';
nodes1(m+3).unic_name = "1+Layer1";
nodes1(m+3).labelSize = '0';

nodes1(m+4).name = "2";
nodes1(m+4).labelPosition = "Top";
nodes1(m+4).labelAdjustX = '0';
nodes1(m+4).labelAdjustY = '0';
nodes1(m+4).isBold = "true";
nodes1(m+4).layer = 'Layer1';
nodes1(m+4).position_x = '0';
nodes1(m+4).position_y = "-500";
nodes1(m+4).position_z = "500";
nodes1(m+4).scale_x = "0";
nodes1(m+4).color = "#000000";
nodes1(m+4).url = '';
nodes1(m+4).descr = '';
nodes1(m+4).unic_name = "2+Layer1";
nodes1(m+4).labelSize = '0';

end

%quadranti
for iter = 1:1

nodes1(m+5).name = "5";
nodes1(m+5).labelPosition = "Bottom";
nodes1(m+5).labelAdjustX = '0';
nodes1(m+5).labelAdjustY = '0';
nodes1(m+5).isBold = "true";
nodes1(m+5).layer = 'Layer1';
nodes1(m+5).position_x = '0';
nodes1(m+5).position_y = "-500";
nodes1(m+5).position_z = "1";
nodes1(m+5).scale_x = "0";
nodes1(m+5).color = "#000000";
nodes1(m+5).url = '';
nodes1(m+5).descr = '';
nodes1(m+5).unic_name = "5+Layer1";
nodes1(m+5).labelSize = '0';

nodes1(m+6).name = "6";
nodes1(m+6).labelPosition = "Bottom";
nodes1(m+6).labelAdjustX = '0';
nodes1(m+6).labelAdjustY = '0';
nodes1(m+6).isBold = "true";
nodes1(m+6).layer = 'Layer1';
nodes1(m+6).position_x = '0';
nodes1(m+6).position_y = "500";
nodes1(m+6).position_z = "1";
nodes1(m+6).scale_x = "0";
nodes1(m+6).color = "#000000";
nodes1(m+6).url = '';
nodes1(m+6).descr = '';
nodes1(m+6).unic_name = "6+Layer1";
nodes1(m+6).labelSize = '0';

nodes1(m+7).name = "7";
nodes1(m+7).labelPosition = "Bottom";
nodes1(m+7).labelAdjustX = '0';
nodes1(m+7).labelAdjustY = '0';
nodes1(m+7).isBold = "true";
nodes1(m+7).layer = 'Layer1';
nodes1(m+7).position_x = '0';
nodes1(m+7).position_y = "0";
nodes1(m+7).position_z = "-500";
nodes1(m+7).scale_x = "0";
nodes1(m+7).color = "#000000";
nodes1(m+7).url = '';
nodes1(m+7).descr = '';
nodes1(m+7).unic_name = "7+Layer1";
nodes1(m+7).labelSize = '0';

nodes1(m+8).name = "8";
nodes1(m+8).labelPosition = "Bottom";
nodes1(m+8).labelAdjustX = '0';
nodes1(m+8).labelAdjustY = '0';
nodes1(m+8).isBold = "true";
nodes1(m+8).layer = 'Layer1';
nodes1(m+8).position_x = '0';
nodes1(m+8).position_y = "0";
nodes1(m+8).position_z = "500";
nodes1(m+8).scale_x = "0";
nodes1(m+8).color = "#000000";
nodes1(m+8).url = '';
nodes1(m+8).descr = '';
nodes1(m+8).unic_name = "8+Layer1";
nodes1(m+8).labelSize = '0';

end

% label quadranti
for iter = 1:1
nodes1(m+9).name = "Lack of Cognitive";
nodes1(m+9).labelPosition = "top";
nodes1(m+9).labelAdjustX = '0';
nodes1(m+9).labelAdjustY = '0';
nodes1(m+9).isBold = "true";
nodes1(m+9).layer = 'Layer1';
nodes1(m+9).position_x = '0';
nodes1(m+9).position_y = "-450";
nodes1(m+9).position_z = "-450";
nodes1(m+9).scale_x = "0";
nodes1(m+9).color = "#ffffff";
nodes1(m+9).url = '';
nodes1(m+9).descr = '';
nodes1(m+9).unic_name = "Lack of Cognitive+Layer1";
nodes1(m+9).labelSize = '20';

nodes1(m+10).name = "Lack_of Affective";
nodes1(m+10).labelPosition = "top";
nodes1(m+10).labelAdjustX = '0';
nodes1(m+10).labelAdjustY = '0';
nodes1(m+10).isBold = "true";
nodes1(m+10).layer = 'Layer1';
nodes1(m+10).position_x = '0';
nodes1(m+10).position_y = "-450";
nodes1(m+10).position_z = "180";
nodes1(m+10).scale_x = "0";
nodes1(m+10).color = "#ffffff";
nodes1(m+10).url = '';
nodes1(m+10).descr = '';
nodes1(m+10).unic_name = "Lack of Affective+Layer1";
nodes1(m+10).labelSize = '20';

nodes1(m+11).name = "Cognitive";
nodes1(m+11).labelPosition = "top";
nodes1(m+11).labelAdjustX = '0';
nodes1(m+11).labelAdjustY = '0';
nodes1(m+11).isBold = "true";
nodes1(m+11).layer = 'Layer1';
nodes1(m+11).position_x = '0';
nodes1(m+11).position_y = "450";
nodes1(m+11).position_z = "-450";
nodes1(m+11).scale_x = "0";
nodes1(m+11).color = "#ffffff";
nodes1(m+11).url = '';
nodes1(m+11).descr = '';
nodes1(m+11).unic_name = "Cognitive+Layer1";
nodes1(m+11).labelSize = '20';

nodes1(m+12).name = "Affective";
nodes1(m+12).labelPosition = "top";
nodes1(m+12).labelAdjustX = '0';
nodes1(m+12).labelAdjustY = '0';
nodes1(m+12).isBold = "true";
nodes1(m+12).layer = 'Layer1';
nodes1(m+12).position_x = '0';
nodes1(m+12).position_y = "450";
nodes1(m+12).position_z = "320";
nodes1(m+12).scale_x = "0";
nodes1(m+12).color = "#ffffff";
nodes1(m+12).url = '';
nodes1(m+12).descr = '';
nodes1(m+12).unic_name = "Affective+Layer1";
nodes1(m+12).labelSize = '20';

nodes1(m+13).name = "Well Being";
nodes1(m+13).labelPosition = "bottom";
nodes1(m+13).labelAdjustX = '0';
nodes1(m+13).labelAdjustY = '0';
nodes1(m+13).isBold = "true";
nodes1(m+13).layer = 'Layer1';
nodes1(m+13).position_x = '0';
nodes1(m+13).position_y = "-450";
nodes1(m+13).position_z = "-450";
nodes1(m+13).scale_x = "0";
nodes1(m+13).color = "#ffffff";
nodes1(m+13).url = '';
nodes1(m+13).descr = '';
nodes1(m+13).unic_name = "Well Being+Layer1";
nodes1(m+13).labelSize = '20';

nodes1(m+14).name = "Well Being";
nodes1(m+14).labelPosition = "bottom";
nodes1(m+14).labelAdjustX = '0';
nodes1(m+14).labelAdjustY = '0';
nodes1(m+14).isBold = "true";
nodes1(m+14).layer = 'Layer1';
nodes1(m+14).position_x = '0';
nodes1(m+14).position_y = "-450";
nodes1(m+14).position_z = "180";
nodes1(m+14).scale_x = "0";
nodes1(m+14).color = "#ffffff";
nodes1(m+14).url = '';
nodes1(m+14).descr = '';
nodes1(m+14).unic_name = "Well Being+Layer1";
nodes1(m+14).labelSize = '20';

nodes1(m+15).name = "Distress";
nodes1(m+15).labelPosition = "bottom";
nodes1(m+15).labelAdjustX = '0';
nodes1(m+15).labelAdjustY = '0';
nodes1(m+15).isBold = "true";
nodes1(m+15).layer = 'Layer1';
nodes1(m+15).position_x = '0';
nodes1(m+15).position_y = "450";
nodes1(m+15).position_z = "-450";
nodes1(m+15).scale_x = "0";
nodes1(m+15).color = "#ffffff";
nodes1(m+15).url = '';
nodes1(m+15).descr = '';
nodes1(m+15).unic_name = "Distress+Layer1";
nodes1(m+15).labelSize = '20';

nodes1(m+16).name = "Distress";
nodes1(m+16).labelPosition = "bottom";
nodes1(m+16).labelAdjustX = '0';
nodes1(m+16).labelAdjustY = '0';
nodes1(m+16).isBold = "true";
nodes1(m+16).layer = 'Layer1';
nodes1(m+16).position_x = '0';
nodes1(m+16).position_y = "450";
nodes1(m+16).position_z = "320";
nodes1(m+16).scale_x = "0";
nodes1(m+16).color = "#ffffff";
nodes1(m+16).url = '';
nodes1(m+16).descr = '';
nodes1(m+16).unic_name = "Distress+Layer1";
nodes1(m+16).labelSize = '20';
end



end

%intelaiatura L2
for iter = 1:1


m = 0;
m = max(size(nodes_unic_names_layer_2));

%vertex
for iter = 1:1

nodes2(m+1).name = "4";
nodes2(m+1).labelPosition = "Top";
nodes2(m+1).labelAdjustX = '0';
nodes2(m+1).labelAdjustY = '0';
nodes2(m+1).isBold = "true";
nodes2(m+1).layer = 'Layer2';
nodes2(m+1).position_x = '0';
nodes2(m+1).position_y = "500";
nodes2(m+1).position_z = "500";
nodes2(m+1).scale_x = "0";
nodes2(m+1).color = "#000000";
nodes2(m+1).url = '';
nodes2(m+1).descr = '';
nodes2(m+1).unic_name = "4+Layer2";
nodes2(m+1).labelSize = '0';

nodes2(m+2).name = "3";
nodes2(m+2).labelPosition = "Top";
nodes2(m+2).labelAdjustX = '0';
nodes2(m+2).labelAdjustY = '0';
nodes2(m+2).isBold = "true";
nodes2(m+2).layer = 'Layer2';
nodes2(m+2).position_x = '0';
nodes2(m+2).position_y = "500";
nodes2(m+2).position_z = "-500";
nodes2(m+2).scale_x = "0";
nodes2(m+2).color = "#000000";
nodes2(m+2).url = '';
nodes2(m+2).descr = '';
nodes2(m+2).unic_name = "3+Layer2";
nodes2(m+2).labelSize = '0';

nodes2(m+3).name = "1";
nodes2(m+3).labelPosition = "Bottom";
nodes2(m+3).labelAdjustX = '0';
nodes2(m+3).labelAdjustY = '0';
nodes2(m+3).isBold = "true";
nodes2(m+3).layer = 'Layer2';
nodes2(m+3).position_x = '0';
nodes2(m+3).position_y = "-500";
nodes2(m+3).position_z = "-500";
nodes2(m+3).scale_x = "0";
nodes2(m+3).color = "#000000";
nodes2(m+3).url = '';
nodes2(m+3).descr = '';
nodes2(m+3).unic_name = "1+Layer2";
nodes2(m+3).labelSize = '0';

nodes2(m+4).name = "2";
nodes2(m+4).labelPosition = "Bottom";
nodes2(m+4).labelAdjustX = '0';
nodes2(m+4).labelAdjustY = '0';
nodes2(m+4).isBold = "true";
nodes2(m+4).layer = 'Layer2';
nodes2(m+4).position_x = '0';
nodes2(m+4).position_y = "-500";
nodes2(m+4).position_z = "500";
nodes2(m+4).scale_x = "0";
nodes2(m+4).color = "#000000";
nodes2(m+4).url = '';
nodes2(m+4).descr = '';
nodes2(m+4).unic_name = "2+Layer2";
nodes2(m+4).labelSize = '0';

end

%quadranti
for iter = 1:1

nodes2(m+5).name = "5";
nodes2(m+5).labelPosition = "Bottom";
nodes2(m+5).labelAdjustX = '0';
nodes2(m+5).labelAdjustY = '0';
nodes2(m+5).isBold = "true";
nodes2(m+5).layer = 'Layer2';
nodes2(m+5).position_x = '0';
nodes2(m+5).position_y = "-500";
nodes2(m+5).position_z = "1";
nodes2(m+5).scale_x = "0";
nodes2(m+5).color = "#000000";
nodes2(m+5).url = '';
nodes2(m+5).descr = '';
nodes2(m+5).unic_name = "5+Layer2";
nodes2(m+5).labelSize = '0';

nodes2(m+6).name = "6";
nodes2(m+6).labelPosition = "Bottom";
nodes2(m+6).labelAdjustX = '0';
nodes2(m+6).labelAdjustY = '0';
nodes2(m+6).isBold = "true";
nodes2(m+6).layer = 'Layer2';
nodes2(m+6).position_x = '0';
nodes2(m+6).position_y = "500";
nodes2(m+6).position_z = "1";
nodes2(m+6).scale_x = "0";
nodes2(m+6).color = "#000000";
nodes2(m+6).url = '';
nodes2(m+6).descr = '';
nodes2(m+6).unic_name = "6+Layer2";
nodes2(m+6).labelSize = '0';

nodes2(m+7).name = "7";
nodes2(m+7).labelPosition = "Bottom";
nodes2(m+7).labelAdjustX = '0';
nodes2(m+7).labelAdjustY = '0';
nodes2(m+7).isBold = "true";
nodes2(m+7).layer = 'Layer2';
nodes2(m+7).position_x = '0';
nodes2(m+7).position_y = "0";
nodes2(m+7).position_z = "-500";
nodes2(m+7).scale_x = "0";
nodes2(m+7).color = "#000000";
nodes2(m+7).url = '';
nodes2(m+7).descr = '';
nodes2(m+7).unic_name = "7+Layer2";
nodes2(m+7).labelSize = '0';

nodes2(m+8).name = "8";
nodes2(m+8).labelPosition = "Bottom";
nodes2(m+8).labelAdjustX = '0';
nodes2(m+8).labelAdjustY = '0';
nodes2(m+8).isBold = "true";
nodes2(m+8).layer = 'Layer2';
nodes2(m+8).position_x = '0';
nodes2(m+8).position_y = "0";
nodes2(m+8).position_z = "500";
nodes2(m+8).scale_x = "0";
nodes2(m+8).color = "#000000";
nodes2(m+8).url = '';
nodes2(m+8).descr = '';
nodes2(m+8).unic_name = "8+Layer2";
nodes2(m+8).labelSize = '0';

end

% label quadranti
for iter = 1:1
nodes2(m+9).name = "Lack_of Cognitive";
nodes2(m+9).labelPosition = "top";
nodes2(m+9).labelAdjustX = '0';
nodes2(m+9).labelAdjustY = '0';
nodes2(m+9).isBold = "true";
nodes2(m+9).layer = 'Layer2';
nodes2(m+9).position_x = '0';
nodes2(m+9).position_y = "-450";
nodes2(m+9).position_z = "-450";
nodes2(m+9).scale_x = "0";
nodes2(m+9).color = "#ffffff";
nodes2(m+9).url = '';
nodes2(m+9).descr = '';
nodes2(m+9).unic_name = "Lack of Cognitive+Layer2";
nodes2(m+9).labelSize = '20';

nodes2(m+10).name = "Lack of Affective";
nodes2(m+10).labelPosition = "top";
nodes2(m+10).labelAdjustX = '0';
nodes2(m+10).labelAdjustY = '0';
nodes2(m+10).isBold = "true";
nodes2(m+10).layer = 'Layer2';
nodes2(m+10).position_x = '0';
nodes2(m+10).position_y = "-450";
nodes2(m+10).position_z = "180";
nodes2(m+10).scale_x = "0";
nodes2(m+10).color = "#ffffff";
nodes2(m+10).url = '';
nodes2(m+10).descr = '';
nodes2(m+10).unic_name = "Lack of Affective+Layer2";
nodes2(m+10).labelSize = '20';

nodes2(m+11).name = "Cognitive";
nodes2(m+11).labelPosition = "top";
nodes2(m+11).labelAdjustX = '0';
nodes2(m+11).labelAdjustY = '0';
nodes2(m+11).isBold = "true";
nodes2(m+11).layer = 'Layer2';
nodes2(m+11).position_x = '0';
nodes2(m+11).position_y = "450";
nodes2(m+11).position_z = "-450";
nodes2(m+11).scale_x = "0";
nodes2(m+11).color = "#ffffff";
nodes2(m+11).url = '';
nodes2(m+11).descr = '';
nodes2(m+11).unic_name = "Cognitive+Layer2";
nodes2(m+11).labelSize = '20';

nodes2(m+12).name = "Affective";
nodes2(m+12).labelPosition = "top";
nodes2(m+12).labelAdjustX = '0';
nodes2(m+12).labelAdjustY = '0';
nodes2(m+12).isBold = "true";
nodes2(m+12).layer = 'Layer2';
nodes2(m+12).position_x = '0';
nodes2(m+12).position_y = "450";
nodes2(m+12).position_z = "320";
nodes2(m+12).scale_x = "0";
nodes2(m+12).color = "#ffffff";
nodes2(m+12).url = '';
nodes2(m+12).descr = '';
nodes2(m+12).unic_name = "Affective+Layer2";
nodes2(m+12).labelSize = '20';

nodes2(m+13).name = "Well Being";
nodes2(m+13).labelPosition = "bottom";
nodes2(m+13).labelAdjustX = '0';
nodes2(m+13).labelAdjustY = '0';
nodes2(m+13).isBold = "true";
nodes2(m+13).layer = 'Layer2';
nodes2(m+13).position_x = '0';
nodes2(m+13).position_y = "-450";
nodes2(m+13).position_z = "-450";
nodes2(m+13).scale_x = "0";
nodes2(m+13).color = "#ffffff";
nodes2(m+13).url = '';
nodes2(m+13).descr = '';
nodes2(m+13).unic_name = "Well Being+Layer2";
nodes2(m+13).labelSize = '20';

nodes2(m+14).name = "Well Being";
nodes2(m+14).labelPosition = "bottom";
nodes2(m+14).labelAdjustX = '0';
nodes2(m+14).labelAdjustY = '0';
nodes2(m+14).isBold = "true";
nodes2(m+14).layer = 'Layer2';
nodes2(m+14).position_x = '0';
nodes2(m+14).position_y = "-450";
nodes2(m+14).position_z = "180";
nodes2(m+14).scale_x = "0";
nodes2(m+14).color = "#ffffff";
nodes2(m+14).url = '';
nodes2(m+14).descr = '';
nodes2(m+14).unic_name = "Well Being+Layer2";
nodes2(m+14).labelSize = '20';

nodes2(m+15).name = "Distress";
nodes2(m+15).labelPosition = "bottom";
nodes2(m+15).labelAdjustX = '0';
nodes2(m+15).labelAdjustY = '0';
nodes2(m+15).isBold = "true";
nodes2(m+15).layer = 'Layer2';
nodes2(m+15).position_x = '0';
nodes2(m+15).position_y = "450";
nodes2(m+15).position_z = "-450";
nodes2(m+15).scale_x = "0";
nodes2(m+15).color = "#ffffff";
nodes2(m+15).url = '';
nodes2(m+15).descr = '';
nodes2(m+15).unic_name = "Distress+Layer2";
nodes2(m+15).labelSize = '20';

nodes2(m+16).name = "Distress";
nodes2(m+16).labelPosition = "bottom";
nodes2(m+16).labelAdjustX = '0';
nodes2(m+16).labelAdjustY = '0';
nodes2(m+16).isBold = "true";
nodes2(m+16).layer = 'Layer2';
nodes2(m+16).position_x = '0';
nodes2(m+16).position_y = "450";
nodes2(m+16).position_z = "320";
nodes2(m+16).scale_x = "0";
nodes2(m+16).color = "#ffffff";
nodes2(m+16).url = '';
nodes2(m+16).descr = '';
nodes2(m+16).unic_name = "Distress  +Layer2";
nodes2(m+16).labelSize = '20';
end


end




%value4 = [nodes1,nodes2];

end


%% field 5 edges
for iter = 1:1
    
field5 = 'edges';


% edges 1 %%%%%%%%%
edges_1 = [];
for iter = 1:max(size(src_names_edges_cross_1))
    edges_1(iter).src = src_names_edges_cross_1(iter);
    edges_1(iter).trg = trg_names_edges_cross_1(iter);
    edges_1(iter).opacity = '0.3';
    edges_1(iter).color = COLORS1(iter);
    edges_1(iter).channel = '';
    %edges_1(iter).size = WEIGHT_ARENA1(iter)/5; %string(cc1(iter)*4);
    edges_1(iter).size = log(WEIGHT_ARENA1(iter)*2)*3;
    edges_1(iter).arrow = "true";
end


m = 0;
m = size((src_names_edges_cross_1),1);

edges_1(m+1).src = "1+Layer1";
edges_1(m+1).trg = "2+Layer1";
edges_1(m+1).opacity = '0.3';
edges_1(m+1).color = color_quadrants;  % #b9bdc4 grigio
edges_1(m+1).channel = '';
edges_1(m+1).size = "2";
edges_1(m+1).arrow = "false";

edges_1(m+2).src = "1+Layer1";
edges_1(m+2).trg = "3+Layer1";
edges_1(m+2).opacity = '0.3';
edges_1(m+2).color = color_quadrants;
edges_1(m+2).channel = '';
edges_1(m+2).size = "2";
edges_1(m+2).arrow = "false";

edges_1(m+3).src = "2+Layer1";
edges_1(m+3).trg = "4+Layer1";
edges_1(m+3).opacity = '0.3';
edges_1(m+3).color = color_quadrants;
edges_1(m+3).channel = '';
edges_1(m+3).size = "2";
edges_1(m+3).arrow = "false";

edges_1(m+4).src = "3+Layer1";
edges_1(m+4).trg = "4+Layer1";
edges_1(m+4).opacity = '0.3';
edges_1(m+4).color = color_quadrants;
edges_1(m+4).channel = '';
edges_1(m+4).size = "2";
edges_1(m+4).arrow = "false";

edges_1(m+5).src = "6+Layer1";
edges_1(m+5).trg = "5+Layer1";
edges_1(m+5).opacity = '0.3';
edges_1(m+5).color = "#000000";
edges_1(m+5).channel = '';
edges_1(m+5).size = "1";
edges_1(m+5).arrow = "false";

edges_1(m+6).src = "7+Layer1";
edges_1(m+6).trg = "8+Layer1";
edges_1(m+6).opacity = '0.3';
edges_1(m+6).color = "#000000";
edges_1(m+6).channel = '';
edges_1(m+6).size = "1";
edges_1(m+6).arrow = "false";


% edges 1_2 %%%%%%%%%
edges_1_2 = [];
for iter = 1:max(size(src_names_edges_longi_1_2))
    edges_1_2(iter).src = src_names_edges_longi_1_2(iter);
    edges_1_2(iter).trg = trg_names_edges_longi_1_2(iter);
    edges_1_2(iter).opacity = '0.3';
    edges_1_2(iter).color = COLORS1_2(iter);
    edges_1_2(iter).channel = '';
    %edges_1_2(iter).size = WEIGHT_ARENA1_2(iter)/5; %string(cc1_2(iter)*4);
    edges_1_2(iter).size = log(WEIGHT_ARENA1_2(iter)*2)*3; 
    edges_1_2(iter).arrow = "true";
end



% edges 2 %%%%%%%%%
edges_2 = [];
for iter = 1:max(size(src_names_edges_cross_2))
    edges_2(iter).src = src_names_edges_cross_2(iter);
    edges_2(iter).trg = trg_names_edges_cross_2(iter);
    edges_2(iter).opacity = '0.3';
    edges_2(iter).color = COLORS2(iter);
    edges_2(iter).channel = '';
    %edges_2(iter).size = WEIGHT_ARENA2(iter)/5; %string(cc2(iter)*4);
    edges_2(iter).size = log(WEIGHT_ARENA2(iter)*2)*3;
    edges_2(iter).arrow = "true";
end


m = 0;
m = size((src_names_edges_cross_2),1);

edges_2(m+1).src = "1+Layer2";
edges_2(m+1).trg = "2+Layer2";
edges_2(m+1).opacity = '0.3';
edges_2(m+1).color = color_quadrants;
edges_2(m+1).channel = '';
edges_2(m+1).size = "2";
edges_2(m+1).arrow = "false";

edges_2(m+2).src = "1+Layer2";
edges_2(m+2).trg = "3+Layer2";
edges_2(m+2).opacity = '0.3';
edges_2(m+2).color = color_quadrants;
edges_2(m+2).channel = '';
edges_2(m+2).size = "2";
edges_2(m+2).arrow = "false";

edges_2(m+3).src = "2+Layer2";
edges_2(m+3).trg = "4+Layer2";
edges_2(m+3).opacity = '0.3';
edges_2(m+3).color = color_quadrants;
edges_2(m+3).channel = '';
edges_2(m+3).size = "2";
edges_2(m+3).arrow = "false";

edges_2(m+4).src = "3+Layer2";
edges_2(m+4).trg = "4+Layer2";
edges_2(m+4).opacity = '0.3';
edges_2(m+4).color = color_quadrants;
edges_2(m+4).channel = '';
edges_2(m+4).size = "2";
edges_2(m+4).arrow = "false";

edges_2(m+5).src = "5+Layer2";
edges_2(m+5).trg = "6+Layer2";
edges_2(m+5).opacity = '0.3';
edges_2(m+5).color = "#000000";
edges_2(m+5).channel = '';
edges_2(m+5).size = "1";
edges_2(m+5).arrow = "false";

edges_2(m+6).src = "7+Layer2";
edges_2(m+6).trg = "8+Layer2";
edges_2(m+6).opacity = '0.3';
edges_2(m+6).color = "#000000";
edges_2(m+6).channel = '';
edges_2(m+6).size = "1";
edges_2(m+6).arrow = "false";

    
%value5 = [edges_1,edges_1_2,edges_2];


end


%% average vector
for iter = 1:1

if plot_avg_vector == 1

label_size_vect = '10';
size_avgv = '5';
size_edge_avgv = '12';

%node
for iter = 1:1
n1 = numel(nodes1);

nodes1(n1+1).name = 'entry_t1';
nodes1(n1+1).labelPosition = 'top';
nodes1(n1+1).labelAdjustX = '0';
nodes1(n1+1).labelAdjustY = '0';
nodes1(n1+1).isBold = 'true';
nodes1(n1+1).layer = 'Layer1';
nodes1(n1+1).position_x = '0';
nodes1(n1+1).position_y = string(coord_avgv_entry_t1(1,2)*1000);
nodes1(n1+1).position_z = string(coord_avgv_entry_t1(1,1)*1000);
nodes1(n1+1).scale_x = size_avgv;
nodes1(n1+1).color = color_t1;
nodes1(n1+1).url = '';
nodes1(n1+1).descr = '';
nodes1(n1+1).unic_name = 'entry_t1+layer1';
nodes1(n1+1).labelSize = label_size_vect;

nodes1(n1+2).name = 'exit_t1';
nodes1(n1+2).labelPosition = 'top';
nodes1(n1+2).labelAdjustX = '0';
nodes1(n1+2).labelAdjustY = '0';
nodes1(n1+2).isBold = 'true';
nodes1(n1+2).layer = 'Layer1';
nodes1(n1+2).position_x = '0';
nodes1(n1+2).position_y = string(coord_avgv_exit_t1(1,2)*1000);
nodes1(n1+2).position_z = string(coord_avgv_exit_t1(1,1)*1000);
nodes1(n1+2).scale_x = size_avgv;
nodes1(n1+2).color = color_longi;
nodes1(n1+2).url = '';
nodes1(n1+2).descr = '';
nodes1(n1+2).unic_name = 'exit_t1+layer1';
nodes1(n1+2).labelSize = label_size_vect;

nodes1(n1+3).name = 'entry_t2';
nodes1(n1+3).labelPosition = 'top';
nodes1(n1+3).labelAdjustX = '0';
nodes1(n1+3).labelAdjustY = '0';
nodes1(n1+3).isBold = 'true';
nodes1(n1+3).layer = 'Layer2';
nodes1(n1+3).position_x = '0';
nodes1(n1+3).position_y = string(coord_avgv_entry_t2(1,2)*1000);
nodes1(n1+3).position_z = string(coord_avgv_entry_t2(1,1)*1000);
nodes1(n1+3).scale_x = size_avgv;
nodes1(n1+3).color = color_longi;
nodes1(n1+3).url = '';
nodes1(n1+3).descr = '';
nodes1(n1+3).unic_name = 'entry_t2+layer1';
nodes1(n1+3).labelSize = label_size_vect;

nodes1(n1+4).name = 'exit_t2';
nodes1(n1+4).labelPosition = 'top';
nodes1(n1+4).labelAdjustX = '0';
nodes1(n1+4).labelAdjustY = '0';
nodes1(n1+4).isBold = 'true';
nodes1(n1+4).layer = 'Layer2';
nodes1(n1+4).position_x = '0';
nodes1(n1+4).position_y = string(coord_avgv_exit_t2(1,2)*1000);
nodes1(n1+4).position_z = string(coord_avgv_exit_t2(1,1)*1000);
nodes1(n1+4).scale_x = size_avgv;
nodes1(n1+4).color = color_t2;
nodes1(n1+4).url = '';
nodes1(n1+4).descr = '';
nodes1(n1+4).unic_name = 'exit_t2+layer1';
nodes1(n1+4).labelSize = label_size_vect;

end

% edges
for iter = 1:1
n1 = numel(edges_1);

edges_1(n1+1).src = "entry_t1+Layer1";
edges_1(n1+1).trg = "exit_t1+Layer1";
edges_1(n1+1).opacity = '0.3';
edges_1(n1+1).color = color_t1;
edges_1(n1+1).channel = '';
edges_1(n1+1).size = size_edge_avgv;
edges_1(n1+1).arrow = "true";

edges_1(n1+2).src = "exit_t1+Layer1";
edges_1(n1+2).trg = "entry_t2+Layer2";
edges_1(n1+2).opacity = '0.3';
edges_1(n1+2).color = color_longi;
edges_1(n1+2).channel = '';
edges_1(n1+2).size = size_edge_avgv;
edges_1(n1+2).arrow = "true";

edges_1(n1+3).src = "entry_t2+Layer2";
edges_1(n1+3).trg = "exit_t2+Layer2";
edges_1(n1+3).opacity = '0.3';
edges_1(n1+3).color = color_t2;
edges_1(n1+3).channel = '';
edges_1(n1+3).size = size_edge_avgv;
edges_1(n1+3).arrow = "true";


end

end


end


value4 = [nodes1,nodes2];
value5 = [edges_1,edges_1_2,edges_2];



%% def struct
clear s;
s = struct(field3,value3,field4,value4,field5,value5);
s1 = s;


end


%% Create structure to plot shortpaths net (s2)
% same procedure used above. The difference is that here the structure is created for the network with all the shortest paths.
for iter = 1:1
    


%% input data for the export
for iter=1:1


%% adj matrix
for iter = 1:max(size(net1))
    net1(iter,iter)= 0;
end

net1(TARGET,SOURCE)=0;
D=1./posweights(net1);
[SPL,hops,Pmat] = distance_wei_floyd(D);

ACROSS_BETW=zeros(max(size(net1)),max(size(net1)));
BTW = zeros(max(size(net1)),1); 
nc = zeros(max(size(net1)),1);


for iter1=1:max(size(SOURCE))
    for iter2=1:max(size(TARGET))
        B=retrieve_shortest_path(SOURCE(iter1),TARGET(iter2),hops,Pmat);
        if min(size(intersect(B,hub_to_plot)))>0
            %adj mat
            ACROSS_BETW2=zeros(max(size(net1)),max(size(net1)));
            for iter3=1:(max(size(B))-1)
                ACROSS_BETW2(B(iter3),B(iter3+1))=1;             
            end
            ACROSS_BETW=ACROSS_BETW+ACROSS_BETW2;
            
            %btw
            BTW2=zeros(max(size(net1)),1);
            BTW2(B)=1;
            BTW=BTW+BTW2;
        end
    end
end


%% node

% node size
for iter = 1:dim_cross
    node_size1(iter) = log(BTW(iter));
    node_size2(iter) = log(BTW(iter+dim_cross));
end

for iter = 1:dim_cross
    node_size1(iter) = log(BTW(iter));
    node_size2(iter) = log(BTW(iter+dim_cross));
end

% node size threshold
for iter = 1:dim_cross
    if node_size1(iter) <= 1
        node_size1(iter) = 1;
    end
    if node_size2(iter) <= 1
        node_size2(iter) = 1;
    end
end


% node unic names
count = 1;
for iter=1:20
        nodes_unic_names_layer_1(count) = strcat(names_3(iter),'+Layer1');
        nodes_unic_names_layer_2(count) = strcat(names_3(iter),'+Layer2');
        count = count + 1;
end
nodes_unic_names_layer_1 = string(nodes_unic_names_layer_1)';
nodes_unic_names_layer_2 = string(nodes_unic_names_layer_2)';


% node color
for iter = 1:max(size(nc))
    if nc(iter) == 1
        nodecolors(iter,:) = color_t1;
    end
    if nc(iter) == 2
        nodecolors(iter,:) = color_longi;
    end
    if nc(iter) == 3
        nodecolors(iter,:) = color_longi;
    end
    if nc(iter) == 4
        nodecolors(iter,:) = color_t2;
    end
    if nc(iter) == 0
        nodecolors(iter,:) = "#E0E0E0";
    end
end



%% edges
% names
NAMES_NEW2 = names_3;

count = 1;
for iter=1:dim_cross
    for iter2=1:dim_cross
        src_names_edges_cross_1(count) = strcat(NAMES_NEW2(iter),'+Layer1');
        trg_names_edges_cross_1(count) = strcat(NAMES_NEW2(iter2),'+Layer1');
        count = count + 1;
    end
end
src_names_edges_cross_1 = src_names_edges_cross_1';
src_names_edges_cross_1 = string(src_names_edges_cross_1);
trg_names_edges_cross_1 = trg_names_edges_cross_1';
trg_names_edges_cross_1 = string(trg_names_edges_cross_1);


count = 1;
for iter=1:dim_cross
    for iter2=dim_cross+1:dim_longi
        src_names_edges_longi_1_2(count) = strcat(NAMES_NEW2(iter),'+Layer1');
        trg_names_edges_longi_1_2(count) = strcat(NAMES_NEW2(iter2-dim_cross),'+Layer2');
        count = count + 1;
    end
end
src_names_edges_longi_1_2 = src_names_edges_longi_1_2';
src_names_edges_longi_1_2 = string(src_names_edges_longi_1_2);
trg_names_edges_longi_1_2 = trg_names_edges_longi_1_2';
trg_names_edges_longi_1_2 = string(trg_names_edges_longi_1_2);

count=1;
for iter=1:dim_cross
    for iter2=1:dim_cross
        src_names_edges_cross_2(count) = strcat(NAMES_NEW2(iter),'+Layer2');
        trg_names_edges_cross_2(count) = strcat(NAMES_NEW2(iter2),'+Layer2');
        count = count + 1;
    end
end
src_names_edges_cross_2 = src_names_edges_cross_2';
src_names_edges_cross_2 = string(src_names_edges_cross_2);
trg_names_edges_cross_2 = trg_names_edges_cross_2';
trg_names_edges_cross_2 = string(trg_names_edges_cross_2);




% weights
W = ACROSS_BETW;


clear WEIGHT_ARENA1
clear COLORS1
count=1;
for iter=1:dim_cross
    for iter2=1:dim_cross
        WEIGHT_ARENA1(count) = W(iter,iter2);
        COLORS1(count) = color_allnet;
        count = count + 1;
    end
end
WEIGHT_ARENA1 = WEIGHT_ARENA1';

clear WEIGHT_ARENA1_2
clear COLORS1_2
count=1;
for iter=1:dim_cross
    for iter2=dim_cross+1:dim_longi
        WEIGHT_ARENA1_2(count) = W(iter,iter2);
        COLORS1_2(count) = color_allnet;
        count = count + 1;
    end
end
WEIGHT_ARENA1_2 = WEIGHT_ARENA1_2';

clear WEIGHT_ARENA2
clear COLORS2
count=1;
for iter=dim_cross+1:dim_longi
    for iter2=dim_cross+1:dim_longi
        WEIGHT_ARENA2(count) = W(iter,iter2);
        COLORS2(count) = color_allnet;
        count = count + 1;
    end
end
WEIGHT_ARENA2 = WEIGHT_ARENA2';

clear a b c d;



% delete 0 or negative corr (negative should not be here)
src_names_edges_cross_1(WEIGHT_ARENA1==0)=[];
trg_names_edges_cross_1(WEIGHT_ARENA1==0)=[];
COLORS1(WEIGHT_ARENA1==0)=[];
WEIGHT_ARENA1(WEIGHT_ARENA1==0)=[];

src_names_edges_cross_1(WEIGHT_ARENA1<0)=[];
trg_names_edges_cross_1(WEIGHT_ARENA1<0)=[];
COLORS1(WEIGHT_ARENA1<0)=[];
WEIGHT_ARENA1(WEIGHT_ARENA1<0)=[];

%%%%%%%%
src_names_edges_longi_1_2(WEIGHT_ARENA1_2==0)=[];
trg_names_edges_longi_1_2(WEIGHT_ARENA1_2==0)=[];
COLORS1_2(WEIGHT_ARENA1_2==0)=[];
WEIGHT_ARENA1_2(WEIGHT_ARENA1_2==0)=[];

src_names_edges_longi_1_2(WEIGHT_ARENA1_2<0)=[];
trg_names_edges_longi_1_2(WEIGHT_ARENA1_2<0)=[];
COLORS1_2(WEIGHT_ARENA1_2<0)=[];
WEIGHT_ARENA1_2(WEIGHT_ARENA1_2<0)=[];

%%%%%%%%
src_names_edges_cross_2(WEIGHT_ARENA2==0)=[];
trg_names_edges_cross_2(WEIGHT_ARENA2==0)=[];
COLORS2(WEIGHT_ARENA2==0)=[];
WEIGHT_ARENA2(WEIGHT_ARENA2==0)=[];

src_names_edges_cross_2(WEIGHT_ARENA2<0)=[];
trg_names_edges_cross_2(WEIGHT_ARENA2<0)=[];
COLORS2(WEIGHT_ARENA2<0)=[];
WEIGHT_ARENA2(WEIGHT_ARENA2<0)=[];


end



%% field 3 layers
for iter = 1:1

field3 = 'layers';

c1.name = 'Layer1';
c1.position_x = '-300';
c1.position_y = '0';
c1.position_z = '0';
c1.last_layer_scale = '1';
c1.rotation_x = '0';
c1.rotation_y = '0';
c1.rotation_z = '0';
c1.floor_current_color = '#777777';
c1.geometry_parameters_width = '1001.90476190476';


c2.name = 'Layer2';
c2.position_x = '300';
c2.position_y = '0';
c2.position_z = '0';
c2.last_layer_scale = '1';
c2.rotation_x = '0';
c2.rotation_y = '0';
c2.rotation_z = '0';
c2.floor_current_color = '#777777';
c2.geometry_parameters_width = '1001.90476190476';



value3 = struct([c1,c2]);

end


%% field 4 nodes
for iter = 1:1

field4 = 'nodes';
nodes1 = [];
nodes2 = [];

% nodes 1 %%%%%%%%%
for iter = 1:dim_cross
    nodes1(iter).name = string(names_3(iter));
    nodes1(iter).labelPosition = "Top";
    nodes1(iter).labelAdjustX = '0';
    nodes1(iter).labelAdjustY = '0';
    nodes1(iter).isBold = "true";
    %nodes1(iter).isBold = bold_no_corr_str(iter);
    nodes1(iter).layer = 'Layer1';
    nodes1(iter).position_x = '0';
    nodes1(iter).position_y = string(XY_CROSS(iter,2)*1000);
    nodes1(iter).position_z = string(XY_CROSS(iter,1)*1000);
    %nodes1(iter).scale_x = string(BTW_size(iter)/25);
    nodes1(iter).scale_x = node_size1(iter);
    nodes1(iter).color = nodecolors(iter);
    nodes1(iter).url = '';
    nodes1(iter).descr = '';
    nodes1(iter).unic_name = nodes_unic_names_layer_1(iter);
    nodes1(iter).labelSize = '8';
end


% nodes 2 %%%%%%%%%
for iter = 1:dim_cross
    nodes2(iter).name = string(names_3(iter));
    nodes2(iter).labelPosition = "Top";
    nodes2(iter).labelAdjustX = '0';
    nodes2(iter).labelAdjustY = '0';
    nodes2(iter).isBold = "true";
    %nodes2(iter).isBold = bold_no_corr_str(iter+20);
    nodes2(iter).layer = 'Layer2';
    nodes2(iter).position_x = '0';
    nodes2(iter).position_y = string(XY_CROSS(iter,2)*1000);
    nodes2(iter).position_z = string(XY_CROSS(iter,1)*1000);
    %nodes2(iter).scale_x = string(BTW_size(iter+20)/25);
    nodes2(iter).scale_x = node_size2(iter);
    nodes2(iter).color = nodecolors(iter+dim_cross);
    nodes2(iter).url = '';
    nodes2(iter).descr = '';
    nodes2(iter).unic_name = nodes_unic_names_layer_2(iter);
    nodes2(iter).labelSize = '8';
end




%intelaiatura L1
for iter = 1:1


m = 0;
m = max(size(nodes_unic_names_layer_1));

%vertex
for iter = 1:1

nodes1(m+1).name = "4";
nodes1(m+1).labelPosition = "Top";
nodes1(m+1).labelAdjustX = '0';
nodes1(m+1).labelAdjustY = '0';
nodes1(m+1).isBold = "true";
nodes1(m+1).layer = 'Layer1';
nodes1(m+1).position_x = '0';
nodes1(m+1).position_y = "500";
nodes1(m+1).position_z = "500";
nodes1(m+1).scale_x = "0";
nodes1(m+1).color = "#000000";
nodes1(m+1).url = '';
nodes1(m+1).descr = '';
nodes1(m+1).unic_name = "4+Layer1";
nodes1(m+1).labelSize = '';

nodes1(m+2).name = "3";
nodes1(m+2).labelPosition = "Top";
nodes1(m+2).labelAdjustX = '0';
nodes1(m+2).labelAdjustY = '0';
nodes1(m+2).isBold = "true";
nodes1(m+2).layer = 'Layer1';
nodes1(m+2).position_x = '0';
nodes1(m+2).position_y = "500";
nodes1(m+2).position_z = "-500";
nodes1(m+2).scale_x = "0";
nodes1(m+2).color = "#000000";
nodes1(m+2).url = '';
nodes1(m+2).descr = '';
nodes1(m+2).unic_name = "3+Layer1";
nodes1(m+2).labelSize = '0';

nodes1(m+3).name = "1";
nodes1(m+3).labelPosition = "Top";
nodes1(m+3).labelAdjustX = '0';
nodes1(m+3).labelAdjustY = '0';
nodes1(m+3).isBold = "true";
nodes1(m+3).layer = 'Layer1';
nodes1(m+3).position_x = '0';
nodes1(m+3).position_y = "-500";
nodes1(m+3).position_z = "-500";
nodes1(m+3).scale_x = "0";
nodes1(m+3).color = "#000000";
nodes1(m+3).url = '';
nodes1(m+3).descr = '';
nodes1(m+3).unic_name = "1+Layer1";
nodes1(m+3).labelSize = '0';

nodes1(m+4).name = "2";
nodes1(m+4).labelPosition = "Top";
nodes1(m+4).labelAdjustX = '0';
nodes1(m+4).labelAdjustY = '0';
nodes1(m+4).isBold = "true";
nodes1(m+4).layer = 'Layer1';
nodes1(m+4).position_x = '0';
nodes1(m+4).position_y = "-500";
nodes1(m+4).position_z = "500";
nodes1(m+4).scale_x = "0";
nodes1(m+4).color = "#000000";
nodes1(m+4).url = '';
nodes1(m+4).descr = '';
nodes1(m+4).unic_name = "2+Layer1";
nodes1(m+4).labelSize = '0';

end

%quadranti
for iter = 1:1

nodes1(m+5).name = "5";
nodes1(m+5).labelPosition = "Bottom";
nodes1(m+5).labelAdjustX = '0';
nodes1(m+5).labelAdjustY = '0';
nodes1(m+5).isBold = "true";
nodes1(m+5).layer = 'Layer1';
nodes1(m+5).position_x = '0';
nodes1(m+5).position_y = "-500";
nodes1(m+5).position_z = "1";
nodes1(m+5).scale_x = "0";
nodes1(m+5).color = "#000000";
nodes1(m+5).url = '';
nodes1(m+5).descr = '';
nodes1(m+5).unic_name = "5+Layer1";
nodes1(m+5).labelSize = '0';

nodes1(m+6).name = "6";
nodes1(m+6).labelPosition = "Bottom";
nodes1(m+6).labelAdjustX = '0';
nodes1(m+6).labelAdjustY = '0';
nodes1(m+6).isBold = "true";
nodes1(m+6).layer = 'Layer1';
nodes1(m+6).position_x = '0';
nodes1(m+6).position_y = "500";
nodes1(m+6).position_z = "1";
nodes1(m+6).scale_x = "0";
nodes1(m+6).color = "#000000";
nodes1(m+6).url = '';
nodes1(m+6).descr = '';
nodes1(m+6).unic_name = "6+Layer1";
nodes1(m+6).labelSize = '0';

nodes1(m+7).name = "7";
nodes1(m+7).labelPosition = "Bottom";
nodes1(m+7).labelAdjustX = '0';
nodes1(m+7).labelAdjustY = '0';
nodes1(m+7).isBold = "true";
nodes1(m+7).layer = 'Layer1';
nodes1(m+7).position_x = '0';
nodes1(m+7).position_y = "0";
nodes1(m+7).position_z = "-500";
nodes1(m+7).scale_x = "0";
nodes1(m+7).color = "#000000";
nodes1(m+7).url = '';
nodes1(m+7).descr = '';
nodes1(m+7).unic_name = "7+Layer1";
nodes1(m+7).labelSize = '0';

nodes1(m+8).name = "8";
nodes1(m+8).labelPosition = "Bottom";
nodes1(m+8).labelAdjustX = '0';
nodes1(m+8).labelAdjustY = '0';
nodes1(m+8).isBold = "true";
nodes1(m+8).layer = 'Layer1';
nodes1(m+8).position_x = '0';
nodes1(m+8).position_y = "0";
nodes1(m+8).position_z = "500";
nodes1(m+8).scale_x = "0";
nodes1(m+8).color = "#000000";
nodes1(m+8).url = '';
nodes1(m+8).descr = '';
nodes1(m+8).unic_name = "8+Layer1";
nodes1(m+8).labelSize = '0';

end

% label quadranti
for iter = 1:1
nodes1(m+9).name = "Lack of Cognitive";
nodes1(m+9).labelPosition = "top";
nodes1(m+9).labelAdjustX = '0';
nodes1(m+9).labelAdjustY = '0';
nodes1(m+9).isBold = "true";
nodes1(m+9).layer = 'Layer1';
nodes1(m+9).position_x = '0';
nodes1(m+9).position_y = "-450";
nodes1(m+9).position_z = "-450";
nodes1(m+9).scale_x = "0";
nodes1(m+9).color = "#ffffff";
nodes1(m+9).url = '';
nodes1(m+9).descr = '';
nodes1(m+9).unic_name = "Lack of Cognitive+Layer1";
nodes1(m+9).labelSize = '20';

nodes1(m+10).name = "Lack_of Affective";
nodes1(m+10).labelPosition = "top";
nodes1(m+10).labelAdjustX = '0';
nodes1(m+10).labelAdjustY = '0';
nodes1(m+10).isBold = "true";
nodes1(m+10).layer = 'Layer1';
nodes1(m+10).position_x = '0';
nodes1(m+10).position_y = "-450";
nodes1(m+10).position_z = "180";
nodes1(m+10).scale_x = "0";
nodes1(m+10).color = "#ffffff";
nodes1(m+10).url = '';
nodes1(m+10).descr = '';
nodes1(m+10).unic_name = "Lack of Affective+Layer1";
nodes1(m+10).labelSize = '20';

nodes1(m+11).name = "Cognitive";
nodes1(m+11).labelPosition = "top";
nodes1(m+11).labelAdjustX = '0';
nodes1(m+11).labelAdjustY = '0';
nodes1(m+11).isBold = "true";
nodes1(m+11).layer = 'Layer1';
nodes1(m+11).position_x = '0';
nodes1(m+11).position_y = "450";
nodes1(m+11).position_z = "-450";
nodes1(m+11).scale_x = "0";
nodes1(m+11).color = "#ffffff";
nodes1(m+11).url = '';
nodes1(m+11).descr = '';
nodes1(m+11).unic_name = "Cognitive+Layer1";
nodes1(m+11).labelSize = '20';

nodes1(m+12).name = "Affective";
nodes1(m+12).labelPosition = "top";
nodes1(m+12).labelAdjustX = '0';
nodes1(m+12).labelAdjustY = '0';
nodes1(m+12).isBold = "true";
nodes1(m+12).layer = 'Layer1';
nodes1(m+12).position_x = '0';
nodes1(m+12).position_y = "450";
nodes1(m+12).position_z = "320";
nodes1(m+12).scale_x = "0";
nodes1(m+12).color = "#ffffff";
nodes1(m+12).url = '';
nodes1(m+12).descr = '';
nodes1(m+12).unic_name = "Affective+Layer1";
nodes1(m+12).labelSize = '20';

nodes1(m+13).name = "Well Being";
nodes1(m+13).labelPosition = "bottom";
nodes1(m+13).labelAdjustX = '0';
nodes1(m+13).labelAdjustY = '0';
nodes1(m+13).isBold = "true";
nodes1(m+13).layer = 'Layer1';
nodes1(m+13).position_x = '0';
nodes1(m+13).position_y = "-450";
nodes1(m+13).position_z = "-450";
nodes1(m+13).scale_x = "0";
nodes1(m+13).color = "#ffffff";
nodes1(m+13).url = '';
nodes1(m+13).descr = '';
nodes1(m+13).unic_name = "Well Being+Layer1";
nodes1(m+13).labelSize = '20';

nodes1(m+14).name = "Well Being";
nodes1(m+14).labelPosition = "bottom";
nodes1(m+14).labelAdjustX = '0';
nodes1(m+14).labelAdjustY = '0';
nodes1(m+14).isBold = "true";
nodes1(m+14).layer = 'Layer1';
nodes1(m+14).position_x = '0';
nodes1(m+14).position_y = "-450";
nodes1(m+14).position_z = "180";
nodes1(m+14).scale_x = "0";
nodes1(m+14).color = "#ffffff";
nodes1(m+14).url = '';
nodes1(m+14).descr = '';
nodes1(m+14).unic_name = "Well Being+Layer1";
nodes1(m+14).labelSize = '20';

nodes1(m+15).name = "Distress";
nodes1(m+15).labelPosition = "bottom";
nodes1(m+15).labelAdjustX = '0';
nodes1(m+15).labelAdjustY = '0';
nodes1(m+15).isBold = "true";
nodes1(m+15).layer = 'Layer1';
nodes1(m+15).position_x = '0';
nodes1(m+15).position_y = "450";
nodes1(m+15).position_z = "-450";
nodes1(m+15).scale_x = "0";
nodes1(m+15).color = "#ffffff";
nodes1(m+15).url = '';
nodes1(m+15).descr = '';
nodes1(m+15).unic_name = "Distress+Layer1";
nodes1(m+15).labelSize = '20';

nodes1(m+16).name = "Distress";
nodes1(m+16).labelPosition = "bottom";
nodes1(m+16).labelAdjustX = '0';
nodes1(m+16).labelAdjustY = '0';
nodes1(m+16).isBold = "true";
nodes1(m+16).layer = 'Layer1';
nodes1(m+16).position_x = '0';
nodes1(m+16).position_y = "450";
nodes1(m+16).position_z = "320";
nodes1(m+16).scale_x = "0";
nodes1(m+16).color = "#ffffff";
nodes1(m+16).url = '';
nodes1(m+16).descr = '';
nodes1(m+16).unic_name = "Distress+Layer1";
nodes1(m+16).labelSize = '20';
end



end

%intelaiatura L2
for iter = 1:1


m = 0;
m = max(size(nodes_unic_names_layer_2));

%vertex
for iter = 1:1

nodes2(m+1).name = "4";
nodes2(m+1).labelPosition = "Top";
nodes2(m+1).labelAdjustX = '0';
nodes2(m+1).labelAdjustY = '0';
nodes2(m+1).isBold = "true";
nodes2(m+1).layer = 'Layer2';
nodes2(m+1).position_x = '0';
nodes2(m+1).position_y = "500";
nodes2(m+1).position_z = "500";
nodes2(m+1).scale_x = "0";
nodes2(m+1).color = "#000000";
nodes2(m+1).url = '';
nodes2(m+1).descr = '';
nodes2(m+1).unic_name = "4+Layer2";
nodes2(m+1).labelSize = '0';

nodes2(m+2).name = "3";
nodes2(m+2).labelPosition = "Top";
nodes2(m+2).labelAdjustX = '0';
nodes2(m+2).labelAdjustY = '0';
nodes2(m+2).isBold = "true";
nodes2(m+2).layer = 'Layer2';
nodes2(m+2).position_x = '0';
nodes2(m+2).position_y = "500";
nodes2(m+2).position_z = "-500";
nodes2(m+2).scale_x = "0";
nodes2(m+2).color = "#000000";
nodes2(m+2).url = '';
nodes2(m+2).descr = '';
nodes2(m+2).unic_name = "3+Layer2";
nodes2(m+2).labelSize = '0';

nodes2(m+3).name = "1";
nodes2(m+3).labelPosition = "Bottom";
nodes2(m+3).labelAdjustX = '0';
nodes2(m+3).labelAdjustY = '0';
nodes2(m+3).isBold = "true";
nodes2(m+3).layer = 'Layer2';
nodes2(m+3).position_x = '0';
nodes2(m+3).position_y = "-500";
nodes2(m+3).position_z = "-500";
nodes2(m+3).scale_x = "0";
nodes2(m+3).color = "#000000";
nodes2(m+3).url = '';
nodes2(m+3).descr = '';
nodes2(m+3).unic_name = "1+Layer2";
nodes2(m+3).labelSize = '0';

nodes2(m+4).name = "2";
nodes2(m+4).labelPosition = "Bottom";
nodes2(m+4).labelAdjustX = '0';
nodes2(m+4).labelAdjustY = '0';
nodes2(m+4).isBold = "true";
nodes2(m+4).layer = 'Layer2';
nodes2(m+4).position_x = '0';
nodes2(m+4).position_y = "-500";
nodes2(m+4).position_z = "500";
nodes2(m+4).scale_x = "0";
nodes2(m+4).color = "#000000";
nodes2(m+4).url = '';
nodes2(m+4).descr = '';
nodes2(m+4).unic_name = "2+Layer2";
nodes2(m+4).labelSize = '0';

end

%quadranti
for iter = 1:1

nodes2(m+5).name = "5";
nodes2(m+5).labelPosition = "Bottom";
nodes2(m+5).labelAdjustX = '0';
nodes2(m+5).labelAdjustY = '0';
nodes2(m+5).isBold = "true";
nodes2(m+5).layer = 'Layer2';
nodes2(m+5).position_x = '0';
nodes2(m+5).position_y = "-500";
nodes2(m+5).position_z = "1";
nodes2(m+5).scale_x = "0";
nodes2(m+5).color = "#000000";
nodes2(m+5).url = '';
nodes2(m+5).descr = '';
nodes2(m+5).unic_name = "5+Layer2";
nodes2(m+5).labelSize = '0';

nodes2(m+6).name = "6";
nodes2(m+6).labelPosition = "Bottom";
nodes2(m+6).labelAdjustX = '0';
nodes2(m+6).labelAdjustY = '0';
nodes2(m+6).isBold = "true";
nodes2(m+6).layer = 'Layer2';
nodes2(m+6).position_x = '0';
nodes2(m+6).position_y = "500";
nodes2(m+6).position_z = "1";
nodes2(m+6).scale_x = "0";
nodes2(m+6).color = "#000000";
nodes2(m+6).url = '';
nodes2(m+6).descr = '';
nodes2(m+6).unic_name = "6+Layer2";
nodes2(m+6).labelSize = '0';

nodes2(m+7).name = "7";
nodes2(m+7).labelPosition = "Bottom";
nodes2(m+7).labelAdjustX = '0';
nodes2(m+7).labelAdjustY = '0';
nodes2(m+7).isBold = "true";
nodes2(m+7).layer = 'Layer2';
nodes2(m+7).position_x = '0';
nodes2(m+7).position_y = "0";
nodes2(m+7).position_z = "-500";
nodes2(m+7).scale_x = "0";
nodes2(m+7).color = "#000000";
nodes2(m+7).url = '';
nodes2(m+7).descr = '';
nodes2(m+7).unic_name = "7+Layer2";
nodes2(m+7).labelSize = '0';

nodes2(m+8).name = "8";
nodes2(m+8).labelPosition = "Bottom";
nodes2(m+8).labelAdjustX = '0';
nodes2(m+8).labelAdjustY = '0';
nodes2(m+8).isBold = "true";
nodes2(m+8).layer = 'Layer2';
nodes2(m+8).position_x = '0';
nodes2(m+8).position_y = "0";
nodes2(m+8).position_z = "500";
nodes2(m+8).scale_x = "0";
nodes2(m+8).color = "#000000";
nodes2(m+8).url = '';
nodes2(m+8).descr = '';
nodes2(m+8).unic_name = "8+Layer2";
nodes2(m+8).labelSize = '0';

end

% label quadranti
for iter = 1:1
nodes2(m+9).name = "Lack_of Cognitive";
nodes2(m+9).labelPosition = "top";
nodes2(m+9).labelAdjustX = '0';
nodes2(m+9).labelAdjustY = '0';
nodes2(m+9).isBold = "true";
nodes2(m+9).layer = 'Layer2';
nodes2(m+9).position_x = '0';
nodes2(m+9).position_y = "-450";
nodes2(m+9).position_z = "-450";
nodes2(m+9).scale_x = "0";
nodes2(m+9).color = "#ffffff";
nodes2(m+9).url = '';
nodes2(m+9).descr = '';
nodes2(m+9).unic_name = "Lack of Cognitive+Layer2";
nodes2(m+9).labelSize = '20';

nodes2(m+10).name = "Lack of Affective";
nodes2(m+10).labelPosition = "top";
nodes2(m+10).labelAdjustX = '0';
nodes2(m+10).labelAdjustY = '0';
nodes2(m+10).isBold = "true";
nodes2(m+10).layer = 'Layer2';
nodes2(m+10).position_x = '0';
nodes2(m+10).position_y = "-450";
nodes2(m+10).position_z = "180";
nodes2(m+10).scale_x = "0";
nodes2(m+10).color = "#ffffff";
nodes2(m+10).url = '';
nodes2(m+10).descr = '';
nodes2(m+10).unic_name = "Lack of Affective+Layer2";
nodes2(m+10).labelSize = '20';

nodes2(m+11).name = "Cognitive";
nodes2(m+11).labelPosition = "top";
nodes2(m+11).labelAdjustX = '0';
nodes2(m+11).labelAdjustY = '0';
nodes2(m+11).isBold = "true";
nodes2(m+11).layer = 'Layer2';
nodes2(m+11).position_x = '0';
nodes2(m+11).position_y = "450";
nodes2(m+11).position_z = "-450";
nodes2(m+11).scale_x = "0";
nodes2(m+11).color = "#ffffff";
nodes2(m+11).url = '';
nodes2(m+11).descr = '';
nodes2(m+11).unic_name = "Cognitive+Layer2";
nodes2(m+11).labelSize = '20';

nodes2(m+12).name = "Affective";
nodes2(m+12).labelPosition = "top";
nodes2(m+12).labelAdjustX = '0';
nodes2(m+12).labelAdjustY = '0';
nodes2(m+12).isBold = "true";
nodes2(m+12).layer = 'Layer2';
nodes2(m+12).position_x = '0';
nodes2(m+12).position_y = "450";
nodes2(m+12).position_z = "320";
nodes2(m+12).scale_x = "0";
nodes2(m+12).color = "#ffffff";
nodes2(m+12).url = '';
nodes2(m+12).descr = '';
nodes2(m+12).unic_name = "Affective+Layer2";
nodes2(m+12).labelSize = '20';

nodes2(m+13).name = "Well Being";
nodes2(m+13).labelPosition = "bottom";
nodes2(m+13).labelAdjustX = '0';
nodes2(m+13).labelAdjustY = '0';
nodes2(m+13).isBold = "true";
nodes2(m+13).layer = 'Layer2';
nodes2(m+13).position_x = '0';
nodes2(m+13).position_y = "-450";
nodes2(m+13).position_z = "-450";
nodes2(m+13).scale_x = "0";
nodes2(m+13).color = "#ffffff";
nodes2(m+13).url = '';
nodes2(m+13).descr = '';
nodes2(m+13).unic_name = "Well Being+Layer2";
nodes2(m+13).labelSize = '20';

nodes2(m+14).name = "Well Being";
nodes2(m+14).labelPosition = "bottom";
nodes2(m+14).labelAdjustX = '0';
nodes2(m+14).labelAdjustY = '0';
nodes2(m+14).isBold = "true";
nodes2(m+14).layer = 'Layer2';
nodes2(m+14).position_x = '0';
nodes2(m+14).position_y = "-450";
nodes2(m+14).position_z = "180";
nodes2(m+14).scale_x = "0";
nodes2(m+14).color = "#ffffff";
nodes2(m+14).url = '';
nodes2(m+14).descr = '';
nodes2(m+14).unic_name = "Well Being+Layer2";
nodes2(m+14).labelSize = '20';

nodes2(m+15).name = "Distress";
nodes2(m+15).labelPosition = "bottom";
nodes2(m+15).labelAdjustX = '0';
nodes2(m+15).labelAdjustY = '0';
nodes2(m+15).isBold = "true";
nodes2(m+15).layer = 'Layer2';
nodes2(m+15).position_x = '0';
nodes2(m+15).position_y = "450";
nodes2(m+15).position_z = "-450";
nodes2(m+15).scale_x = "0";
nodes2(m+15).color = "#ffffff";
nodes2(m+15).url = '';
nodes2(m+15).descr = '';
nodes2(m+15).unic_name = "Distress+Layer2";
nodes2(m+15).labelSize = '20';

nodes2(m+16).name = "Distress";
nodes2(m+16).labelPosition = "bottom";
nodes2(m+16).labelAdjustX = '0';
nodes2(m+16).labelAdjustY = '0';
nodes2(m+16).isBold = "true";
nodes2(m+16).layer = 'Layer2';
nodes2(m+16).position_x = '0';
nodes2(m+16).position_y = "450";
nodes2(m+16).position_z = "320";
nodes2(m+16).scale_x = "0";
nodes2(m+16).color = "#ffffff";
nodes2(m+16).url = '';
nodes2(m+16).descr = '';
nodes2(m+16).unic_name = "Distress  +Layer2";
nodes2(m+16).labelSize = '20';
end


end




%value4 = [nodes1,nodes2];

end


%% field 5 edges
for iter = 1:1
    
field5 = 'edges';




% edges 1 %%%%%%%%%
edges_1 = [];
for iter = 1:max(size(src_names_edges_cross_1))
    edges_1(iter).src = src_names_edges_cross_1(iter);
    edges_1(iter).trg = trg_names_edges_cross_1(iter);
    edges_1(iter).opacity = '0.3';
    edges_1(iter).color = COLORS1(iter);
    edges_1(iter).channel = '';
    %edges_1(iter).size = WEIGHT_ARENA1(iter)/5; %string(cc1(iter)*4);
    edges_1(iter).size = 1; %log(WEIGHT_ARENA1(iter)*2);
    edges_1(iter).arrow = "true";
end


m = 0;
if exist('edges_1','var')
    m = numel(edges_1);
end

edges_1(m+1).src = "1+Layer1";
edges_1(m+1).trg = "2+Layer1";
edges_1(m+1).opacity = '0.3';
edges_1(m+1).color = color_quadrants;  % #b9bdc4 grigio
edges_1(m+1).channel = '';
edges_1(m+1).size = "2";
edges_1(m+1).arrow = "false";

edges_1(m+2).src = "1+Layer1";
edges_1(m+2).trg = "3+Layer1";
edges_1(m+2).opacity = '0.3';
edges_1(m+2).color = color_quadrants;
edges_1(m+2).channel = '';
edges_1(m+2).size = "2";
edges_1(m+2).arrow = "false";

edges_1(m+3).src = "2+Layer1";
edges_1(m+3).trg = "4+Layer1";
edges_1(m+3).opacity = '0.3';
edges_1(m+3).color = color_quadrants;
edges_1(m+3).channel = '';
edges_1(m+3).size = "2";
edges_1(m+3).arrow = "false";

edges_1(m+4).src = "3+Layer1";
edges_1(m+4).trg = "4+Layer1";
edges_1(m+4).opacity = '0.3';
edges_1(m+4).color = color_quadrants;
edges_1(m+4).channel = '';
edges_1(m+4).size = "2";
edges_1(m+4).arrow = "false";

edges_1(m+5).src = "6+Layer1";
edges_1(m+5).trg = "5+Layer1";
edges_1(m+5).opacity = '0.3';
edges_1(m+5).color = "#000000";
edges_1(m+5).channel = '';
edges_1(m+5).size = "1";
edges_1(m+5).arrow = "false";

edges_1(m+6).src = "7+Layer1";
edges_1(m+6).trg = "8+Layer1";
edges_1(m+6).opacity = '0.3';
edges_1(m+6).color = "#000000";
edges_1(m+6).channel = '';
edges_1(m+6).size = "1";
edges_1(m+6).arrow = "false";


% edges 1_2 %%%%%%%%%
edges_1_2 = [];
for iter = 1:max(size(src_names_edges_longi_1_2))
    edges_1_2(iter).src = src_names_edges_longi_1_2(iter);
    edges_1_2(iter).trg = trg_names_edges_longi_1_2(iter);
    edges_1_2(iter).opacity = '0.3';
    edges_1_2(iter).color = COLORS1_2(iter);
    edges_1_2(iter).channel = '';
    %edges_1_2(iter).size = WEIGHT_ARENA1_2(iter)/5; %string(cc1_2(iter)*4);
    edges_1_2(iter).size = 1; % log(WEIGHT_ARENA1_2(iter)*2); 
    edges_1_2(iter).arrow = "true";
end



% edges 2 %%%%%%%%%
edges_2 = [];
for iter = 1:max(size(src_names_edges_cross_2))
    edges_2(iter).src = src_names_edges_cross_2(iter);
    edges_2(iter).trg = trg_names_edges_cross_2(iter);
    edges_2(iter).opacity = '0.3';
    edges_2(iter).color = COLORS2(iter);
    edges_2(iter).channel = '';
    %edges_2(iter).size = WEIGHT_ARENA2(iter)/5; %string(cc2(iter)*4);
    edges_2(iter).size = 1; % log(WEIGHT_ARENA2(iter)*2);
    edges_2(iter).arrow = "true";
end


m = 0;
if exist('edges_2','var')
    m = numel(edges_2);
end

edges_2(m+1).src = "1+Layer2";
edges_2(m+1).trg = "2+Layer2";
edges_2(m+1).opacity = '0.3';
edges_2(m+1).color = color_quadrants;
edges_2(m+1).channel = '';
edges_2(m+1).size = "2";
edges_2(m+1).arrow = "false";

edges_2(m+2).src = "1+Layer2";
edges_2(m+2).trg = "3+Layer2";
edges_2(m+2).opacity = '0.3';
edges_2(m+2).color = color_quadrants;
edges_2(m+2).channel = '';
edges_2(m+2).size = "2";
edges_2(m+2).arrow = "false";

edges_2(m+3).src = "2+Layer2";
edges_2(m+3).trg = "4+Layer2";
edges_2(m+3).opacity = '0.3';
edges_2(m+3).color = color_quadrants;
edges_2(m+3).channel = '';
edges_2(m+3).size = "2";
edges_2(m+3).arrow = "false";

edges_2(m+4).src = "3+Layer2";
edges_2(m+4).trg = "4+Layer2";
edges_2(m+4).opacity = '0.3';
edges_2(m+4).color = color_quadrants;
edges_2(m+4).channel = '';
edges_2(m+4).size = "2";
edges_2(m+4).arrow = "false";

edges_2(m+5).src = "5+Layer2";
edges_2(m+5).trg = "6+Layer2";
edges_2(m+5).opacity = '0.3';
edges_2(m+5).color = "#000000";
edges_2(m+5).channel = '';
edges_2(m+5).size = "1";
edges_2(m+5).arrow = "false";

edges_2(m+6).src = "7+Layer2";
edges_2(m+6).trg = "8+Layer2";
edges_2(m+6).opacity = '0.3';
edges_2(m+6).color = "#000000";
edges_2(m+6).channel = '';
edges_2(m+6).size = "1";
edges_2(m+6).arrow = "false";

    
%value5 = [edges_1,edges_1_2,edges_2];


end


value4 = [nodes1,nodes2];
value5 = [edges_1,edges_1_2,edges_2];



%% def struct
clear s;
s = struct(field3,value3,field4,value4,field5,value5);
s2 = s;


end


%% Create final json that combine the s1 and s2
% This section creates the final structure that will then be exported to a json file. This structure is the union of the previous two
for iter = 1:1

clear s
s = s1;

c = numel(s1.edges) + 1;
for iter1 = 1:numel(s2.edges)
    a = 0;
    for iter2 = 1:numel(s1.edges)
        if strcmp(s1.edges(iter2).src, s2.edges(iter1).src) & strcmp(s1.edges(iter2).trg, s2.edges(iter1).trg)
            a = 1;
            break
        end
    end
    if a == 0
        s.edges(c) = s2.edges(iter1)
        c = c + 1;
    end
end

end


%% create jsonfile
j = jsonencode(s);
file = prettyjson(j)


%% save json
for iter = 1:1

% this section saves the json files at each interaction of the loop

% filename = strcat("json", "_", "22q", "_", string(iter_cluster), ".json");
% 
% %save file json
% fid = fopen(filename,'w');
% fprintf(fid,'%s',file);
% fclose(fid);

end


%% open mlnetwork
% this section causes a browser (we recommend chrome) to automatically open to the MLNetwork page with the figure displayed.
for iter = 1:1

httpsUrl = 'https://dev.mlnetwork-diplab.ch/cache.php';
options = weboptions('RequestMethod','POST');
responseEmployee = webwrite(httpsUrl, file, options);
web(responseEmployee);

end


