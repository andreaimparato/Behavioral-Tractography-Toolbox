clc
clear all
close all


%% load EMA data and names
for iter = 1:1

% load data set
load EMA_data.mat

% selecting only controls with IQ evaluated
try_table(not(strcmp(try_table.Diagnosis,'Normal_Control')),:) = [];
try_table(isnan(try_table.FSIQ),:) = [];
try_table((try_table.FSIQ==0),:) = [];


%Load symptoms / nodes names
load names_3_v4    

end


%% Create Variables 
for iter = 1:1


Table_names=try_table.Properties.VariableNames';

Symptom_Identifier=[12:31];
Symptoms=try_table{:,Symptom_Identifier};
Names_Symptoms=Table_names(Symptom_Identifier)';
Names_Symptoms=Names_Symptoms';
Notification_No=try_table{:,'NotificationNo'};
ParticipantID=try_table{:,'ParticipantID'};
NotificationNo=try_table{:,'NotificationNo'};
Groupe=try_table{:,'Diagnosis'};
Session_Istance=try_table{:,'SessionInstance'};
Responded=try_table{:,'Responded'};
Completed=try_table{:,'CompletedSession'};


%% Filters or completed time-points
Symptoms=Symptoms(Completed==1,:);
Notification_No=Notification_No(Completed==1,:);
ParticipantID=ParticipantID(Completed==1,:);
NotificationNo=NotificationNo(Completed==1,:);
Groupe=Groupe(Completed==1,:);
Session_Istance=Session_Istance(Completed==1,:);
Responded=Responded(Completed==1,:);
Completed=Completed(Completed==1,:);

% Reverse code positive symptoms
Symptoms_to_reverse_code=[1,4,7,9,15,16,18,19];
Reverse_code=Names_Symptoms(Symptoms_to_reverse_code);
Symptoms(:,Symptoms_to_reverse_code)=-Symptoms(:,Symptoms_to_reverse_code);

% Modify names of symptoms that are reverse coded
for iter=1:max(size(Symptoms_to_reverse_code))
    b=(Names_Symptoms(Symptoms_to_reverse_code(iter),1));
    Names_Symptoms(Symptoms_to_reverse_code(iter),1)=strcat('NOT_',b);    
end


end


%% valid tp for crosssectional
for iter = 1:1
% keep only timepoints that have a consecutive or previous timepoint

if abs(NotificationNo(1) - NotificationNo(2)) == 1 
    vtp(1) = 1;
else
    vtp(1) = 0;
end

for iter = 2:(max(size(NotificationNo))-1)
    if ( abs(NotificationNo(iter) - NotificationNo(iter-1)) == 1 | abs(NotificationNo(iter) - NotificationNo(iter+1)) == 1) == 1
        vtp(iter) = 1;
    else
        vtp(iter) = 0;
    end
end
        
if abs(NotificationNo(max(size(NotificationNo))) - NotificationNo(max(size(NotificationNo))-1)) == 1 
    vtp(max(size(NotificationNo))) = 1;
else
    vtp(max(size(NotificationNo))) = 0;
end    

vtp = vtp';

Symptoms = Symptoms(find(vtp),:);
ParticipantID = ParticipantID(find(vtp));
Groupe = Groupe(find(vtp));
Notification_No = Notification_No(find(vtp));
Responded = Responded(find(vtp));
Completed = Completed(find(vtp));


end


%% chose groupe for matrix
for iter = 1:1

Chose_Group = "Normal_Control";
max(size(unique(ParticipantID(strcmp(Groupe,Chose_Group)))))

end


%% compute cross matrix
for iter = 1:1

F1 = find(strcmp(Groupe,Chose_Group));
[R2,P2] = cross_corr_c(Symptoms,F1,ParticipantID);


% Make adjecency matrix perfectly simmetrical
VECTOR_R=jUpperTriMatToVec(R2,1);
WHICH_NETWORK=Contraire_jUpperT(VECTOR_R,max(size(R2)),1);

VECTOR_P=jUpperTriMatToVec(P2,1);
WHICH_P_VALUE=Contraire_jUpperT(VECTOR_P,max(size(P2)),1);

% Threshold Network
[~,P_Tresh]=fdr_bh(VECTOR_P);
WHICH_NETWORK(WHICH_P_VALUE>P_Tresh)=0;
WHICH_NETWORK(WHICH_NETWORK<0)=0;

network_cross = WHICH_NETWORK;

end


%% compute longi matrix
for iter = 1:1


% Identify timepoints for time-lagged correlations

UniqueID=unique(ParticipantID);

Potential_T1=[1:7,9:15,17:23,25:31,33:39,41:47];
Potential_T2=Potential_T1+1;

FIND_T1=[];
FIND_T2=[];

for iter= 1:max(size(UniqueID))
    for iter2=1:max(size(Potential_T1))      
       if Responded(ParticipantID==UniqueID(iter)&Notification_No==Potential_T1(iter2),1)==1 & Responded(ParticipantID==UniqueID(iter)&Notification_No==Potential_T2(iter2))==1        
           find_t1=find(ParticipantID==UniqueID(iter)&Notification_No==Potential_T1(iter2));
           find_t2=find(ParticipantID==UniqueID(iter)&Notification_No==(Potential_T1(iter2)+1));
           FIND_T1=[FIND_T1;find_t1];
           FIND_T2=[FIND_T2;find_t2];
       end  
    end    
end

FIND_T1_DEF=FIND_T1(FIND_T1>0&FIND_T2>0,:);
FIND_T2_DEF=FIND_T2(FIND_T1>0&FIND_T2>0,:);

% Extract longitudinal time-points for time lagged correlations
Id_longi_1=ParticipantID(FIND_T1_DEF);
Id_longi_2=ParticipantID(FIND_T2_DEF);
Session_Istance_T1=Session_Istance(FIND_T1_DEF);
Session_Istance_T2=Session_Istance(FIND_T2_DEF);
Group_Longi=Groupe(FIND_T1_DEF);
Symptoms_T1=Symptoms(FIND_T1_DEF,:);
Symptoms_T2=Symptoms(FIND_T2_DEF,:);

% Compute multilayer correlation matrix
Longitudinal_data=[Symptoms_T1,Symptoms_T2];
Chose_group=Chose_Group;


[ML_NETWORK,ML_PVAL]=MIXED_CORR_MATRIX_c(Longitudinal_data(strcmp(Group_Longi,Chose_Group),:),Id_longi_1(strcmp(Group_Longi,Chose_Group),:));


% Make Matrix Perfectly Simmetrical
VECTOR_R=jUpperTriMatToVec(ML_NETWORK,1);
WHICH_NETWORK=Contraire_jUpperT(VECTOR_R,max(size(ML_NETWORK)),1);


% Treshold network 
P_Vector=jUpperTriMatToVec(ML_PVAL,1);
[~,P_Tresh]=fdr_bh(P_Vector);
WHICH_NETWORK(ML_PVAL>P_Tresh)=0;
% WHICH_NETWORK(WHICH_NETWORK<0)=0

network = WHICH_NETWORK;


end


%% set variables
for iter = 1:1

    % set variables
    dim_cross = max(size(network))/2;  % Calculate the size of 'network_cross'
    dim_longi = max(size(network));     % Calculate the size of 'network'

    SOURCE = 1:dim_cross;                 % Define SOURCE as indices from 1 to dim_cross
    TARGET = (dim_cross + 1) : dim_longi; % Define TARGET as indices from dim_cross+1 to dim_longi

    % names
    names_3_t1 = string(names_3(1:dim_cross));  % Extract names for the first half
    names_3_t2 = string(names_3(1:dim_cross));  % Extract names for the second half

    % color quadrants
    color_quadrants = "#000000";  % Define color for quadrants (currently set to black)

end


%% creating final adj matrix
for iter = 1:1


% remove all negative edges
network_cross(network_cross<0) = 0;   % Set all negative values in 'network_cross' to 0
network(network<0) = 0;               % Set all negative values in 'network' to 0

% use the cross network for layer_1 and layer_2
network(1:dim_cross, 1:dim_cross) = network_cross;                        % Assign 'network_cross' to the upper left quadrant of 'network'
network(dim_cross + 1:dim_longi, dim_cross + 1:dim_longi) = network_cross; % Assign 'network_cross' to the lower right quadrant of 'network'


end


%% choose symptom (use to plot only one sympt)
for iter = 1:1

% sympt = 18;  % Define the value of sympt
% 
% for iter2 = 1:40   % Loop over values from 1 to 40
%     for iter3 = 1:40   % Nested loop over values from 1 to 40
%         if (iter2 == sympt | iter3 == sympt) == 1   % Check if iter2 or iter3 is equal to sympt
%             network(iter2,iter3) = network(iter2,iter3);  % If true, keep the original value of network(iter2,iter3)
%         else
%             network(iter2,iter3) = 0;  % If false, set network(iter2,iter3) to 0
%         end
%     end
% end

end


%% Shortest paths
for iter = 1:1
% This section computes the shortest paths in the network and counts the number of shortest paths that pass through each node



% remove the diagonal
for iter = 1:max(size(network))
    network(iter,iter)= 0;  % Set diagonal elements of 'network' to 0
end

% remove connections going back in time
network(TARGET,SOURCE)=0;  % Set connections from TARGET to SOURCE to 0

% compute shortest paths (SP)
D=1./posweights(network);  % Compute the distance matrix from inverse of positive weights
[SPL,hops,Pmat] = distance_wei_floyd(D);  % Compute shortest paths using Floyd's algorithm

ACROSS_BETW=zeros(max(size(network)),max(size(network)));  % Initialize matrix for counting shortest paths between nodes
BTW = zeros(max(size(network)),1);  % Initialize vector for counting shortest paths passing through each node

for iter1=1:max(size(SOURCE))  % Loop over SOURCE nodes
    for iter2=1:max(size(TARGET))  % Loop over TARGET nodes
        B=retrieve_shortest_path(SOURCE(iter1),TARGET(iter2),hops,Pmat);  % Retrieve the shortest path between SOURCE and TARGET nodes
        % count shortest paths (sp)
        ACROSS_BETW2=zeros(max(size(network)),max(size(network)));  % Initialize matrix for counting shortest paths between nodes along a specific path
        for iter3=1:(max(size(B))-1)  % Loop over elements of the shortest path B
            ACROSS_BETW2(B(iter3),B(iter3+1))=1;  % Mark edges along the shortest path in ACROSS_BETW2
        end
        ACROSS_BETW=ACROSS_BETW+ACROSS_BETW2;  % Accumulate shortest path counts in ACROSS_BETW
        % count betweenness (btw)
        BTW2=zeros(max(size(network)),1);  % Initialize vector for counting shortest paths passing through nodes in the shortest path
        BTW2(B)=1;  % Mark nodes in the shortest path
        BTW=BTW+BTW2;  % Accumulate shortest path counts passing through each node
    end
end

net_sp = ACROSS_BETW;  % Assign the matrix of shortest path counts to net_sp



end


%% Compute BTW
for iter = 1:1


% This section perform a randomization test to assess the significance of betweenness centrality values in the network.

% Clear existing variable P_T1_BETW_2TAIL
clear P_T1_BETW_2TAIL

% Compute betweenness centrality for edges from SOURCE to TARGET
ACROSS_BETW=BTW_SOURCE_TARGET_DIR(network,SOURCE,TARGET);
% ACROSS_BETW=BTW_SOURCE_TARGET(WHICH_NETWORK,SOURCE,TARGET);

n_rand=1000;  % Define the number of randomizations

clear VECT  % Clear VECT variable
clear VECTI  % Clear VECTI variable
clear ACROSS_BETW_RAND  % Clear ACROSS_BETW_RAND variable

for iter=1:n_rand  % Loop for the number of randomizations
    iter;  % Display the current iteration number
    
    % Generate a random permutation of upper triangular part of the network
    VECT=jUpperTriMatToVec(network,1);
    VECT=VECT(randperm(max(size(VECT))));
    VECTI(:,iter)=VECT;  % Store the random permutation in VECTI
    
    % Construct matrix from the random permutation
    MATRIX=Contraire_jUpperT(VECT,max(size(network)),1);
    
    % Compute betweenness centrality for edges from SOURCE to TARGET using the random permutation
    ACROSS_BETW_RAND(:,iter)=BTW_SOURCE_TARGET_DIR(MATRIX,SOURCE,TARGET);
    % ACROSS_BETW_RAND(:,iter)=BTW_SOURCE_TARGET(MATRIX,SOURCE,TARGET);
end

% Compute p-values for betweenness centrality
for iter=1:max(size(ACROSS_BETW))
    % Calculate the proportion of random betweenness centralities greater than the observed value
    P_T1_BETW_2TAIL(iter,1)=max(size(find(abs(ACROSS_BETW_RAND(iter,:))>abs(ACROSS_BETW(iter)))))/(n_rand+1);
end

% Perform false discovery rate (FDR) correction
[~,P_Tresh]=fdr_bh(P_T1_BETW_2TAIL);


end


%% bold labels
for iter = 1:1
% This section calculates boolean vectors indicating significant correlations
% after FDR correction and correlations with p-value < 0.05 without FDR correction,
% and then converts them into string representations.

% Initialize vectors to store bold results
bold_corr = zeros(dim_longi,1);   % Vector to indicate significant BTW after FDR correction
bold_no_corr = zeros(dim_longi,1); % Vector to indicate BTW with p-value < 0.05 without FDR correction
p_value_btw = P_T1_BETW_2TAIL;     % Store the p-values for betweenness centrality

% Identify significant BTW after FDR correction
bold_corr(P_T1_BETW_2TAIL <= P_Tresh) = 1;

% Identify BTW with p-value < 0.05 without FDR correction
bold_no_corr(P_T1_BETW_2TAIL < 0.05) = 1;

% Initialize string vectors to store boolean values
bold_corr_str = strings(dim_longi,1);   % String vector to indicate true/false for significant BTW after FDR correction
bold_no_corr_str = strings(dim_longi,1); % String vector to indicate true/false for BTW with p-value < 0.05 without FDR correction

% Convert boolean vectors to strings
for iter = 1:max(size(bold_corr))
    if bold_corr(iter) == 1
        bold_corr_str(iter) = "true";  % Assign "true" if the BTW is significant after FDR correction
    else
        bold_corr_str(iter) = "false"; % Assign "false" otherwise
    end
end

for iter = 1:max(size(bold_no_corr))
    if bold_no_corr(iter) == 1
        bold_no_corr_str(iter) = "true";  % Assign "true" if the BTW has p-value < 0.05 without FDR correction
    else
        bold_no_corr_str(iter) = "false"; % Assign "false" otherwise
    end
end


end


%% label size
for iter = 1:1
% This section assigns label sizes based on whether the btw is significant (bold_no_corr(iter) == 1)


% Initialize string vector to store label sizes
label_size = strings(dim_longi,1);

% Loop through each element in bold_no_corr
for iter = 1:max(size(bold_no_corr))
    % Check if the btw is significant
    if bold_no_corr(iter) == 1
        % Check if the index is beyond the cross network
        if iter > dim_cross
            label_size(iter) = "20";  % If beyond the cross network, assign label size of 20
        else
            label_size(iter) = "25";  % If within the cross network, assign label size of 25
        end
    else
        % If correlation is not significant
        % Check if the index is beyond the cross network
        if iter > dim_cross
            label_size(iter) = "12";  % If beyond the cross network, assign label size of 12
        else
            label_size(iter) = "15";  % If within the cross network, assign label size of 15
        end
    end
end



end


%% input data for the export
for iter=1:1

% This code generates the data for a visual representation of a network
% It involves processing node and edge data, adjusting coordinates,
% computing colors and sizes, and setting up arrows for visualization.

%% nodes
% PCA
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(network(1:dim_cross,1:dim_cross));  % Perform PCA on the cross network
XY=[COEFF(:,1),COEFF(:,2)];  % Extract XY coordinates from the principal components
XY_CROSS = XY;  % Store the XY coordinates for the cross network

% Adjusting coordinates
max_x = max(XY_CROSS(:,1));
max_y = max(XY_CROSS(:,2));
min_x = min(XY_CROSS(:,1));
min_y = min(XY_CROSS(:,2));
centre_x = (max_x + min_x)/2;  % Calculate the center x-coordinate
centre_y = (max_y + min_y)/2;  % Calculate the center y-coordinate

% Adjust XY_CROSS coordinates
aaaa = XY_CROSS;
XY_CROSS(:,1) = XY_CROSS(:,1) - centre_x;
XY_CROSS(:,2) = XY_CROSS(:,2) - centre_y;
XY_CROSS(:,1) = XY_CROSS(:,1)*1.2;  % Adjust x-coordinates
XY_CROSS(:,2) = XY_CROSS(:,2)*1.2;  % Adjust y-coordinates

% Label position top down
for iter = 1:max(size(XY_CROSS(:,2)))
    if XY_CROSS(iter,2) < 0
        labelPosition(iter) = "bottom";  % Set label position as 'bottom' if y-coordinate is negative
    else
        labelPosition(iter) = "top";  % Set label position as 'top' if y-coordinate is non-negative
    end
end

% Label length and position
f = 130/11;

for iter = 1:dim_cross
    label_lenght(iter) = strlength(names_3(iter));  % Compute label length for each node name
end

for iter = 1:max(size(XY_CROSS(:,1)))
    if XY_CROSS(iter,1) < 0
        labelAdjustX(iter) = -label_lenght(iter)*f;  % Adjust label X-coordinate for nodes with negative X-coordinate
    else
        labelAdjustX(iter) = "0";  % Set label X-coordinate as 0 for other nodes
    end
end

labelAdjustX(11) = "0"; % Special adjustment for certain nodes
labelAdjustX(12) = "0"; % Special adjustment for certain nodes
labelAdjustX(4) = -(label_lenght(4)*f+100); % Special adjustment for certain nodes
labelAdjustX = [labelAdjustX,labelAdjustX];
labelAdjustX(24) = 0; % Special adjustment for certain nodes

% Adjusting coordinates
coordinates_layer_1 = string(XY_CROSS*1000);  % Scale and convert XY_CROSS coordinates to string
coordinates_layer_2 = string(XY_CROSS*1000);  % Scale and convert XY_CROSS coordinates to string

% Node size
nodes_size_layer_1 = string(sum(network(1:20,1:20)));  % Compute node size for layer 1
nodes_size_layer_2 = string(sum(network(21:40,21:40)));  % Compute node size for layer 2

% Hub scaling based on p-value
max_btw = max(abs(log(p_value_btw)));  % Compute maximum p-value
for iter = 1:max(size(bold_no_corr))
    rgb_value = round(255 - (255*(abs(log(p_value_btw(iter)))/max_btw)));  % Compute RGB value based on p-value
    rgb_y = round(100 - (100*(abs(log(p_value_btw(iter)))/max_btw)));  % Compute RGB value based on p-value
    rgb_g = round(100 - (100*(abs(log(p_value_btw(iter)))/max_btw)));  % Compute RGB value based on p-value
    if iter <= 20
        if bold_no_corr(iter) == 1
            COLORS(iter,:) = rgb2hex([255 200 0]); % Set color to yellow for significant correlations
        else
            COLORS(iter,:) = rgb2hex([(255-rgb_value) (200-rgb_y) 0]); % Set color from green to yellow based on p-value
        end
    else
        if bold_no_corr(iter) == 1
            COLORS(iter,:) = rgb2hex([0 100 0]); % Set color to green for significant correlations
        else
            COLORS(iter,:) = rgb2hex([rgb_value (100+rgb_g) 0]); % Set color from yellow to green based on p-value
        end
    end
end

COLORS = string(COLORS);  % Convert COLORS to string

nodes_color_layer_1 = (COLORS(1:dim_cross));  % Set node colors for layer 1
nodes_color_layer_2 = (COLORS(dim_cross+1:dim_longi));  % Set node colors for layer 2

% Node names
nodes_names_layer_1 = string(names_3(1:dim_cross));  % Set node names for layer 1
nodes_names_layer_2 = string(names_3(1:dim_cross));  % Set node names for layer 2

% Nodes unique names
count = 1;
for iter=1:20
    nodes_unic_names_layer_1(count) = strcat(names_3(iter),'+Layer1');  % Append layer information to node names for layer 1
    nodes_unic_names_layer_2(count) = strcat(names_3(iter),'+Layer2');  % Append layer information to node names for layer 2
    count = count + 1;
end
nodes_unic_names_layer_1 = string(nodes_unic_names_layer_1)';  % Convert to string
nodes_unic_names_layer_2 = string(nodes_unic_names_layer_2)';  % Convert to string

%% edges
% Names
count = 1;
for iter=1:dim_cross
    for iter2=1:dim_cross
        src_names_edges_cross_1(count) = strcat(names_3(iter),'+Layer1');  % Set source names for cross network layer 1
        trg_names_edges_cross_1(count) = strcat(names_3(iter2),'+Layer1');  % Set target names for cross network layer 1
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
        src_names_edges_longi_1_2(count) = strcat(names_3(iter),'+Layer1');  % Set source names for longitudinal network layer 1
        trg_names_edges_longi_1_2(count) = strcat(names_3(iter2-dim_cross),'+Layer2');  % Set target names for longitudinal network layer 2
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
        src_names_edges_cross_2(count) = strcat(names_3(iter),'+Layer2');  % Set source names for cross network layer 2
        trg_names_edges_cross_2(count) = strcat(names_3(iter2),'+Layer2');  % Set target names for cross network layer 2
        count = count + 1;
    end
end
src_names_edges_cross_2 = src_names_edges_cross_2'; 
src_names_edges_cross_2 = string(src_names_edges_cross_2);
trg_names_edges_cross_2 = trg_names_edges_cross_2';
trg_names_edges_cross_2 = string(trg_names_edges_cross_2);

% Weights
W = network;  % Set weights for the network

count=1;
for iter=1:dim_cross
    for iter2=1:dim_cross
        WEIGHT_ARENA1(count) = W(iter,iter2);  % Set weights for cross network layer 1
        sp1(count) = net_sp(iter,iter2);  % Set shortest path info for cross network layer 1
        count = count + 1;
    end
end
WEIGHT_ARENA1 = WEIGHT_ARENA1';

count=1;
for iter=1:dim_cross
    for iter2=dim_cross+1:dim_longi
        WEIGHT_ARENA1_2(count) = W(iter,iter2);  % Set weights for longitudinal network layer 1 to layer 2
        sp1_2(count) = net_sp(iter,iter2);  % Set shortest path info for longitudinal network layer 1 to layer 2
        count = count + 1;
    end
end
WEIGHT_ARENA1_2 = WEIGHT_ARENA1_2';

count=1;
for iter=dim_cross+1:dim_longi
    for iter2=dim_cross+1:dim_longi
        WEIGHT_ARENA2(count) = W(iter,iter2);  % Set weights for cross network layer 2
        sp2(count) = net_sp(iter,iter2);  % Set shortest path info for cross network layer 2
        count = count + 1;
    end
end
WEIGHT_ARENA2 = WEIGHT_ARENA2';

% Size edge
sp_size_l1 = sp1;  
sp_size_l1_2 = sp1_2;  
sp_size_l2 = sp2;  

m1 = min(sp1(sp1>0)); 
m1_2 = min(sp1_2(sp1_2>0)); 
m2 = min(sp2(sp2>0));

sp_size_l1(sp_size_l1==0) = 1; 
sp_size_l1_2(sp_size_l1_2==0) = 1;
sp_size_l2(sp_size_l2==0) = 1;

% Colors edges
max_weight_long = max(log(sp_size_l1_2 + 1));  % Compute maximum weight for longitudinal network layer 1 to layer 2

a = max(log(sp_size_l1 + 1));
b = max(log(sp_size_l2 + 1));
d = [a;b];
max_weight_cross = max(d);  % Compute maximum weight for cross networks

max_weight = max(max_weight_cross,max_weight_long);  % Compute maximum weight

clear a b c d;  % Clear variables

% Set edge colors and arrows for longitudinal network layer 1 to layer 2
for iter = 1:max(size(WEIGHT_ARENA1_2))
    rgb_value = round(255 - (255*(log(sp_size_l1_2(iter)+1)/max_weight_long)));  % Compute RGB value based on weight
    rgb_value2 = round(180 - (180*(log(sp_size_l1_2(iter)+1)/max_weight_long)));  % Compute RGB value based on weight
    rgb_value = string(rgb_value);  % Convert to string
    rgb_value2 = string(rgb_value2);  % Convert to string
    if sp1_2(iter) > 0
        COLORS1_2(iter) = strcat("rgba","(",rgb_value2,",",rgb_value2,",",rgb_value2,",","1",")");
    else
        COLORS1_2(iter) = strcat("rgba","(",rgb_value,",",rgb_value,",",rgb_value,",","1",")");
    end
    COLORS1_2=string(COLORS1_2);  % Convert to string
end

% Set edge colors and arrows for cross network layer 1
for iter = 1:max(size(WEIGHT_ARENA1))
    rgb_value = round(255 - (255*(log(sp_size_l1(iter)+1)/max_weight_cross)));  % Compute RGB value based on weight
    rgb_value2 = round(180 - (180*(log(sp_size_l1(iter)+1)/max_weight_cross)));  % Compute RGB value based on weight
    rgb_value = string(rgb_value);  % Convert to string
    rgb_value2 = string(rgb_value2);  % Convert to string
    if sp1(iter) > 0
        COLORS1(iter) = strcat("rgba","(",rgb_value2,",",rgb_value2,",",rgb_value2,",","1",")");
    else
        COLORS1(iter) = strcat("rgba","(",rgb_value,",",rgb_value,",",rgb_value,",","1",")");
    end
    COLORS1=string(COLORS1);  % Convert to string
end

% Set edge colors and arrows for cross network layer 2
for iter = 1:max(size(WEIGHT_ARENA2))
    rgb_value = round(255 - (255*(log(sp_size_l2(iter)+1)/max_weight_cross)));  % Compute RGB value based on weight
    rgb_value2 = round(180 - (180*(log(sp_size_l2(iter)+1)/max_weight_cross)));  % Compute RGB value based on weight
    rgb_value = string(rgb_value);  % Convert to string
    rgb_value2 = string(rgb_value2);  % Convert to string
    if sp2(iter) > 0
        COLORS2(iter) = strcat("rgba","(",rgb_value2,",",rgb_value2,",",rgb_value2,",","1",")");
    else
        COLORS2(iter) = strcat("rgba","(",rgb_value,",",rgb_value,",",rgb_value,",","1",")");
    end
    COLORS2=string(COLORS2);  % Convert to string
end

% Set arrows for cross network layer 1
for iter = 1:max(size(sp1))
    if sp1(iter) > 0
        arrows_l1(iter) = "true";  % Set arrows as 'true' for shortest paths
    else
        arrows_l1(iter) = "false";  % Set arrows as 'false' for non shortest paths
    end
end

% Set arrows for longitudinal network layer 1 to layer 2
for iter = 1:max(size(sp1_2))
    if sp1_2(iter) > 0
        arrows_l1_2(iter) = "true";  % Set arrows as 'true' for shortest paths
    else
        arrows_l1_2(iter) = "false";  % Set arrows as 'false' for non shortest paths
    end
end

% Set arrows for cross network layer 2
for iter = 1:max(size(sp2))
    if sp1_2(iter) > 0
        arrows_l2(iter) = "true";  % Set arrows as 'true' for shortest paths
    else
        arrows_l2(iter) = "false";  % Set arrows as 'false' for significant shortest paths
    end
end

% Delete 0 or negative correlations
src_names_edges_cross_1(WEIGHT_ARENA1==0)=[];
trg_names_edges_cross_1(WEIGHT_ARENA1==0)=[];
COLORS1(WEIGHT_ARENA1==0)=[];
arrows_l1(WEIGHT_ARENA1==0)=[];
WEIGHT_ARENA1(WEIGHT_ARENA1==0)=[];

src_names_edges_cross_1(WEIGHT_ARENA1<0)=[];
trg_names_edges_cross_1(WEIGHT_ARENA1<0)=[];
COLORS1(WEIGHT_ARENA1<0)=[];
arrows_l1(WEIGHT_ARENA1<0)=[];
WEIGHT_ARENA1(WEIGHT_ARENA1<0)=[];

%%%%%%%%
src_names_edges_longi_1_2(WEIGHT_ARENA1_2==0)=[];
trg_names_edges_longi_1_2(WEIGHT_ARENA1_2==0)=[];
COLORS1_2(WEIGHT_ARENA1_2==0)=[];
arrows_l1_2(WEIGHT_ARENA1_2==0)=[];
WEIGHT_ARENA1_2(WEIGHT_ARENA1_2==0)=[];

src_names_edges_longi_1_2(WEIGHT_ARENA1_2<0)=[];
trg_names_edges_longi_1_2(WEIGHT_ARENA1_2<0)=[];
COLORS1_2(WEIGHT_ARENA1_2<0)=[];
arrows_l1_2(WEIGHT_ARENA1_2<0)=[];
WEIGHT_ARENA1_2(WEIGHT_ARENA1_2<0)=[];

%%%%%%%%
src_names_edges_cross_2(WEIGHT_ARENA2==0)=[];
trg_names_edges_cross_2(WEIGHT_ARENA2==0)=[];
COLORS2(WEIGHT_ARENA2==0)=[];
arrows_l2(WEIGHT_ARENA2==0)=[];
WEIGHT_ARENA2(WEIGHT_ARENA2==0)=[];

src_names_edges_cross_2(WEIGHT_ARENA2<0)=[];
trg_names_edges_cross_2(WEIGHT_ARENA2<0)=[];
COLORS2(WEIGHT_ARENA2<0)=[];
arrows_l2(WEIGHT_ARENA2<0)=[];
WEIGHT_ARENA2(WEIGHT_ARENA2<0)=[];




end


%% json
% This section creates the structure that will then be exported as a json. 
for iter = 1:1



%% field 3 layers
for iter = 1:1

field3 = 'layers';

c1.name = 'Layer1';
c1.position_x = '-200';
c1.position_y = '0';
c1.position_z = '0';
c1.last_layer_scale = '1';
c1.rotation_x = '0';
c1.rotation_y = '0';
c1.rotation_z = '0';
c1.floor_current_color = '#777777';
c1.geometry_parameters_width = '1001.90476190476';


c2.name = 'Layer2';
c2.position_x = '400';
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
for iter = 1:dim_cross
    nodes1(iter).name = nodes_names_layer_1(iter);
    %nodes1(iter).labelPosition = "Top";
    nodes1(iter).labelPosition = labelPosition(iter);
    %nodes1(iter).labelAdjustX = '0';
    nodes1(iter).labelAdjustX = labelAdjustX(iter);
    nodes1(iter).labelAdjustY = '0';
    nodes1(iter).isBold = bold_no_corr_str(iter);
    nodes1(iter).layer = 'Layer1';
    nodes1(iter).position_x = '0';
    nodes1(iter).position_y = coordinates_layer_1(iter,2);
    nodes1(iter).position_z = coordinates_layer_1(iter,1);
    nodes1(iter).scale_x = nodes_size_layer_1(iter);
    nodes1(iter).color = nodes_color_layer_1(iter);
    nodes1(iter).url = '';
    nodes1(iter).descr = '';
    nodes1(iter).unic_name = nodes_unic_names_layer_1(iter);
    nodes1(iter).labelSize = label_size(iter);
end

% nodes 2 %%%%%%%%%
for iter = 1:dim_cross
    nodes2(iter).name = nodes_names_layer_2(iter);
    %nodes2(iter).labelPosition = "Top";
    nodes2(iter).labelPosition = labelPosition(iter);
    %nodes2(iter).labelAdjustX = '0';
    nodes2(iter).labelAdjustX = labelAdjustX(iter+dim_cross);
    nodes2(iter).labelAdjustY = '0';
    nodes2(iter).isBold = bold_no_corr_str(iter+dim_cross);
    nodes2(iter).layer = 'Layer2';
    nodes2(iter).position_x = '0';
    nodes2(iter).position_y = coordinates_layer_2(iter,2);
    nodes2(iter).position_z = coordinates_layer_2(iter,1);
    nodes2(iter).scale_x = nodes_size_layer_2(iter);
    nodes2(iter).color = nodes_color_layer_2(iter);
    nodes2(iter).url = '';
    nodes2(iter).descr = '';
    nodes2(iter).unic_name = nodes_unic_names_layer_2(iter);
    nodes2(iter).labelSize = label_size(iter+dim_cross);
end

% figure border and quadrants L1
for iter = 1:1


m = 0;
m = max(size(nodes_names_layer_1));

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


% figure border and quadrants L2
for iter = 1:1


m = 0;
m = max(size(nodes_names_layer_2));

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


value4 = [nodes1,nodes2];

end


%% field 5 edges
for iter = 1:1
    
field5 = 'edges';

% edges 1 %%%%%%%%%
for iter = 1:max(size(src_names_edges_cross_1))
    edges_1(iter).src = src_names_edges_cross_1(iter);
    edges_1(iter).trg = trg_names_edges_cross_1(iter);
    edges_1(iter).opacity = '0.3';
    edges_1(iter).color = COLORS1(iter);
    edges_1(iter).channel = '';
    edges_1(iter).size = string(WEIGHT_ARENA1(iter)*15);
    %edges_1(iter).size = string(log(sp_size_l1(iter)+1));
    edges_1(iter).arrow = arrows_l1(iter);
end

% figure border and quadrants L1
m = 0;
m = max(size(src_names_edges_cross_1));

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
for iter = 1:max(size(src_names_edges_longi_1_2))
    edges_1_2(iter).src = src_names_edges_longi_1_2(iter);
    edges_1_2(iter).trg = trg_names_edges_longi_1_2(iter);
    edges_1_2(iter).opacity = '0.3';
    edges_1_2(iter).color = COLORS1_2(iter);
    edges_1_2(iter).channel = '';
    edges_1_2(iter).size = string(WEIGHT_ARENA1_2(iter)*15);
    %edges_1_2(iter).size = string(log(sp_size_l1_2(iter)+1));
    edges_1_2(iter).arrow = arrows_l1_2(iter);
end


% edges 2 %%%%%%%%%
for iter = 1:max(size(src_names_edges_cross_2))
    edges_2(iter).src = src_names_edges_cross_2(iter);
    edges_2(iter).trg = trg_names_edges_cross_2(iter);
    edges_2(iter).opacity = '0.3';
    edges_2(iter).color = COLORS2(iter);
    edges_2(iter).channel = '';
    edges_2(iter).size = string(WEIGHT_ARENA2(iter)*15);
    %edges_2(iter).size = string(log(sp_size_l2(iter)+1));
    edges_2(iter).arrow = arrows_l2(iter);
end

% figure border and quadrants L2
m = 0;
m = max(size(src_names_edges_cross_2));

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


value5 = [edges_1,edges_1_2,edges_2];


end


s = struct(field3,value3,field4,value4,field5,value5);
%s = struct(field3,value3,field4,value4,field5,value5);
j = jsonencode(s);
file = prettyjson(j)



end


%% open mlnetwork
% this section causes a browser (we recommend chrome) to automatically open to the MLNetwork page with the figure displayed.
for iter = 1:1

httpsUrl = 'https://dev.mlnetwork-diplab.ch/cache.php';
options = weboptions('RequestMethod','POST');
responseEmployee = webwrite(httpsUrl, file, options);
web(responseEmployee);

end