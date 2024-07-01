function ACROSS_BETW=BTW_SOURCE_TARGET_DIR(WHICH_NETWORK,SOURCE,TARGET)
WHICH_NETWORK(TARGET,SOURCE)=0;
D=1./posweights(WHICH_NETWORK);
[SPL,hops,Pmat] = distance_wei_floyd(D);
ACROSS_BETW=zeros(max(size(WHICH_NETWORK)),1);
for iter1=1:max(size(SOURCE))
    for iter2=1:max(size(TARGET))
        B=retrieve_shortest_path(SOURCE(iter1),TARGET(iter2),hops,Pmat);
        ACROSS_BETW2=zeros(max(size(WHICH_NETWORK)),1);
        ACROSS_BETW2(B)=1;
        ACROSS_BETW=ACROSS_BETW+ACROSS_BETW2;
    end
end
        
        
        