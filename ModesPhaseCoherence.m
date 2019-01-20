%% Main Modes of Phase Coherence

V1=zeros(size(Phases_Save));

% Get the eigenvectors of Phase Coherence
for t=1:size(PC,1)
    [V1(:,t), ~]=eigs(squeeze(PC(t,:,:)),1);
    if t>2 && dot(V1(:,t),V1(:,t-1))<0
        V1(:,t)=-V1(:,t);
    end   
end

% Select Number of Clusters
% Here number of Eigenvectors needed to represent 70% of the Covariance
Mat=cov(V1');
[V, D]=eig(Mat);
PropEigV= sort(diag(D),'descend')/sum(diag(D));
K=find(cumsum(PropEigV)>0.60,1);

figure
imagesc(0:dt_save:(size(Phases_Save,2)-1)*dt_save,1:N,V1(Order,:))
ylabel('Brain area number')
xlabel('Time (seconds')
title('Leading Eigenvector of Phase Coherence')

mink=7;
maxk=7;
rangeK=mink:maxk;

% Set the parameters for Kmeans clustering
Kmeans_results=cell(size(rangeK));

for k=1:length(rangeK)
    disp(['- ' num2str(rangeK(k)) ' PC states'])
    [IDX, C, SUMD, D]=kmeans(V1',rangeK(k),'Distance','cosine','Replicates',20,'MaxIter',200,'Display','final','Options',statset('UseParallel',0));
    [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
    [~,idx_sort]=sort(ind_sort,'ascend');
    Kmeans_results{k}.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectors
    Kmeans_results{k}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
    Kmeans_results{k}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
    Kmeans_results{k}.D=D(ind_sort);       % Distance from each point to every centroid   end
end

load aal_cog.txt
[x,y,z] = sphere;
scale=3;
x=scale*x;
y=scale*y;
z=scale*z;

figure

for c=1:K
    V_cluster=Kmeans_results{1, 1}.C(c,:);
    V_cluster=V_cluster/max(abs(V_cluster(:)));
    
    subplot(1,K,c)
    hold on
    for n=1:N
        if V_cluster(n)>=0
                surf(x+aal_cog(n,1), y+aal_cog(n,2),z+aal_cog(n,3),'FaceColor','b','EdgeColor','none','FaceAlpha',0.7);
        else
                surf(x+aal_cog(n,1), y+aal_cog(n,2),z+aal_cog(n,3),'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
        end
    end
end



