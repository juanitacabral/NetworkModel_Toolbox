
%%%%%%%%%%
%
% Code to compare simulation results with neuroimaging data:
% MEG: Correlation matrix of frequency-specific power envelopes 
% fMRI: the BOLD signal correlation matrix (static FC)
%
% By Joana Cabral joanacabral@med.uminho.pt
%
%%%%%%%%%

X=sin(Phases_Save); 
[N, w] = size(X);

%% Compare with MEG frequency-specific Envelope FC

% First filter the signals in a given frequency band and
% Obtain the Envelopes using the absolute of the Hilbert transform

freq_range_MEG=[6 15];
X_env_band=zeros(size(X));
for n=1:N
    X_band(n,:) = bandpass(X(n,:)-mean(X(n,:)), freq_range_MEG, 1/dt_save, 0);
    X_env_band(n,:)=abs(hilbert(X_band(n,:)));
end

% Compute the envelope correlation matrix
FC_simu=corrcoef(X_env_band');

figure
colormap(jet)
imagesc(FC_simu,[-1 1])
title(['Simulated ' num2str(freq_range_MEG(1)) '-' num2str(freq_range_MEG(2)) ' Envelope FC'])

% Compare with Empirical FC matrix
% load FCempirical.mat FC_emp
% FC_MEG=squeeze(mean(FC_emp(2:4,:,:),1));
% FC_MEG=FC_MEG(Order,Order);
% clear FC_emp
% fit_env_FC=corr2(FC_simu(Isubdiag),FC_MEG(Isubdiag));


%% Compare with static BOLD FC (NEEDS LONG SIMULATION TIMES)

% First transform simulations into BOLD signal using the Balloon-Windkessel
% model
% Transform into BOLD signal
w_cut=w-(20/dt_save)+2; % it will cut the first 20 seconds in the code
BOLD_X=zeros(N,w_cut);
for n=1:N
    BOLD_X(n,:)=BOLDs(w*dt_save,dt_save,X(n,:));   
end

% Compute the correlatin matrix
FC_BOLD_X=corrcoef(BOLD_X');

figure
imagesc(FC_BOLD_X)
title('Simulated BOLD FC')

% Compare with empirical data
% load AAL_matrices.mat FC_OLF
% FC_emp=FC_OLF;
% clear FC_OLF
% Isubdiag = find(tril(ones(N),-1));
% fitFC=corr2(FC_BOLD_X(Isubdiag),FC_emp(Isubdiag));

%%