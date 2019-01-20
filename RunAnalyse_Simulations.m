%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to run dynamical simulations on a SpaceTime Network Model
% and visualize the results.
%
% Generates figure with:
% - timeseries, power spectra and phase order over time
% - matrices of coupling wights (C), distances (D) and Phase Coherence (PC) 
% 
% Joana Cabral January 2019 joanacabral@med.uminho.pt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define here the Network Spatial and Temporal structure

% Define the strenght of all existing connections between the N nodes
% as a NxN matrix 
load AAL_matrices.mat C

% Number of coupled Units
N=size(C,1); 

% Normalize such that the mean of all non-diagonal elements is 1. 
C=C/mean(C(ones(N)-eye(N)>0));

% Distance between areas
load AAL_matrices.mat D
D=D/1000; % Distance matrix in meters

% Change the order of units for visualization (optional)
Order=[1:2:N N:-2:2]; % This sorts the AAL areas into left and right hemispheres
C=C(Order,Order);
D=D(Order,Order);

%% Define here the parameters of the network model to be manipulated

    % Node natural frequency in Hz 
          f=40;    % (i.e. f=40)):          

    % Mean Delay in seconds 
          MD=0.02;  %(i.e. MD=0.01)             
       
    % Global Coupling strength
          K=5; % 
          
%% Call here the function of the Network Model

[Phases_Save, dt_save] = Kuramoto_Delays_Run(C,D,f,K,MD);


%% Process the simulated data and generate figures

figure('Name',[' K= ' num2str(K) ' and mean delay ' num2str(MD*1e3) ' ms'])

% Plot simulated time series
subplot(3,4,1:3)
plot(0:dt_save:(size(Phases_Save,2)-1)*dt_save,sin(Phases_Save'))
xlabel('Time (seconds)')
ylabel('sin(\theta)')
title(['Coupled ' num2str(f) 'Hz oscillators with for K= ' num2str(K) ' and mean delay ' num2str(MD*1e3) ' ms'])

% Power Spectrum
subplot(3,4,4)
fbins=5000;
freqZ = (0:fbins-1)/(dt_save*fbins);
ind100Hz=find(freqZ==100);
Fourier = fft(sin(Phases_Save),fbins,2);     
PSD=abs(Fourier).^2;
plot(freqZ(1:round(fbins/2)),PSD(:,1:round(fbins/2)))
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('Power')

% Plot the Order Parameter over time
subplot(3,4,5:7)
OP=abs(mean(exp(1i*Phases_Save),1));
plot(0:dt_save:(size(Phases_Save,2)-1)*dt_save,OP)
ylim([0 1])
xlabel('Time (seconds)')
ylabel('Order Parameter')

% Power Spectrum of the Order Parameter
subplot(3,4,8)
Fourier = fft(OP-mean(OP),fbins,2); %% Fourier of Z (complex)    
PSD=abs(Fourier).^2;
PSD=PSD(1:round(fbins/2))/sum(PSD(1:round(fbins/2)));
plot(freqZ(1:round(fbins/2)),PSD)
xlim([0 10])
xlabel('Frequency (Hz)')
ylabel('Power')

% Plot Connectivity and Distance Matrices and mean Phase Coherence

subplot(4,3,10)
colormap(jet)
imagesc(C,'AlphaData',~C==0)
xlabel('Nodes')
ylabel('Nodes')
title('Coupling Matrix')
axis square
colorbar

subplot(4,3,11)
imagesc(D)
xlabel('Nodes')
ylabel('Nodes')
title('Distance Matrix')
axis square
colorbar

subplot(4,3,12)
for n=1:N
    for p=1:N
        PC(:,n,p)=cos(Phases_Save(n,:)-Phases_Save(p,:));
    end
end
imagesc(squeeze(mean(PC,1)),[-1 1])
xlabel('Nodes')
ylabel('Nodes')
title('Mean Phase Coherence Matrix')
axis square
colorbar

%% Video Showing the Phase Coherence Matrix
figure
colormap jet

B0=zeros(N,1);

for t=1:size(PC,1)
    
    subplot(1,2,1)
    imagesc(squeeze(PC(t,:,:)),[-1 1])
    xlabel('Nodes')
    ylabel('Nodes')
    title(['Phase Coherence Matrix t=' num2str(t*dt_save) 's'])
    axis square

    subplot(1,2,2)
    cla
    hold on
    B1=Phases_Save(:,t);
    plot([B0 cos(B1)]',[B0 sin(B1)]','b')
    ylabel('Imag')
    xlabel('Real')
    xlim([-1 1])
    ylim([-1 1])
    axis square
    title('BOLD Phases, \Theta_n(t)')

    pause(0.01)
end

