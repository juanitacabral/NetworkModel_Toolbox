function [Phases_Save, dt_save] = Kuramoto_Delays_Run_AAL(C,D,f,K,MD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to run simulations of coupled systems using the
% KURAMOTO MODEL OF COUPLED OSCILLATORS WITH TIME DELAYS
%
% - Each node is represented by a Phase Oscillator 
% - with Natural Frequency f
% - coupled according to a  connectivity matrix C
% - with a distance matrix D (scaled to delays by the mean delay MD)
%  
% All units in Seconds, Meters, Hz
%
% Simulation code written by Joana Cabral, 
% January 2019, joana.cabral@psych.ox.ac.uk
%
%%%%%%%%%%%%%%%%%%%%

%% MODEL PARAMETERS

% % Define simulation Parameters
% f=40; % Node natural frequency in Hz (i.e. f=40))
% MD=0.02; % Mean Delay in seconds (i.e. MD=0.01)
% K=5;
dt=1e-4; % Resolution of model integration in seconds (i.e. dt=1e-4)
tmax=5; % Total simulation time in seconds
t_prev=0; % Total preview time in seconds
dt_save=2e-3; % Resolution of saved simulations in seconds
noise=0;

% Normalize parameters
N=size(C,1);   % Number of units
Omegas=2*pi*f*ones(N,1)*dt; % Frequency of all units in radians/second
kC=K*C*dt; % Scale matrix C with 'K' and 'dt' to avoid doing it at each step
dsig = sqrt(dt)*noise; % normalize std of noise with dt

% Set a matrix of Delays containing the number of time-steps between nodes
% Delays are integer numbers, so make sure the dt is much smaller than the
% smallest delay.
if MD==0
    Delays=zeros(N); % Set all delays to 0.
else
    Delays=(round((D/mean(D(C>0))*MD)/dt));
end
Delays(C==0)=0;

Max_History=max(Delays(:))+1;
% Revert the Delays matrix such that it contains the index of the History 
% that we need to retrieve at each dt
Delays=Max_History-Delays;


%% Initialization

% History of Phases is needed at dt resolution for as long as the longest 
% delay. The system is initialized all desinchronized and uncoupled
% oscillating at their natural frequencies

    Phases_History = 2*pi*rand(N,1)+Omegas.*ones(N,Max_History).*(1:Max_History);
    Phases_History = mod(Phases_History,2*pi); % all values between 0 and 2pi
    
% This History matrix will be continuously updated (sliding history).
% figure; plot(sin(Phases_History'))   % To see the History

% Simulated activity will be saved only at dt_save resolution
Phases_Save=zeros(N,tmax/dt_save); 
sumz=zeros(N,1);

%% Run simulations
disp(['Now running for K=' num2str(K) ', mean Delay = ' num2str(MD*1e3) 'ms'])

tic % to count the time of each run
nt=0;

for t=0:dt:t_prev+tmax  % We only start saving after t_prev
    
    Phase_Now=Phases_History(:,end); % The last collumn of Z is 'now'
    
    % Input from coupled units
    for n=1:N
        sumzn=0; % Intitalize total coupling received into node n
        for p=1:N
            if kC(n,p) % Calculate only input from connected nodes (faster)             
                sumzn = sumzn + kC(n,p) * sin(Phases_History(p,Delays(n,p))-Phase_Now(n));
            end  
        end
        sumz(n)=sumzn;
    end
    
    if MD>0 % Slide History only if the delays are >0
        Phases_History(:,1:end-1)=Phases_History(:,2:end);
    end
    
    % Update the end of History with evolution in the last time step           
    Phases_History(:,end)= Phase_Now + Omegas + sumz + dsig*randn(N,1);
    
    % Save dt_save resolution after t_prev
    if ~mod(t,dt_save) && t>t_prev
        nt=nt+1;
        Phases_Save(:,nt)=Phases_History(:,end);
    end
end

disp(['Finished, lasted ' num2str(round(toc)) ' secs for real ' num2str(t_prev+tmax) ' seconds'])
disp(['Simu speed ' num2str(round(toc/(t_prev+tmax))) ' Realseconds/SimuSecond'])

