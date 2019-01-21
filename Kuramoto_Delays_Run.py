import numpy as np
def Kuramoto_Delays_Run_AAL(C, D, f, K, MD):
    ############################%%%#
    #
    # Function to run simulations of coupled systems using the
    # KURAMOTO MODEL OF COUPLED OSCILLATORS WITH TIME DELAYS
    #
    # - Each node is represented by a Phase Oscillator
    # - with Natural Frequency f
    # - coupled according to a  connectivity matrix C
    # - with a distance matrix D (scaled to delays by the mean delay MD)
    #
    # All units in Seconds, Meters, Hz
    #
    # Simulation code written in Matlab by Joana Cabral,
    # January 2019, joana.cabral@psych.ox.ac.uk
    # 
    # Adapted to Python by Laurent Perrinet Jan 2019
    #
    ##############%%%%%#

    # MODEL PARAMETERS

    # # Define simulation Parameters
    # f=40 # Node natural frequency in Hz (i.e. f=40))
    # MD=0.02 # Mean Delay in seconds (i.e. MD=0.01)
    # K=5;
    dt = 1.e-4 # Resolution of model integration in seconds (i.e. dt=1e-4)
    tmax = 1. # Total simulation time in seconds
    t_prev = 0. # Total preview time in seconds
    dt_save = 2.e-3 # Resolution of saved simulations in seconds
    noise = 0

    # Normalize parameters
    N = C.shape[0] # Number of units
    Omegas = 2*np.pi*f*np.ones(N)*dt # Frequency of all units in radians/second
    kC = K*C*dt # Scale matrix C with 'K' and 'dt' to avoid doing it at each step
    dsig = np.sqrt(dt)*noise # normalize std of noise with dt

    # Set a matrix of Delays containing the number of time-steps between nodes
    # Delays are integer numbers, so make sure the dt is much smaller than the
    # smallest delay.
    if MD==0:
        Delays = np.zeros_like(C) # Set all delays to 0.
    else:
        Delays = D / D[C>0].mean() * MD / dt
        Delays = Delays.astype(np.int)

    Delays[C==0] = 0.

    Max_History = np.max(Delays.ravel()) + 1

    # Revert the Delays matrix such that it contains the index of the History
    # that we need to retrieve at each dt
    Delays = Max_History - Delays


    # Initialization

    # History of Phases is needed at dt resolution for as long as the longest
    # delay. The system is initialized all desinchronized and uncoupled
    # oscillating at their natural frequencies

    Phases_History = 2*np.pi*np.random.rand(N, Max_History)
    Phases_History += Omegas[:, None] #* np.ones(N, Max_History)#*np.arange(Max_History)[:, None]
    Phases_History = np.mod(Phases_History, 2*np.pi) # all values between 0 and 2pi

    # This History matrix will be continuously updated (sliding history).
    # figure plot(sin(Phases_History'))   # To see the History

    # Simulated activity will be saved only at dt_save resolution
    Phases_Save = np.zeros((N, int(tmax/dt_save)))
    sumz = np.zeros(N)

    # Run simulations
    print('Now running for K=',  K,  ', mean Delay = ', MD*1e3, 'ms')
    print('Max_History=',  Max_History, 'steps')


    import time
    tic = time.time() # to count the time of each run
    nt = 0
    import tqdm
    for i_t in tqdm.tqdm(range(int((t_prev+tmax)/dt))):
        # We only start saving after t_prev
        Phase_Now = Phases_History[:, -1] # The last collumn of Z is 'now'
        # Input from coupled units
        for n in range(N):
            sumzn = 0 # Intitalize total coupling received into node n
            for p in range(N):
                if kC[n, p]>0: # Calculate only input from connected nodes (faster)
                    phase = Phases_History[p, Delays[n, p]]-Phase_Now[n]
                    sumzn += kC[n,p] * np.sin(phase)
            sumz[n] = sumzn

        if MD>0: # Slide History only if the delays are >0
            Phases_History[:, :-2] = Phases_History[:, 1:-1]

        # Update the end of History with evolution in the last time step
        Phases_History[:, -1] = Phase_Now
        Phases_History[:, -1] += Omegas + sumz + dsig*np.random.randn(N)

        # Save dt_save resolution after t_prev
        if np.mod(i_t, int(dt_save/dt)) == 0 and i_t*dt>t_prev:
            Phases_Save[:, nt] = Phases_History[:, -1]
            nt += 1

    toc = time.time() - tic
    print('Finished, lasted %.3f ' % toc , ' secs for real %.3f ', t_prev+tmax , ' seconds')

    print('Simu speed ratio %.3f =' % (toc/(t_prev+tmax)), ' Realseconds/SimuSecond')

    return Phases_Save, dt_save
