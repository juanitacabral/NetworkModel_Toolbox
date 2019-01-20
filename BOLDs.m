function b = BOLDs(T,dt,r)

% Transform the simulated LFP activity into BOLD signal
%
% used in Cabral, Hugues, Sporns and Deco (2011)
% 
% b = BOLDs(T,dt,r)
%
% T       : total time (s)
% dt      : time step (s)
% r       : Simulated LFP activity
%
% b       : BOLD signal

% BOLD model parameters

taus   = 0.65; % 0.8;    % time unit (s)
tauf   = 0.41; % 0.4;    % time unit (s)
tauo   = 0.98; % 1;      % mean transit time (s)
alpha  = 0.32; % 0.2;    % a stiffness exponent
itaus  = 1/taus;
itauf  = 1/tauf;
itauo  = 1/tauo;
ialpha = 1/alpha;
Eo     = 0.34; % 0.8;    % resting oxygen extraction fraction
vo     = 0.02;
k1     = 7*Eo; 
k2     = 2; 
k3     = 2*Eo-0.2;

% Time integration

n_t    = round(T/dt)+1;
x      = zeros(n_t,4);
x(1,:) = [0 1 1 1];

for n = 1:n_t-1
    x(n+1,1) = x(n,1) + dt*( r(n)-itaus*x(n,1)-itauf*(x(n,2)-1) );
    x(n+1,2) = x(n,2) + dt*x(n,1);
    x(n+1,3) = x(n,3) + dt*itauo*(x(n,2)-x(n,3)^ialpha);
    x(n+1,4) = x(n,4) + dt*itauo*(x(n,2)*(1-(1-Eo)^(1/x(n,2)))/Eo - (x(n,3)^ialpha)*x(n,4)/x(n,3));
end

b  = 100/Eo*vo*( k1.*(1-x(:,4)) + k2*(1-x(:,4)./x(:,3)) + k3*(1-x(:,3)) );

clear x;

t_min = 20;
n_min = round(t_min/dt);

b     = b(n_min:end);
