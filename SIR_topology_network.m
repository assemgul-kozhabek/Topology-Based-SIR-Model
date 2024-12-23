% Dr. Wei Koong Chai, Dr. Assemgul Kozhabek

function [Sfrac, Ifrac, Rfrac,V, t] = SIR_topology_network(time, adj, beta, delta, i0)

%%%%%% Initialisation %%%%%%%
N=length(adj);            % number of nodes

% Creating initial infected unique nodes (initial conditions)
Io=zeros(1,N);

% Seed nodes
% MetrLa Rho0.4, io=2 
%Io(15)=1;
%Io(39)=1;

Ro = zeros(1,N);
So = ones(1,N) - Io;
Vo = [So, Io, Ro];

%set an error
options=odeset('RelTol',1e-4); 

%timespan
tspan = [0,time]; 

%call the solver
[t,V] = ode15s(@TestFunction, tspan, Vo, options, beta, delta, adj, N); 

%Scale
%%ids=[];

%normalisation of each line from 0 to 1
%di=[]; 
%dr=[]; 
%ds=[];

%plot the results
figure 
hold on
V=V';

xlswrite('outputAdjMetrLa.xlsx', V);

Sfrac = sum(V(1:N, 1:end))/N;
Ifrac = sum(V(N+1:2*N, 1:end))/N;
Rfrac = sum(V(2*N+1:3*N, 1:end))/N;

plot(t, Sfrac, 'b^', 'LineWidth', 0.1);
plot(t, Ifrac, 'rx', 'LineWidth', 0.1);
plot(t, Rfrac, 'go', 'LineWidth', 0.1);
plot(ids, ds, 'b--', 'LineWidth', 3);
plot(ids, di, 'r-', 'LineWidth', 3);
plot(ids, dr, 'g.', 'LineWidth', 4);


set(gca,'XTickLabel',{'06:00','07:30','09:00','10:30','12:00','13:30'})


xlabel('Time')
ylabel('Fraction of the nodes')

legend('S_{model}', 'I_{model}', 'R_{model}', 'S_{data}', 'I_{data}', 'R_{data}');

return

function [dV_dt]= TestFunction(~, V, beta, delta, adj, N) 

S = V(1:N);  % dimension: Nx1
I = V(N+1:2*N);
R = V(2*N+1:3*N);

dS_dt = -beta*diag(adj*I)*S; 
dI_dt = beta*diag(adj*I)*S -delta*I; 
dR_dt = delta*I;

dV_dt = [dS_dt; dI_dt; dR_dt];

return

