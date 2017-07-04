%% toy model oscillatory network to illustrate main points
% based on Matsuoka 1985,1987 Biol Cybern
%
% Mark Humphries 16/3/16

clear all; close all

% model parameters - per neuron
N = 3*10; % one triplet per example (original, IC2, IC3, noise 1, noise 2, noise 3, transient perturb, sustained perturb, decaying input, no connections) 
pars.Ta = zeros(N,1) + 0.025; % neuron time-constant
pars.Ty = zeros(N,1) + 0.2;  % adaptation
pars.b = zeros(N,1) + 2;
pars.c = zeros(N,1) + 3; % constant input   
pars.D = zeros(N,1);  % noise amplitude

pars.D(10:12) = 200;  % pair with low noise
pars.D(13:15) = 1000;  % pair with medium noise
pars.D(16:18) = 2000;  % pair with highnoise


% a chain of reciprocal inhibition
% weights = zeros(N-1,1) - 1.5;
% weights(3:3:N-1) = 0;
% pars.w = zeros(N) + diag(weights,1) + diag(weights,-1); 

% a loop
weights = zeros(3) - 1.5;
pars.w = blkdiag(weights,weights,weights,weights,weights,weights,weights,weights,weights,zeros(3));
pars.w(eye(N)==1) = 0;

pars.d = 0.001;  % lag in milliseconds
sim.dt = 0.001;  % 1 millisecond

sim.T = 3; 
nsteps = round(sim.T/sim.dt);
lagsteps = round(pars.d/sim.dt);


sim.tsperturb = 1500:1600;

% sim.tsD = repmat(pars.D,1,nsteps);
% sim.tsD(9:10,round(nsteps/3):round(nsteps/3)+200) = 2000;
sim.tsC = repmat(pars.c,1,nsteps);
sim.tsC(19:21,sim.tsperturb) = 5;  % transient perturbation
sim.tsC(22,round(nsteps/2):end) = 5;  % permenant perturbation
sim.tsC(23,round(nsteps/2):end) = 1;  % permenant perturbation

sim.tsC(25:27,:) = repmat(linspace(5,0,nsteps),3,1); % decaying input


%% initialise model
% giving initial conditions that are NOT EQUAL
a = zeros(N,nsteps); r = a; y = a;
ICs = [1 2 0, 0 0.5 3,5 2 1,1 2 0,1 2 0,1 2 0,1 2 0,1 2 0,1 2 0,1 2 0];  % one pair per simulated network
a(:,1:lagsteps) = repmat(ICs',1,lagsteps);

%% run model
for t = lagsteps+1:nsteps
    s(:,t) = sim.tsC(:,t) + pars.w * r(:,t-lagsteps) + pars.D.*sim.dt.*randn(N,1);  % input: constant, other neurons, noise (sim.tsD(:,t))
    a(:,t) = a(:,t-1) + sim.dt .* (-a(:,t-1) + s(:,t) - pars.b.*y(:,t-1))./pars.Ta;  % activation update
    y(:,t) = y(:,t-1) + sim.dt .* (-y(:,t-1) + r(:,t-1)) ./ pars.Ty;
    r(:,t) = a(:,t) .* (a(:,t) > 0);
end


%% decompose into different examples
plotlag = 0;
clr1 = [0.8 0.6 0.2];
clr2 = [0.6 0.3 0.8];
clr3 = [0.3 0.8 0.6];

time = sim.dt:sim.dt:sim.T;

% 1: basic attractor
coeffs1 = pca(r(1:3,:)');
pc1 = r(1:3,:)' * coeffs1(:,1);  % reconstruct PC1
pc2 = r(1:3,:)' * coeffs1(:,2);  % reconstruct PC2

figure
plot(time,r(1,:),'Color',clr1); hold on
plot(time,r(2,:),'Color',clr2);
plot(time,r(3,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot3(r(1,:),r(2,:),r(3,:),'k')
grid on

figure
plot(pc1(2:end),pc2(2:end),'k'); hold on
plot(pc1(2),pc2(2),'k.','Markersize',10)
xlabel('PC 1')
ylabel('PC 2')


% COMPUTE THE MANIFOLD


% 2: different initial conditions
coeffs2 = pca(r(4:6,:)');
pc1IC2 = r(4:6,:)' * coeffs2(:,1);  % reconstruct PC1
pc2IC2 = r(4:6,:)' * coeffs2(:,2);  % reconstruct PC2

coeffs3 = pca(r(7:9,:)');
pc1IC3 = r(7:9,:)' * coeffs3(:,1);  % reconstruct PC1
pc2IC3 = r(7:9,:)' * coeffs3(:,2);  % reconstruct PC2

figure
plot(pc1(2:end),pc2(2:end),'k'); hold on
plot(pc1(2),pc2(2),'k.','Markersize',10)
plot(pc1IC2(2:end),pc2IC2(2:end),'r'); hold on
plot(pc1IC2(2),pc2IC2(2),'r.','Markersize',10)
plot(pc1IC3(2:end),pc2IC3(2:end),'b'); hold on
plot(pc1IC3(2),pc2IC3(2),'b.','Markersize',10)
xlabel('PC 1')
ylabel('PC 2')

% figure
% plot(r(1,2:end-plotlag),r(2,plotlag+2:end),'k'); hold on
% plot(r(1,2),r(2,plotlag+2),'k.','Markersize',10)
% plot(r(3,2:end-plotlag),r(4,plotlag+2:end),'r'); 
% plot(r(3,2),r(4,plotlag+2),'r.','Markersize',10)
% plot(r(5,2:end-plotlag),r(6,plotlag+2:end),'b'); 
% plot(r(5,2),r(6,plotlag+2),'b.','Markersize',10)
% xlabel('Neuron 1')
% ylabel('Neuron 2')

% 3: noise
coeffsN1 = pca(r(10:12,:)');
pc1N1 = r(10:12,:)' * coeffsN1(:,1);  % reconstruct PC1
pc2N1 = r(10:12,:)' * coeffsN1(:,2);  % reconstruct PC2

figure
plot(time,r(10,:),'Color',clr1); hold on
plot(time,r(11,:),'Color',clr2);
plot(time,r(12,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot(pc1N1(2:end),pc2N1(2:end),'k'); hold on
plot(pc1N1(2),pc2N1(2),'k.','Markersize',10)
xlabel('PC 1')
ylabel('PC 2')

coeffsN2 = pca(r(13:15,:)');
pc1N2 = r(13:15,:)' * coeffsN2(:,1);  % reconstruct PC1
pc2N2 = r(13:15,:)' * coeffsN2(:,2);  % reconstruct PC2

figure
plot(time,r(13,:),'Color',clr1); hold on
plot(time,r(14,:),'Color',clr2);
plot(time,r(15,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot(pc1N2(2:end),pc2N2(2:end),'k'); hold on
plot(pc1N2(2),pc2N2(2),'k.','Markersize',10)
xlabel('PC 1')
ylabel('PC 2')

coeffsN3 = pca(r(16:18,:)');
pc1N3 = r(16:18,:)' * coeffsN3(:,1);  % reconstruct PC1
pc2N3 = r(16:18,:)' * coeffsN3(:,2);  % reconstruct PC2

figure
plot(time,r(16,:),'Color',clr1); hold on
plot(time,r(17,:),'Color',clr2);
plot(time,r(18,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot(pc1N3(2:end),pc2N3(2:end),'k'); hold on
plot(pc1N3(2),pc2N3(2),'k.','Markersize',10)
xlabel('PC 1')
ylabel('PC 2')



% 4: transient perturbation
coeffsP = pca(r(19:21,:)');
pc1P = r(19:21,:)' * coeffsP(:,1);  % reconstruct PC1
pc2P = r(19:21,:)' * coeffsP(:,2);  % reconstruct PC2

figure
plot(time,r(19,:),'Color',clr1); hold on
plot(time,r(20,:),'Color',clr2);
plot(time,r(21,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot(pc1P(2:sim.tsperturb(1)),pc2P(2:sim.tsperturb(1)),'k'); hold on
plot(pc1P(2),pc2P(2),'k.','Markersize',10)
plot(pc1P(sim.tsperturb),pc2P(sim.tsperturb),'r')
plot(pc1P(sim.tsperturb(end):end),pc2P(sim.tsperturb(end):end),'Color',[0 0 1]); hold on
xlabel('PC 1')
ylabel('PC 2')

% persistent pertubation
coeffsPP = pca(r(22:24,:)');
pc1PP = r(22:24,:)' * coeffsPP(:,1);  % reconstruct PC1
pc2PP = r(22:24,:)' * coeffsPP(:,2);  % reconstruct PC2

figure
plot(time,r(22,:),'Color',clr1); hold on
plot(time,r(23,:),'Color',clr2);
plot(time,r(24,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot(pc1PP(2:round(nsteps/2)),pc2PP(2:round(nsteps/2)),'k'); hold on
plot(pc1PP(2),pc2PP(2),'k.','Markersize',10)
plot(pc1PP(round(nsteps/2)+1:end),pc2PP(round(nsteps/2)+1:end),'Color',[0 0 1]); hold on
xlabel('PC 1')
ylabel('PC 2')

% decaying input
coeffsDec = pca(r(19:21,:)');
pc1Dec = r(25:27,:)' * coeffsDec(:,1);  % reconstruct PC1
pc2Dec = r(25:27,:)' * coeffsDec(:,2);  % reconstruct PC2

figure
plot(time,r(25,:),'Color',clr1); hold on
plot(time,r(26,:),'Color',clr2);
plot(time,r(27,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')

figure
plot(pc1Dec(2:end),pc2Dec(2:end),'k'); hold on
plot(pc1Dec(2),pc2Dec(2),'k.','Markersize',10)
xlabel('PC 1')
ylabel('PC 2')

% no weights
figure
plot(time,r(28,:),'Color',clr1); hold on
plot(time,r(29,:),'Color',clr2);
plot(time,r(30,:),'Color',clr3);
xlabel('Time (s)')
ylabel('Activity')
title('No connections')


save Toy3Dmodel sim pars r




