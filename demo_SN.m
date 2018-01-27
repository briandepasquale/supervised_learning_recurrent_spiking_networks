%This script will train a network of LIF neurons based on the method
%described in "Using firing rate dynamics to train recurrent networks of
%spiking model neurons" (https://arxiv.org/abs/1601.07620)

%1. PCs of the firing rate dyanmics are computed, so that only the most
%dominate signals are learned in the spiking network.
%2. The mean input into each spiking neuron is computed and subtracted, to
%control aphysiological firing rates.
%3. RLS is used to train the inputs into each spiking neuron.
%4. Training performance is tested.

global N Ntilde n Nout Nin ...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr
    
%% Number of neurons

%N: number of spiking neurons
%Ntilde: number of rate units

N = 500;
Ntilde = 500;

%% How long to run each routine and which ones to run

%init: number of trials for dynamics to settle from initial condition
%PCA: how many trials to compute PCs
%demean: how many trials to compute mean
%RLS: how many trials to do RLS training for
%test: how many trials to test for, after training

Tinit = 50;
TPCA = 300;
Tdemean = 300;
TRLS = 500;
Ttest = 100;

%% Gain parameters

%gtilde: gain of recurrent connectivity
%gfout: gain of output feedback
%gfin: gain of input signal
%gz: learned feedback gain of spiking net
%muf: global inhibition
%gs: gain of slow random synapses
%gf: gain of fast random synpases

gtilde = 1.4;
gfout = 1.0;
gfin = 1.0;
gz = 6;
muf = -0.3;
%gs = 0.11; 
%gf = 0.13;

gs = 0;
gf = 0;

%% task

%task: Pick the task
%out: dimension of output
%in: dimension of input
%n: Number of s's targets to use to construct z

task = 'osc4';
Nout = 1;
Nin = 1;
n = round(0.4 * N);

%% Parameters

%fracPCA: used to compute number of PCs to use
fracPCA = 0.999;

%DTRLS: number of time steps between RLS updates
%dt: timestep
%taux: decay time for rate net
%tauf: decay for fast synapses
%taus: decay for slow synapses
%taum: decay for voltage in spiking net
%Vth: spike threshold
%Vr: reset voltage

DTRLS = 2; 
dt = 1e-3;
taux = 1e-2;
taus = 1e-1;
tauf = 5e-3;
taum = 1e-2;
Vth = 1e-16;
Vr = -10;

%precompute for doing Euler integration
etaux = exp(-dt/taux);
etaus = exp(-dt/taus);
etauf = exp(-dt/tauf);
etaum = exp(-dt/taum);

%% Random connections

%Jtilde: recurrent firing rate connections
%ufout: input connections into rate network for fout
%ufin: input connecdtions into rate network for fin
%uJ: connections from low D learned signals into spiking network
%Jmu: global recurrent inhibition
%Jf: fast random synapses
%Js: slow random synapses
%s_train: which synapses will be trained in the spiking network

randseed = 6;
rng(randseed); %set ran seed, for ran matrcies
Jtilde = gtilde * 1/sqrt(Ntilde) * randn(Ntilde);
u_tilde_fout  = gfout * (-1 + 2 * rand(Ntilde,Nout));
u_tilde_fin = gfin * (1 + 2 * rand(Ntilde,Nin));

uJ = gz * (orth((-1 + 2 * rand(N,Ntilde))));

Jmu = muf * (1/tauf) * 1/(N) * ones(N);
Jf = gf * (1/tauf) * 1/sqrt(N) * randn(N);
Js = gs * (1/taus) * 1/sqrt(N) * randn(N);

s_train = false(1,N);
s_train(1,randsample(N,n)) = true;

%% Compute PCs of rate model
uPC = run_model('PCA',TPCA,task);

%% demean spiking inputs
[Imu,mFR] = run_model('demean',Tdemean,task,uPC);

%% Train with RLS
[dJ_PC,w] = run_model('RLStrain',TRLS,task,uPC,Imu);

%% Test RLS solution
ERR_RLS = run_model('RLStest',Ttest,task,uPC,Imu,dJ_PC,w);

