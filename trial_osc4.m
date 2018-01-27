function [fin,fout,TTrial,cond,ttrialb,ITIb,ITI] = trial_osc4(dt)
%%

%timing of event sequence, for each trial
event = [0.05,0.95] * (1/dt);

%length of a trial
TTrial = sum(event);

%inter-trial interval basline
ITIb = 0;     
ITI = 0;

%trial type
cond = 1;
%start and stop of trial for collecting data
ttrialb = [1,TTrial];

%freq of oscillations in target
freq = [1,2,3,5];

%input pulse
fin1 = 0.5 * ones(1,event(1));
fin2 = zeros(1,sum(event(2)));

%target output
fout1 = 1/max(sum(sin(2*pi*freq'*linspace(0,sum(event(1:2))-1,sum(event(1:2)))*dt))) * ...
    sum(sin(2*pi*freq'*linspace(0,sum(event(1:2))-1,sum(event(1:2)))*dt));

%concatenate various elements into fin and fout
fin = [fin1, fin2];
fout = fout1;

