function varargout = run_model(mode,T,task,varargin)

global N Ntilde n Nout ...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr ...

%% Figure

fh = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
    'NumberTitle','off','ToolBar','none');

ah1 = axes(fh,'LineWidth',2,'FontSize',12,'Position',[0.1 0.7 0.8 0.2],...
    'ylim',[-1.1 1.1]);
ah2 = axes(fh,'LineWidth',2,'FontSize',12,'Position',[0.1 0.1 0.8 0.5]);
xlabel(ah2,'time (s)');

delete(findobj(ah1,'Type','Line'));
delete(findobj(ah2,'Type','Line'));

lh(1) = line(ah1,'Color','b','LineWidth',1,'Marker','none','LineStyle','-');
lh(2) = line(ah1,'Color','k','LineWidth',1,'Marker','none','LineStyle','-');

switch mode
    
    case 'PCA'
        
        nplt = 5;
        for i = 1:nplt
            lhr(i) = line(ah2,'Color','r','LineWidth',2,'Marker','none','LineStyle','-');
        end
        
    case {'demean','RLStrain','RLStest'}
        
        nplt = 5;
        for i = 1:nplt
            lhV(i) = line(ah2,'Color','k','LineWidth',2,'Marker','none','LineStyle','-');
        end
        lh(3) = line(ah1,'Color','r','LineWidth',1,'Marker','none','LineStyle','-');
        
end

%% Unpack varargin

if any(strcmp(mode,{'demean','RLStrain','RLStest'}))
    
    uPC = varargin{1};
    
    if any(strcmp(mode,{'RLStrain','RLStest'}))
        
        Imu = varargin{2};
        
        switch mode
            
            case 'RLStrain'
                
                if numel(varargin) > 2
                    
                    dJ_PC = varargin{3};
                    w = varargin{4};
                    P = varargin{5};
                    
                end
                
            case 'RLStest'
                
                dJ_PC = varargin{3};
                w = varargin{4};
            
        end
        
    end
    
end

%% Initialize matrices for saving data

if strcmp(mode,'PCA')
    
    %for saving covariance of firing rate network inputs
    COV_fJ = zeros(Ntilde);
    
elseif strcmp(mode,'demean')
    
    %computed mean of input into each spiking neuron due to learned signal
    sum_Imu = zeros(N,1);
    Imu = zeros(N,1);
    
    %mean firing rate, saved for each trial (in demean mode)
    mFR = NaN(T,1);
    
elseif any(strcmp(mode,{'RLStrain','RLStest'}))
    
    %nMSE, saved for each trial
    nMSE = NaN(T,1);   
    
    if any(strcmp(mode,{'RLStrain'}))
        
        %dimension of PC space
        ntilde = size(uPC,2);
        
        if strcmp(mode,'RLStrain') && numel(varargin) < 5
       
            %learned spike net matrix
            dJ_PC = zeros(ntilde,2*n);
            w = zeros(Nout,2*n); %output matrix
            P = 1 * eye(2*n); %inverse covariance matrix
        end
        
    end
    
end

%% Initalize state variables

%place network into the same state
rng(3);

x = 1e0 * randn(Ntilde,1); %teaching network state
%time in current trial
ttrial = inf;
%total time in current trial
TTrial = 0;
%trial number
t_trial_num = 0;
%time across all trials
ttime = 1;
%flag to quit loop
go = true;

if any(strcmp(mode,{'demean','RLStrain','RLStest'}))
    
    %product of PCs to Ntilde, Ntilde to N. For going directory from PCs to
    %spiking network
    uJ_uPC = uJ * uPC;
    
    %spike net state
    V = 1e-3 * randn(N,1);
    
    %learned input current
    zJ = zeros(N,1);
    %initial slow input current
    z0s = zeros(N,1);
    %inital fast input current
    z0f = zeros(N,1);
    %initial fast inhibition
    muglo = zeros(N,1);
        
    if any(strcmp(mode,{'RLStrain','RLStest'}))
        
        %slow presynaptic current when using RLS
        ss = zeros(n,1);
        %fast presynaptic current when using RLS
        sf = zeros(n,1);
      
    end
    
end

%% Run simulation

while go
    
    %generate trial data
    if ttrial > TTrial
        
        %reset time for this trial to 1
        ttrial = 1;  
        
        %increase trial count by 1
        t_trial_num = t_trial_num + 1;
        
        switch task
            case 'osc4'
                [fin,fout,TTrial] = trial_osc4(dt);
        end
        
        %fout input into each rate unit
        u_tilde_fout_fout = u_tilde_fout * fout;
        %fin input into each rate unit
        u_tilde_fin_fin = u_tilde_fin * fin;
        
        if any(strcmp(mode,{'demean','RLStrain','RLStest'}))
            %fin input into each spiking neuron
            uJ_utilde_fin_fin = uJ * u_tilde_fin_fin;
        end
        
        if t_trial_num > Tinit
            
            %Make temporary empty arrays for saving things, or trial specific
            %inputs
            switch mode
                
                case 'PCA'
                    
                    %for collecting fJ of the rate network, for doing PCA
                    fJs = zeros(Ntilde,TTrial);
                    
                case {'demean','RLStrain','RLStest'}
                    
                    %for collecting z(t), for plotting
                    zs = zeros(Nout,TTrial);
                    %for collecting V(t), for plotting
                    Vs = zeros(nplt,TTrial);
                    
                    if strcmp(mode,'demean')
                        %for collecting zJ,z0f and z0s, for mean
                        %computation
                        all_zs = zeros(N,TTrial);
                        %for collecting the spikes, for computing the mean
                        %firing rate
                        Ss = zeros(N,TTrial);
                        
                    end
                    
            end
            
        end
        
    end
    
    %rate model
    fJ = Jtilde * tanh(x) + u_tilde_fout_fout(:,ttrial);
    xinf = fJ + u_tilde_fin_fin(:,ttrial);
    x = xinf + (x - xinf) * etaux;
    
    %spiking model
    if any(strcmp(mode,{'demean','RLStrain','RLStest'}))
        
        %project firing rate network inputs down to the PC basis
        fPC = uPC' * fJ;
        
        Vinf = zJ + z0s + z0f - Imu ...
            + muglo + uJ_utilde_fin_fin(:,ttrial);
        V = Vinf + (V - Vinf) * etaum;
        
        S = V >= Vth;
        V(S) = Vth + Vr;
        
        z0s = z0s * etaus + sum(Js(:,S),2);
        z0f = z0f * etauf + sum(Jf(:,S),2);
        muglo = muglo * etauf + sum(Jmu(:,S),2);
        
        %make presynaptic currents and concatenate
        if any(strcmp(mode,{'RLStrain','RLStest'}))
            
            ss = etaus * ss + S(s_train);
            sf = etauf * sf + S(s_train);
            
            s = [ss;sf];
            
        end
        
        %generate output and learned feedback output
        if t_trial_num > Tinit
                        
            switch mode
                
                case 'demean'

                    %use target as feedback 
                    zPC = fPC;
                    %project from PC to spiking basis
                    zJ = uJ_uPC * fPC;
                    
                    %output is just the target
                    z = fout(:,ttrial);
                    
                case {'RLStrain','RLStest'}
                    
                    %feedback currents in PC basis
                    zPC = dJ_PC * s;
                    %project from PC to spiking
                    zJ = uJ_uPC * zPC;
                    
                    %output
                    z = w * s;
             
            end
            
        else
            
            %feedback is targets in PC basis
            zPC = fPC;
            %project from PC basis to spiking basis
            zJ = uJ_uPC * zPC;
            
            %output
            z = fout(:,ttrial);
            
        end
        
    end
    
    %after initial period
    if t_trial_num > Tinit
        
        %after initial period, save certain results for different computations
        switch mode
            
            case 'PCA'
                
                %save for plotting and computing cov matrix
                fJs(:,ttrial) = fJ;
                
            case {'RLStrain','demean','RLStest'}
                
                %save for plotting and computing nMSE
                Vs(:,ttrial) = V(1:nplt);
                Vs(S(1:nplt),ttrial) = 2 * (Vth - Vr);
                zs(:,ttrial) = z;
                
                if strcmp(mode,'RLStrain')
                    
                    %do RLS
                    if rand < 1/DTRLS
                        
                        xP = P * s;
                        k = (1 + s' * xP)\xP';
                        
                        P = P - xP * k;
                        dJ_PC = dJ_PC - (zPC - fPC) * k;
                        w = w - (z - fout(:,ttrial)) * k;
                        
                    end
                    
                elseif strcmp(mode,'demean')
                    
                    %save for mean computation, project targets from PC
                    %basis to spiking basis, and add initial slow and fast
                    %currents
                    all_zs(:,ttrial) = uJ_uPC * fPC + z0s + z0f;
                    %save for mean firing rate computation
                    Ss(:,ttrial) = S;
                       
                end
                
        end
        
        % text output, plot things, perform computations
        if ttrial == TTrial
            
            %perform computations
            switch mode
                               
                case 'PCA'
                    
                    %update covariance matrix
                    COV_fJ = COV_fJ + fJs * fJs';
                    
                    %compute PCs from cov matrix, and arrange them from
                    %largest to smallest
                    [uPC,d] = eig(COV_fJ);
                    uPC = fliplr(uPC);
                    d = flipud(diag(d));
                    
                    %pick the number of dimensions that accounts for fracPCA of variance
                    ntilde = find(cumsum(d)/sum(d) > fracPCA,1, 'first');
                    uPC = uPC(:,1:ntilde);
                    
                case 'demean'
                    
                    %add sum of current trial's currents to running sum
                    sum_Imu = sum_Imu + sum(all_zs,2);
                    %divde by total elapsed time to compute mean
                    Imu = sum_Imu/ttime;
                    
                    %sum all spikes on this trial, and normalize by number
                    %of neurons and elapsed time on this trial to compute
                    %mean population firing rate on this trial
                    mFR(t_trial_num-Tinit) = sum(Ss(:))/(N*TTrial*dt);
                    
                case {'RLStrain','RLStest'}
                    
                    %compute normalized output error on this trial
                    nMSE(t_trial_num-Tinit) = sum(diag((zs - fout) * (zs - fout)'))/...
                        sum(diag(fout * fout'));
                    
            end
            
            %printing and plotting
            clc;
            
            %plot fout and fin
            set(lh(1),'XData',dt*[ttrial-TTrial+1:ttrial],'YData',fout(1,:));
            set(lh(2),'XData',dt*[ttrial-TTrial+1:ttrial],'YData',fin(1,:));
            
            switch mode
                
                case 'PCA'
                    
                    fprintf('%s, %g PCs, %g trials of %g \n', mode, ntilde, t_trial_num-Tinit, T);
                    
                    %project down to the PCs
                    temp = uPC' * fJs;
                    max_temp = 0;
                    %plot a few of them
                    for i = 1:nplt
                        max_temp = max_temp + abs(min(temp(i,:)));
                        set(lhr(i),'XData',dt*[ttrial-TTrial+1:ttrial],'YData', temp(i,:) + max_temp);
                        max_temp = max_temp + max(temp(i,:));
                    end
                    axis tight
                    
                case {'demean','RLStrain','RLStest'}
                    
                    if strcmp(mode,{'demean'})
                        
                        fprintf('%s, %g trials of %g \n', mode, t_trial_num-Tinit, T);
                        
                    else
                        
                        %print median (across trials) nMSE
                        fprintf('%s, %g Error, %g trials of %g \n', ...
                            mode, nanmedian(nMSE), t_trial_num-Tinit, T);                       
                        
                    end
                    
                    %plot generated output
                    set(lh(3),'XData',dt*[ttrial-TTrial+1:ttrial],'YData',zs(1,:));
                    
                    %plot some voltage trajectories
                    maxV = 0;
                    for i = 1:nplt
                        maxV = maxV + abs(min(Vs(i,:)));
                        set(lhV(i),'XData',dt*[ttrial-TTrial+1:ttrial],'YData', Vs(i,:) + maxV);
                        maxV = maxV + max(Vs(i,:));
                    end
                    axis tight
                    
            end
            
            drawnow;
                          
        end
        
        %counter of number of timesteps that have passed in total, over all
        %trials, only starts counting after initial period is passed
        ttime = ttime + 1;
        
    end
    
    %quit simulation loop
    if t_trial_num == T+Tinit && ttrial == TTrial
        
        %quit loop
        go = false;
        
    end
    
    %counter for number of timesteps that have passed in THIS trial
    ttrial = ttrial + 1;
    
end

%% output
switch mode
    
    case 'PCA'
        
        varargout{1} = uPC;
        
    case 'demean'
        
        varargout{1} = Imu;
        varargout{2} = mFR;
        
    case 'RLStrain'
        
        varargout{1} = dJ_PC;
        varargout{2} = w;
        varargout{3} = P;
        
    case {'RLStest'}
        
        varargout{1} = nMSE;
        
end

close(fh);
