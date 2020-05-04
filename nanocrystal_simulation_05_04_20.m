clear, clc, close all % clear variables, clear command window, close figures
tic % start of timing
rng(2) % random seed for consistent results

% initial parameters
dir = 'C:\Users\jcavi\Box Sync\Research\Rohan\Projects\Nanocrystals\MATDAT'; % directory for output 
N_NC = 10%200; % number of nanocrystals
N_sites = 20%20; % number of sites per nanocrystal
N_time =  600000; % maximum time-steps
N_plots = 10; % number of individual trajectories to plot/write; usually keep this as a small number (10)
tol = .999; % fraction of total successful swaps necessary to exit
fnuc = 2/10; % fraction of swaps needed to initialize nucleation (make sure to test that this is working how you want)

k1 = 3/10*10^-1; % scaling factor corresponding to concentration
krange = 1%[1, 2, 4, 8]; % the range of k/k0 used
k1s = k1*krange; % the k values

% bool options
isWrite = false; % determines whether data will be written to a "MATDAT" file
isWriteXL = false; % determines whether concentration-dependent statistics will be written to an excel file 
isRejection = false; % determines whether to use the "exception rule"; adds a multiplicative factor to the probability that makes later swaps more difficult 
isNE = false; % "isNon-Equilibrium" determines whether dG<0 means p=1.00

isPlotK = numel(krange)>1; % if there is a range of k values, timing statistics can be plotted vs k

% important statistics
STmeds = []; % switching time medians
STerrs = []; % switching time standard deviations
WTmeds = []; % waiting time medians
WTerrs = []; % waiting time standard deviations

% which models to use
modeltype = [3,6]; % array of models to use, e.g. diffusion limited and Phase + nucleation
modeltype = 6
% legend for modeltype
%     0: Sudden Nucleation
%     1: Old Nucleation
%     2: Phase + Positive Cooperativity
%     3: Diffusion Limited
%     4: Positive Cooperativity
%     5: New Nucleation
%     6: Phase + Nucleation
%     7: exchange_density
    
if isRejection
    imgdir = [dir,'\rejection'];
    mkdir(imgdir)
else
    imgdir = [dir,'\normal'];
    mkdir(imgdir)
end
dir = imgdir;

% figure parameters
figpos = [1,1,5,4];
fs = 14;

% making figures
Gfig = figure;
pffig = figure;
psfig = figure;

leglabels = {};
for nt = modeltype % loop over models
    for k1=k1s % loop over concentrations
        
        dG0 = -log(k1); %initial free energy from k1
        
        %initialize cell arrays
        labels = cell(0);
        data = cell(0);
        
        % set k and G functions for different models
        if nt==0
            plot_title = 'Sudden Nucleation';
            kf = -3; %final k value
            k = @(i) k1*(i<=N_sites*fnuc+2)+kf*(i>N_sites*fnuc+2); % the +2's here were fudge factors because I never got fnuc to work exactly; just make sure the transition happens where you want
            G = @(i) -log(k(i));
        elseif nt==1
            plot_title = 'Old Nucleation';
            % this is some sort of curve-fitting procedure to make the
            % curve intersect G0 at i=1 and 0 at i=Nc
            Nc = N_sites/5;
            a = (99/k1)^(1/(N_sites/2));
            beta(2) = a;
            beta(1) = Nc^(1/3)*a;
            G = @(i) beta(1)*(i-1).^(2/3)-beta(2)*(i-1)+dG0;
            k = @(i) exp(-(G(i)));
        elseif nt==2
            plot_title = 'Phase + Positive Cooperativity';
            ddG = -1.3154; % slope of dG
            G = @(i) ddG*(i-(N_sites*fnuc+1)).*(i>N_sites*fnuc+1)+dG0;
            k = @(i) exp(-(G(i)));
        elseif nt==3
            plot_title = 'Diffusion Limited';
            G = @(i) dG0*ones(1,numel(i));
            k = @(i) exp(-(G(i)));
        elseif nt==4
            plot_title = 'Positive Cooperativity';
            ddG = -1.3154/3; % slope of dG
            G = @(i) ddG*(i-1)+dG0;
            k = @(i) exp(-(G(i)));
        elseif nt==5
            plot_title = 'New Nucleation';
            nc = N_sites*fnuc*2;
            beta = dG0*(1/(nc^(1/3)-1));
            alpha = dG0 *(1/(-nc^(-1/3)+1));
            G = @(i) -beta+alpha*i.^(-1/3);
            k = @(i) exp(-(G(i)));
        elseif nt==6
            plot_title = 'Phase + Nucleation';
            nc = N_sites*fnuc+1;
            nn = N_sites*fnuc+-1;
            beta = dG0*(1/(nc^(1/3)-1));
            alpha = dG0 *(1/(-nc^(-1/3)+1));

            G = @(i) (alpha*(i-nn+1e-10).^(-1/3)-beta).*(i>nn)+dG0.*(i<=nn);
            k = @(i) exp(-(G(i)));
        elseif nt == 7
            plot_title = 'Exchange-density';
            dE = (4)/N_sites;
            Gf = dG0-dE*N_sites;
            N = N_sites;
            G = @(i) G_parabola( dG0, Gf, N, i );
            k = @(i) exp(-(G(i)));
        end
        type = plot_title;
        if isRejection
            type = [type,'-rejection'];
        end
        type = [type,'-k=',num2str(k1)];
        label = plot_title;
        leglabels{end+1} = type;
        fp = dircat(dir,['mc_05-04-2020',type,'.matdat']); % create output data file name
        
        
        
        p = @(i) k(i)./(1+k(i)); %convert k into p
        
        if isNE % NE condition
            p = @(i) p(i).*(G(i)>0)+1.*(G(i)<=0);
        end
        if isRejection % rejection condition
            p = @(i) p(i).*(1-(i-1)/20);
        end
        
        is = 1:N_sites; % vector of sites #'s

        % add data to cell arrays
        labels = {labels{:},'swap'};
        data = {data{:},is};
        labels = {labels{:},'dG'};
        data = {data{:},G(is)};
        labels = {labels{:},'1-p'};
        data = {data{:},1-p(is)};
        labels = {labels{:},'p'};
        data = {data{:},p(is)};

        % plot free energy vs swap
        figure(Gfig)
        hold on
        plot(is,G(is),'o-')
        title('Free Energy','FontSize',fs)
        ylabel('\Delta G (kT)','FontSize',fs)
        xlabel('nth ion swap','FontSize',fs)
        ylim([-5,6])
        legend(leglabels,'Location','southwest')


        % plot probability of failure vs swap
        figure(pffig)
        hold on
        semilogy(is,(1-p(is)),'o-')
        title('Probability of Failure','FontSize',fs)
        ylabel('probability','FontSize',fs)
        xlabel('nth ion swap','FontSize',fs)
        %ylim([1e-6,1])
        set(gca,'yscale','log')
        
        legend(leglabels,'Location','southwest')
        
        % plot probability of success vs swap
        figure(psfig)
        hold on
        plot(is,(p(is)),'o-')
        title('Probability of Success','FontSize',fs)
        ylabel('probability','FontSize',fs)
        xlabel('nth ion swap','FontSize',fs)
        ylim([0,1])
        
        legend(leglabels,'Location','northwest')


        filename = strrep(plot_title,' ','_');
        t_range = 1:N_time;
        
        nano_sites = zeros(N_NC,N_time);
        WTs = zeros(N_NC,1); % vector of waiting times
        STs = zeros(N_NC,1); % vector of switching times
        Starts = zeros(N_NC,1); % vector to keep track of when switching starts 

        j=0;
        for t = t_range
            j=j+1;
            NC = randi(N_NC);
            nano_sites(:,t+1)=nano_sites(:,t);
            
            if nano_sites(NC,t) < N_sites && p(nano_sites(NC,t)+1)>rand% && ~reject
                nano_sites(NC,t+1) = nano_sites(NC,t+1)+1;
            end
            if nano_sites(NC,t+1)==N_sites/2 && WTs(NC)==0 % write to ST if NC is half-swapped
                WTs(NC)=t+1;
            end
            if nano_sites(NC,t+1)==N_sites/4 && STs(NC)==0 % write to start if NC is 1/4 swapped
                STs(NC)=t+1;
                Starts(NC) = 1;
            end
            if nano_sites(NC,t+1)==3*N_sites/4 && Starts(NC)==1 % write to ST if NC is 3/4 swapped AND ST hasn't already been recorded already
                STs(NC)=(t+1-STs(NC))/2; % scale by 1/2
                Starts(NC) = 2; % avoid double counting
            end
            if mean(nano_sites(:,t+1))>tol*N_sites % tolerence condition for breaking the loop
                break
            end

        end

        WTs=WTs(WTs~=0); %keep WTs that are non-zero
        %STs=STs(STs~=0);
        nano_sites = nano_sites(:,1:j);
        Y = mean(nano_sites,1);

        % plot ensemble trajectory
        figure
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', figpos)
        plot(t_range(1:j),Y)
        labels = {labels{:},['time']};
        labels = {labels{:},['ions_ens']};
        data = {data{:},t_range(1:j)};
        data = {data{:},Y};
        ylabel('Average Incorporated Ions','FontSize',fs)
        xlabel('Time(a.u.)','FontSize',fs)
        title(strcat(type, '- Ensemble Trajectory'),'FontSize',fs)
        print(char(strcat(imgdir,'\ensemble_',filename)),'-dpng')
        
        % plot individual trajectories
        figure
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', figpos)

        for i=1:N_plots
            hold on 
            Y_i = nano_sites(i,:);
            plot(t_range(1:j),Y_i)
            ylabel('Incorporated Ions','FontSize',fs)
            xlabel('Time(a.u.)','FontSize',fs)
            title(strcat(type,'- Individual Trajectories'),'FontSize',fs)

        end

        % write important times of NCs to cell arrays
        labels = {labels{:},['ions_indv']};
        data = {data{:},nano_sites(1:N_plots,:)};

        labels = {labels{:},['waittimes']};
        data = {data{:},WTs'};
        
        labels = {labels{:},['switchtimes']};
        data = {data{:},STs'};

        % plot waiting time histogram
        hold off
        print(char(strcat(imgdir,'\individual_',filename)),'-dpng')
        figure
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', figpos)
        histogram(WTs,30)
        title(strcat(type,' - Wait Times Histogram'),'FontSize',fs)
        xlabel('Time(a.u.)','FontSize',fs)
        ylabel('# NCs in Bin','FontSize',fs)
        print(char(strcat(imgdir,'\waittimes_',filename)),'-dpng')
        
        % plot switching time histogram
        figure
        STs=STs(Starts==2); % only keep switched NCs
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', figpos)
        histogram(STs,30)
        title(strcat(type,' - Switch Times Histogram'),'FontSize',fs)
        xlabel('Time(a.u.)','FontSize',fs)
        ylabel('# NCs in Bin','FontSize',fs)
        print(char(strcat(imgdir,'\switchtimes_',filename)),'-dpng')
        
        % write timing statistics for ensemble to vector
        STmeds(end+1) = median(STs);
        WTmeds(end+1) = median(WTs);
        STerrs(end+1) = std(STs);
        WTerrs(end+1) = std(WTs);

        if isWrite %write MATDAT file
            write_MATDAT( fp,type,labels,data )
        end
    end
end

if isPlotK
    % plot timing statistics vs "concentration"
    % these additionaly plots only get made for the last model simulated above

    % Switching time median vs k
    figure
    hold on
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', figpos)
    errorbar(krange,STmeds,STerrs,'o')
    F = @(x,xdata) x(1)*xdata.^(-1)+x(2);
    x0 = [1*10^5  1];
    [x_ST,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,krange,STmeds);
    [x_WT,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,krange,WTmeds);

    xx = 1:.05:krange(end);
    [ b, yy ] = inversefit( krange,STmeds, xx,2 );
    %plot(xx,yy)

    plot(xx,F(x_ST,xx))
    xlabel('Concentration (k_0)')
    ylabel('Median Switching Time (a.u.)')
    print(char(dircat(imgdir,strcat(type,'-switching_times.png'))),'-dpng')

    % Waiting time median vs k
    figure
    hold on
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', figpos)
    errorbar(krange,WTmeds,WTerrs,'o')

    [ b, yy ] = inversefit( krange,WTmeds, xx,1 );
    %plot(xx,yy)
    plot(xx,F(x_WT,xx))
    xlabel('Concentration (k_0)')
    ylabel('Median Wait Time (a.u.)')
    print(char(dircat(imgdir,strcat(type,'-waiting_times.png'))),'-dpng')
    
    % write excel file of times vs concentration 
    xlfn = char(strcat(imgdir,filename,'03-23-2020.xlsx'));
    if isWriteXL
        write_NC_2_xlsx( xlfn, x_WT, x_ST, krange, WTerrs, STerrs, WTmeds, STmeds )
    end

end

% save ... vs swap plots
figure(Gfig)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', figpos)
plot([0,20],[0,0])
print(char(dircat(imgdir,'free_energy')),'-dpng')

figure(pffig)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', figpos)
print(char(dircat(imgdir,'probf')),'-dpng')
    
figure(psfig)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', figpos)
print(char(dircat(imgdir,'probs')),'-dpng')

toc % end timing






