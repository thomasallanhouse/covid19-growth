%% Figure_S5
% Copyright (c) 2020 Lorenzo Pellis
% See LICENCE for licensing information
%
% Random infectivity profile and Lotka-Euler Equation relating growth rate and R0
% Figure S5 and table S1 in the electronic supplementary material of:
% 
% Pellis L, Scarabel F, Stage HB, Overton CE, Chappell LHK, Fearon E, Bennett E, 
% University of Manchester COVID-19 Modelling Group, Lythgoe KA, House TA and Hall I, 
% "Challenges in control of COVID-19: short doubling time and long delay to effect of interventions", 
% Philosophical Transactions of the Royal Society B (2021)

close all;
clearvars;

rng(7); % Fix random seed for reproducibility -- Fig_S5 and Table_S1 obtained with rng(7)
smax = 25; ds = 0.01; vs = 0:ds:smax; ls = length(vs);
np = 10000;
% profiles = {'TVI','TVIlate', 'SEIR','supspr'}; nprof = length(profiles);
profiles = {'SEIR','TVI', 'supspr', 'TVIlate'}; nprof = length(profiles);
vr = [ 0.1; 0.15; 0.2; 0.25; 0.3 ]; lr = length(vr); % Vector of values for the growth rate
doubt = log(2)./vr; % Doubling time
R0 = NaN(lr,nprof);
gentime = NaN(lr,nprof);
Ip = 7; % Baseline duration of infectious period
Dp = 3; % Duration of prodromal period: manually change values 0,1,2,3
Ps = 1; % Proportion of symptomatic infections. Try 1 and 0.5
Fa = 0.5; % Fraction of transmission if asymptomatic
ninfprofplot = 10; % Plot only some infectivity profiles
iinfprofplot = round(linspace(1,np,ninfprofplot)); 
iinfprofplot(end) = round( iinfprofplot(end)-(iinfprofplot(end)-iinfprofplot(end-1))/20 ); % Looks better
cplot = jet(ninfprofplot); % Colours for the different infectivity profiles
set(groot,'defaultAxesColorOrder',cplot);
for iprof = 1:nprof
    ptab = zeros(np,5); % Duration of: Latent (E), Prodromal (P), Infectious (I) period, total infectiousness (R), symptomatic indicator (S)
    omega = zeros(np,ls); % Generation time distribution (scaled down by relative transmission if asymptomatic)
    omega_sympt = zeros(np,ls); % Generation time distribution
    symptoms = NaN(np,ls);
    switch profiles{iprof}
        case 'SEIR'
            Em = 4.85-Dp; Ea = 1; Eb = Em/Ea; % Reconstruct parameters a and b from mean
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration latent period
            Im = 3; Ia = 1; Ib = Im/Ia; % Reconstruct parameters a and b from mean
            ptab(:,3) = random('Gamma',Ia,Ib,[np,1]); % Duration infectious period
            ptab(:,4) = 1; % No further variation in total infectiousness
        case 'TVI' % Gamma-shaped infectivity with mean 2 and std 1.5
            Em = 4.85-Dp; Es = 2.79; Ea = (Em/Es)^2; Eb = Es^2/Em; % Reconstruct parameters a and b from mean and standard deviation 
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration of latent period
            ptab(:,3) = 7;
            Om = 2; Os = 1.5; Oa = (Om/Os)^2; Ob = Os^2/Om;  % Reconstruct parameters a and b from mean and standard deviation 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            ptab(:,4) = 1; % No variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        case 'supspr' % Gamma-shaped infectivity with mean 2 and std 1.5, but variability in total infectivity
            Em = 4.85-Dp; Es = 2.79; Ea = (Em/Es)^2; Eb = Es^2/Em; % Reconstruct parameters a and b from mean and standard deviation 
            ptab(:,1) = sort(sort(random('Gamma',Ea,Eb,[np,1]))); % Duration of latent period
            ptab(:,3) = Ip;
            Om = 2; Os = 1.5; Oa = (Om/Os)^2; Ob = Os^2/Om; % Reconstruct parameters a and b from mean and standard deviation 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            k = 0.25; Km = 1; Ks = 1/sqrt(k); Ka = (Km/Ks)^2; Kb = Ks^2/Km; 
            ptab(:,4) = random('Gamma',Ka,Kb,[np,1]); % Variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        case 'TVIlate' % Gamma-shaped infectivity with mean 3 and std 1.5
            Em = 4.85-Dp; Es = 2.79; Ea = (Em/Es)^2; Eb = Es^2/Em; % Reconstruct parameters a and b from mean and standard deviation 
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration of latent period
            ptab(:,3) = 7;
            Om = 3; Os = 2; Oa = (Om/Os)^2; Ob = Os^2/Om; % Reconstruct parameters a and b from mean and standard deviation 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            ptab(:,4) = 1; % No variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        otherwise
            error('No good profile specified');
    end
    ptab(:,5) = (rand(np,1) < Ps); % 1 if symptomatic, 0 if asymptomatic

    
    for ip = 1:np
        ptab(ip,2) = Dp; % Duration prodromal period
        if strcmp(profiles{iprof},'SEIR')
            vO = ones(size(0:ds:ptab(ip,3)))/ptab(ip,3); lO = length(vO);
        end
        iE = ceil(ptab(ip,1)/ds); % index at which the infectious period starts
        iD = ceil((ptab(ip,1)+ptab(ip,2))/ds); % index at which symptoms appear (if symptomatic)
        iR = ceil((ptab(ip,1)+ptab(ip,3))/ds); % index at which symptoms stop (if symptomatic)
        if (iE+lO) <= ls
            omega(ip,(iE+1):(iE+lO)) = vO;
        else
            omega(ip,(iE+1):ls) = vO(1:(ls-iE));
        end
        if ptab(ip,5) == 1
            omega_sympt(ip,:) = ptab(ip,4) * omega(ip,:);
        else
            omega_sympt(ip,:) = Fa * ptab(ip,4) * omega(ip,:);
        end
        symptoms(ip,(iD+1):min(iR,ls)) = ptab(ip,5);
    end
    figure(iprof); clf; hold on;
    plot(vs,omega_sympt(iinfprofplot,:),'Linewidth',1);
    avomega = mean(omega,1);
    avomega_sympt = mean(omega_sympt,1);
    gt = ds*trapz(vs.*avomega_sympt);
    symptoms(isnan(symptoms)) = 0;
    avsymptoms = mean(symptoms,1);
    xlabel('Time (days)');
    ylabel('Density');
    plot(vs,avomega_sympt,'k','Linewidth',3);


    for ir = 1:lr
        r = vr(ir);
        R0(ir,iprof) = 1 / (ds*trapz(avomega.*exp(-r*vs)));
        gentime(ir,iprof) = gt;
    end

end
output = [vr, doubt, R0, gentime];
% Figure_S5 manually saved from figure(2) (i.e. TVI) with following options: 
% rng(7), Dp = 1 (i.e. 1-day prodromal period), Ps = 1, Fa = 0.5

% Table_S1 values obtained from variable "output" setting rng(7), and manually
% changing Dp to values 0, 1, 2 and 3, and copying in Table_S1.xlsx
% the correct columns (i.e. 1, 2, 4, 6 for Table_S1A and 1, 2, 3, 5 for
% Table_S1B). The last columns (gentime) are independent of growth rates
% and are copied at the bottom of the table (columns 8 and 10 for Table_S1A
% and columns 7 and 9 for Table_S1B).