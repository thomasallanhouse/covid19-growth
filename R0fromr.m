%%% Random infectivity profile and Lotka Euler Equation
close all;
clearvars;

rng(7);
smax = 25; ds = 0.01; vs = 0:ds:smax; ls = length(vs);
np = 10000;
profiles = {'SEIR','TVI', 'supspr', 'TVIlate'}; nprof = length(profiles);
vr = [ 0.1; 0.15; 0.2; 0.25; 0.3 ]; lr = length(vr);
doubt = log(2)./vr;
R0 = NaN(lr,nprof);
Ip = 7; % Baseline duration of infectious period
Dp = 2; % Duration of prodromal period: manually change values 0,1,2,3
Ps = 1; % Proportion of symptomatic infections. Try 1 and 0.5
Fa = 0.5; % Fraction of transmission if anymptomatic
ninfprofplot = 10;
iinfprofplot = round(linspace(1,np,ninfprofplot));
cplot = jet(ninfprofplot);
set(groot,'defaultAxesColorOrder',cplot);
for iprof = 1:nprof
    ptab = zeros(np,5); % Duration of: Latent (E), Prodromal (P), Infectious (I) period, total infectiousness (R), symptomatic? (S)
    omega = zeros(np,ls);
    omega_sympt = zeros(np,ls);
    symptoms = NaN(np,ls);
    switch profiles{iprof}
        case 'SEIR'
            Em = 4.84-Dp; Ea = 1; Eb = Em/Ea;
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration latent period
            Im = 3; Ia = 1; Ib = Im/Ia; %%%% Check mean inf period
            ptab(:,3) = random('Gamma',Ia,Ib,[np,1]); % Duration infectious period
            ptab(:,4) = 1; % No further variation in total infectiousness
        case 'shortSI' %%% To check!
            Em = 2; Es = 1.5; Ea = (Em/Es)^2; Eb = Es^2/Em; 
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration of latent period
            ptab(:,3) = Ip;
            Om = 2; Os = 1.5; Oa = (Om/Os)^2; Ob = Os^2/Om; 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            ptab(:,4) = 1; % No variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        case 'TVI'
            Em = 4.84-Dp; Es = 2.79; Ea = (Em/Es)^2; Eb = Es^2/Em; 
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration of latent period
            ptab(:,3) = 7;
            Om = 2; Os = 1.5; Oa = (Om/Os)^2; Ob = Os^2/Om; 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            ptab(:,4) = 1; % No variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        case 'supspr'
            Em = 4.84-Dp; Es = 2.79; Ea = (Em/Es)^2; Eb = Es^2/Em; 
            ptab(:,1) = sort(sort(random('Gamma',Ea,Eb,[np,1]))); % Duration of latent period
            ptab(:,3) = Ip;
            Om = 2; Os = 1.5; Oa = (Om/Os)^2; Ob = Os^2/Om; 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            k = 0.25; Km = 1; Ks = 1/sqrt(k); Ka = (Km/Ks)^2; Kb = Ks^2/Km; 
            ptab(:,4) = random('Gamma',Ka,Kb,[np,1]); % Variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        case 'TVIlate'
            Em = 4.84-Dp; Es = 2.79; Ea = (Em/Es)^2; Eb = Es^2/Em; 
            ptab(:,1) = sort(random('Gamma',Ea,Eb,[np,1])); % Duration of latent period
            ptab(:,3) = 7;
            Om = 3; Os = 2; Oa = (Om/Os)^2; Ob = Os^2/Om; 
            vO = pdf('Gamma',0:ds:ptab(1,3),Oa,Ob); lO = length(vO);
            ptab(:,4) = 1; % No variation in total infectiousness
            appintO = (ds*trapz(vO)); vO = vO / appintO; intO = ds*trapz(vO);
            disp(['Original truncation = ',num2str(1-appintO),'; Corrected truncation = ',num2str(1-intO)])
        otherwise
            error('No good profile specified');
    end
    ptab(:,5) = (rand(np,1) < Ps); % 1 if you are symptomatic, 0 if asymptomatic

    
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
