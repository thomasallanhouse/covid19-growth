function Figure_3()
% Copyright (c) 2020 Lorenzo Pellis
% See LICENCE for licensing information
%
% Growth rate fit to confirmed cases in UK and Italy 
% using MLE with negative binomial likelihood
% Figure 3 in the main text of:
% 
% Pellis L, Scarabel F, Stage HB, Overton CE, Chappell LHK, Fearon E, Bennett E, 
% University of Manchester COVID-19 Modelling Group, Lythgoe KA, House TA and Hall I, 
% "Challenges in control of COVID-19: short doubling time and long delay to effect of interventions", 
% Philosophical Transactions of the Royal Society B (2021)

    savefigs = 0; % save fig files if turned on
    addpi = 1; % computes prediction intervals if turned on

    dat = importdata('WHO_data_18.csv');
    data = diff([zeros(1,18);dat.data]);

    % select country (UK or IT)
    country = 'UK';
%     country = 'IT';

    startdata = datenum('09-Jan-2020');

    switch country
        case 'UK'
            inccol = 18; 
            datafitstart = datenum('23-Feb-2020');
            datafitend = datenum('12-Mar-2020');
            predstart = datafitstart;
            predend = datafitend + 15;
            datateststart = datafitend + 1;
            datatestend = datenum('30-Mar-2020');
            countryintitle = 'the UK';
            x0 = [100,0.1,10];
            ymax = 4000; 
        case 'IT'
            inccol = 8; 
            datafitstart = datenum('23-Feb-2020');
            datafitend = datenum('4-Mar-2020');
            predstart = datafitstart;
            predend = datafitend + 15;
            datateststart = datafitend + 1;
            datatestend = datenum('30-Mar-2020'); 
            countryintitle = 'Italy';
            x0 = [100,0.1,10];
            ymax = 8000;
        otherwise
            error('No valid country selected');
    end
    xdata = (datafitstart:datafitend)';
    dt = 0.1;
    xpred = (predstart:dt:predend)';
    xtest = (datateststart:datatestend)';
    % relative indexing
    ixdatastart = datafitstart - startdata + 1;
    ixdataend = datafitend - startdata + 1;
    ixdata = (ixdatastart:ixdataend)';
    ixteststart = datateststart - startdata + 1;
    ixtestend = datatestend - startdata + 1;
    ixtest = (ixteststart:ixtestend)';

    ydata = data(ixdata,inccol); % daily confirmed cases
    ytest = data(ixtest,inccol);

    % To perform the MLE fit
    f = @(x)nll(x(1),x(2),x(3),xdata,ydata); % Negative log-likelihood
    xhat = fminsearch(f,x0);

    % To produce a 95% confidence interval around the MLE of r (assuming asymptotic normality)
    h = 1e-3; % Small step for finite difference of second derivative (next line)
    sig = h/sqrt(f(xhat + [0.0, h, 0.0]) - 2.0*f(xhat) + f(xhat - [0.0, h, 0.0])); % Standard deviation sigma

    % To calculate a prediction region from the model and a mean
    yl = 0*xpred;
    ym = 0*xpred;
    yu = 0*xpred;
    a = xhat(1);
    r = xhat(2);
    y0 = xhat(3);
    for i=1:length(xpred)
        yl(i) = nbininv(0.025,y0*exp(r*(xpred(i)-xdata(1)))/(a-1),1/a);
        ym(i) = y0*exp(r*(xpred(i)-xdata(1)));
        yu(i) = nbininv(0.975,y0*exp(r*(xpred(i)-xdata(1)))/(a-1),1/a);    
    end

    if addpi
        legentry = {'Predicted','95% PI','Data (fit)','Data (test)'};
    else
        legentry = {'Predicted','Data (fit)','Data (test)'};
    end

    % Plot fitting daily cases
    figure(5); clf; hold on
    plot( xpred, ym, 'r-', 'Linewidth', 2 );

    if addpi
        plot( xpred, yl, 'r:', 'Linewidth', 2 );
    end
    plot( xdata, ydata, 'ko', 'Linewidth', 2, 'MarkerSize', 7 );
    plot( xtest, ytest, 'kx', 'Linewidth', 2, 'MarkerSize', 9 );
    if addpi
        plot( xpred, yu, 'r:', 'Linewidth', 2 );
    end
    limx = [min(xdata) max(xtest)];
    tickxfull = [xdata;xtest];
    set(gca,'Xlim',limx,'XTick',tickxfull(mod(tickxfull,7)==2));%,'XMinorTick','on');
    set(gca,'Ylim',[0,ymax]);
    limy = get(gca,'Ylim');
    legend(legentry,'Location','NorthWest','Fontsize',11)

    title({
        ['Daily cases in ',countryintitle],...
        ['growth rate = ',num2str(xhat(2),'%1.2f'),'(',num2str(xhat(2)-1.96*sig,'%1.2f'),',',num2str(xhat(2)+1.96*sig,'%1.2f'),') day^{-1}; doubling time = ',...
        num2str(log(2)/(xhat(2)),'%1.2f'),'(',num2str(log(2)/(xhat(2)+1.96*sig),'%1.2f'),',',num2str(log(2)/(xhat(2)-1.96*sig),'%1.2f'),') days']
        })
    xlabel('Time (days) - labelled days are Sundays')
    ylabel('Daily case counts')
    datetick('x',19,'keeplimits','keepticks')
    set(gca,'Fontsize',11);

    % To save figures
    if savefigs
        figname = [country,'_daily_cases_NBfit.eps'];
        exportgraphics(gca,figname);
    end

end

function N = nll(a,r,y0,xdata,ydata)
% Negative negative binomial log-likelihood

    N = 0;
    for i=1:length(ydata)
        N = N - log(nbinpdf(ydata(i),y0*exp(r*(xdata(i)-xdata(1)))/(a-1),1/a));
    end

end

