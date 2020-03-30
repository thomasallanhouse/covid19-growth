function Covid19growthNBfit()

% Comment out which line is not needed
country = 'UK'; inccol = 15;
% country = 'IT'; inccol = 8;
addci = 1;

dat = importdata('who_data.csv');
data = diff([zeros(1,15);dat.data]);

startdaynum = datenum('10-Jan-2020'); % absolute plot start
startdata = datenum('09-Jan-2020');
[ndd,~] = size(data); % number of days of total data
enddata = startdata + ndd - 1;
switch country
    case 'UK'
        datafitstart = datenum('01-Mar-2020');
        datafitend = datenum('19-Mar-2020');
        predstart = datafitstart;
        predend = datafitend + 8;
        datateststart = datafitend + 1;
        datatestend = datafitend + 6;
        countryintitle = 'the UK';
    case 'IT'
        datafitstart = datenum('21-Feb-2020');
        datafitend = datenum('28-Feb-2020');
        predstart = datafitstart;
        predend = datafitend + 8;
        datateststart = datafitend + 1;
        datatestend = datafitend + 7;
        countryintitle = 'Italy';
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

ydata = data(ixdata,inccol);
ytest = data(ixtest,inccol);


% This does the MLE fit
f=@(x)nll(x(1),x(2),x(3),xdata,ydata);
x0=[2,0.1,10];
xhat=fminsearch(f,x0);

% This works out a 95% CI on r
h=1e-3;
sig = h/sqrt(f(xhat + [0.0, h, 0.0]) - 2.0*f(xhat) + f(xhat - [0.0, h, 0.0]));

% This calculates a prediction region from the model and a mean
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


if addci
    legentry = {'Predicted','95% CI','Data (fit)','Data (test)'};
else
    legentry = {'Predicted','Data (fit)','Data (test)'};
end

% Plot fitting daily cases
figure(5); clf; hold on;
plot( xpred, ym, 'r-', 'Linewidth', 2 );
if addci
    plot( xpred, yl, 'r:', 'Linewidth', 2 );
end
plot( xdata, ydata, 'ko', 'Linewidth', 2 );
plot( xtest, ytest, 'kx', 'Linewidth', 2, 'MarkerSize', 8 );
if addci
    plot( xpred, yu, 'r:', 'Linewidth', 2 );
end
limx = [min(xpred) max(xpred)];
tickxfull = [xdata;xtest];
set(gca,'Xlim',limx,'XTick',tickxfull(mod(tickxfull,7)==2));
set(gca,'Ylim',[0,2500]);
limy = get(gca,'Ylim');
legend(legentry,'Location','NorthWest')

title({
    ['Daily cases in ',countryintitle],
    ['growth rate = ',num2str(xhat(2),'%1.2f'),'(',num2str(xhat(2)-2*sig,'%1.2f'),',',num2str(xhat(2)+2*sig,'%1.2f'),') day^{-1}; doubling time = ',...
    num2str(log(2)/(xhat(2)),'%1.2f'),'(',num2str(log(2)/(xhat(2)+2*sig),'%1.2f'),',',num2str(log(2)/(xhat(2)-2*sig),'%1.2f'),') days']
    })
xlabel('Time (days) - labelled days are Sundays')
ylabel('Numbers')
datetick('x',19,'keeplimits','keepticks')

end

function N=nll(a,r,y0,xdata,ydata)

N = 0;
for i=1:length(ydata)
    N = N - log(nbinpdf(ydata(i),y0*exp(r*(xdata(i)-xdata(1)))/(a-1),1/a));
end

end

