%{
DelayDistributions_MLE.m - 
To run this file, first ensure that there is a data file in the directory.
Two formats are considered.
(1) If the input data is exactly observed, the data should have two columns. 
(2) If the first event date is interval censored, the file should have three columns.
The code will run the appropriate model based on the input data, and is currently set up for an underlying gamma
distribution, but this can easily be amended as required. 

Maximum likelihood estimates are returned for the distribution parameters.

Before running the code, the epidemic growth rate (r) and the truncation date (T) should be amended as required.

The following data files are provided:

IncubationPeriod.xlsx - File containing the incubation period data for cases with travel history related to Wuhan.
The first sheet contains the MATLAB input data. The second sheet contains the dates.
When event dates are equal, we have amended the MATLAB input by half a day,
since it is unlikely the events occurred simulataneously.

HospitalisationDelay_HongKong.xlsx and HospitalisationDelay_Singapore.xlsx- Files
containing the onset to hospitalisation delaydata for cases in Hong Kong (Singapore).
Similarly to above, first sheet is MATLAB input, second sheet contains dates,
which are editted if they coincide to produce the MATLAB input.
%}

%load input data: uncomment relevant line
% I = xlsread('IncubationPeriod.xlsx');
% I = xlsread('HospitalisationDelay_HongKong.xlsx');
I = xlsread('HospitalisationDelay_Singapore.xlsx');

T=max(I(:))+1; %set the truncation date, data set variable name 'I'
r=0.25; %epidemic growth rate
x0=[0.1,0.1]; %initial conditions for fminsearch
f=@(x)loglikelihood(I,x(1),x(2),r,T).^(-1); %loglikelihood function
out=fminsearch(f,x0); %output,
a=out(1);b=out(2); % a shape parameter, b scale parameter

%__________________________________________________________________________
function [logL] = loglikelihood(I,a,b,r,T)

for k=1:length(I)
    if size(I,2)==2
        %%%% both events observed exactly__________________________________________
        event1=I(k,1); %first event time
        event2=I(k,2); %second event time
        f1=@(x) gampdf(event2-x,a,b).*exp(r.*(x));
        f2=@(x) gamcdf(T-x,a,b).*exp(r.*(x));
        l(k)=f1(event1)./f2(event1);
        %%%%_______________________________________________________________________
    else
        %%%% first event interval censored_________________________________________
        event1L=I(k,1); %first event time - left censor
        event1R=I(k,2); %first event time - right censor
        event2=I(k,3); %second event time
        f1=@(x) gampdf(event2-x,a,b).*exp(r.*(x));
        f2=@(x) gamcdf(T-x,a,b).*exp(r.*(x));
        l(k)=integral(@(x) f1(x),event1L,event1R)./integral(@(x)f2(x),event1L,event1R);
        %%%%_______________________________________________________________________
    end
end
ll=log(l);
logL=sum(ll);

end