clear all
last=100;  % how much of the time trace will be considered for the fit to streached exponential, 100 us for example

load 'ESEEM_results.mat' %%%% the saved worksapce of the "ESEEM_easyspin1e_7p.m" 
[val, index]= min( abs( data11(:,1) - last));

%%%%% fitting to streched ecponential %%%%%%
for i=4 % at a time I took one of the eseem traces, i=2 is -3.5mT data, 3=-2.3mT and so on,  i goes from 2 to 8
  x=data11(1:index,1);
  y=data11(1:index,i);
  %plot(x,y)

myfittype = fittype('a*exp(-(x/T)^b)+c', 'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c', 'T'});
[f] = fit(x,y,myfittype,'StartPoint', [.24, 3, .05, 8.6]);
disp (f);
figure (1000)
plot(f, x,y,'-')
xlabel('Time (microsec)')
ylabel('ESEEM signal')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npoints=10240; % 10240 points
S= x(2)-x(1); % time spacing  

expo=f(x); 
y=y-expo;  %%%%%% subtracting exponential
y=y-(sum(y(end-4:end))/5);  %%%% making DC zero, approximately, just by taking an average of last 5 points

%%%%%% apodisation with an exponential function

EPO=0.1; % this value kind of make a flat deacy of the osciilation by 20 us

for i=1:length(x)
fun(i)=exp(-x(i)*EPO);
end

Ap=1;win=0; % for a comparison of the use of apodiation versus a proper windowing function in the FFT

if Ap==1
y=y.*fun'; % Apodisaiton
end

%%%% windowing fucntion %%%
kk=length(x);
%ww = chebwin(kk, 40);
ww = hamming(kk,'symmetric');
%ww = hamming(kk,'periodic');
%ww = hann(kk,'symmetric');
figure (2000); plot (ww);

if win==1
    y=y .*ww ;
end
%%%%%%%%%%%%%%%%%%%
figure (1001)
hold on

if Ap==1 && win==0 
    plot ( x, y,'r');
elseif Ap==0 && win==1
    plot ( x,y,'b');
elseif Ap==0 && win==0
    plot ( x,y, 'k');
end

xlabel('Time (microsec)')
ylabel('ESEEM - Fitting Func')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FFT %%%%%%%%%%%%%%%%%%%
spectra_final=abs(fftshift(fft(y, npoints))); %%%magnitude 
freq =(1/S).*(-npoints/2:npoints/2-1)./npoints; % creating frequency axis
figure(1002)
hold on
if Ap==1 && win==0 
    plot ( freq, abs(spectra_final),'r');
elseif Ap==0 && win==1
    plot ( freq, 10*abs(spectra_final),'b');
elseif Ap==0 && win==0
    plot ( freq, abs(spectra_final),'k');
end

xlim([.1 5])

xlabel('Frequency(MHz)') 
ylabel('Relative Intensity(Arb Unit)') 
spectra_final=spectra_final';
%%%%% For ploting, it is desirable to plot the frequency axis start from 0, that is why it starts from (npoints+1)/2
data22(:,1)=freq((npoints+1)/2:1:npoints); data22(:,2)=spectra_final((npoints+1)/2:1:npoints); 

