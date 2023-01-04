clear all
clc
dummy=1;

offset=[0,1,2,3,5,10,15,20,25,30,50,100,200]; % in mT
offset=[20];



for delB=1:length(offset)
%%
tic
realisation=1;  %%%% this loop is actually not in use
for RR=1:1:realisation  %%%%% this loop is not in use




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temperature=5; %in K
kb=1.3806503e-23; %%% Boltzmann constant: m^2 kg / s^2 K
hplank=6.626068e-34;  %%% Planck constant m^2 kg / s  
beta=hplank/(kb*temperature)*10^6; % usec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining the system
Sys1=[1,  1/2];

%% using easyspin definig the operators
Identity=sop(Sys1,'e1');
x1 = sop(Sys1,'x1'); y1 = sop(Sys1,'y1'); z1 = sop(Sys1,'z1');  % electron operator

x2 = sop(Sys1,'x2'); y2 = sop(Sys1,'y2'); z2 = sop(Sys1,'z2');  % nucleus operators



%%
Ne=1;
Nn=1;
%%%% Spin operators
N=Ne+Nn;
Iz=zeros(length(x1),length(x1),N);
Ix=Iz; Iy=Iz;
%reassinging the operators as in future we have write the Hamiltonian in a loop
    Ix(:,:,1)=x1;Iy(:,:,1)=y1;Iz(:,:,1)=z1;     % electron operator
    
    Ix(:,:,2)=x2;Iy(:,:,2)=y2;Iz(:,:,2)=z2;     % nucleus operators
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Mag=[0:0.01:200]; 
shift=23.5; % to move the CT away from zero, first peak of Ho is also at 23.5 mT
detuning=offset(delB); % detuning from CT, mT % 1.2, 2.3,3.5mT is used in experiment
Magfield=shift+detuning; % in mT;


                
                we=(Magfield-shift)*28.024 ; % electron Larmour , MHz/mT
                wn=Magfield*0.04258; % proton larmour ,  MHz/mT, 
                

                S=1; D=-45000; E= 4500; % MHz, Anisoropy parameters, makes CT 9 GHz

                H0=  we*Iz(:,:,1); %  Electron zeeman
                H1=   D*(Iz(:,:,1)*Iz(:,:,1) - Identity*S*(S+1)/2) + E*(Ix(:,:,1)*Ix(:,:,1)-Iy(:,:,1)*Iy(:,:,1)); %electron anisotropy
             
                % e-n hyperfine intearction
                hp=1;  % in MHz, experimentally also we have seen similar value
             
               
                          
                Hhp=zeros(length(H0), length(H0));
           
                Hhp= Hhp +  wn*Iz(:,:,2)+ hp*Iz(:,:,1)*Iz(:,:,2) + (hp/2)*Iz(:,:,1)*Ix(:,:,2); 
              
 
               
                       
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Total Hamiltonian
               H=H0+H1+Hhp;
             
            
              %%%%%%%%%%%%%%%%%%%
              % diagonalisation for fun, not required for ESEEM
              [V, Hd] = eig(H);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%

              %%%%%%% calculating equlibrium density matrix %%%%%%%%%%
              rho_eq=expm(-H*beta)./trace(expm(-H*beta));
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
              % Forming the time axis
                   STEP=5000; % 1000
                   jump=0.02; % increament of tau ,Step of 100ns, as maximum frequency we want to capture is less than 10 Mhz
                 % jump*STEP will decide the total time evoulution
                 
                    Sz1=0.5*(Ix(:,:,1)*Iy(:,:,1)+Iy(:,:,1)*Ix(:,:,1)); %The pulse propagator in +1 -1 subspace
                    Sz1D=Sz1;
              % some propagator that need NOT to be calculated for each tau
              
                    Pulse = expm(-1i*(pi/2)*Sz1D); % 90 degree pulse
                    Dens1 = Pulse*(rho_eq)*Pulse';  % forming the density matrix after 90 degree pulse   
                    Pulse2 = expm(-1i*pi*Sz1D);  % 180 degree pulse
                    freevol = expm(-2i*pi*jump*H); % the smallest unit of free evolution
              % Time evolution  % parallel computing
              Signal=zeros(1,STEP);
              Signal_sz=zeros(1,STEP);
              
                    for i=1:STEP
                        
                        tau=jump+(i-1)*jump;
                        time(i)=tau; % forming the time axis
                       
                        % Propagation operators

                        %Dens1 = Pulse*(H)*Pulse';       % density matric after 90 pulse, but it is already calculated
                        
                        TauEvol = freevol^i;             % free evoution propagator for each tau
                        
                        Dens2=TauEvol*Dens1*TauEvol';    % density matrix after first free evolution              
                        Dens3=Pulse2*Dens2*Pulse2';      % density matric after 180 pulse
                        Dens4=TauEvol*Dens3*TauEvol';    % density matrix after second free evolution
                        
                        Sig(i)=real(trace(Dens4*(Sz1D)));  % if we observe SxSy + SySx, equvalent to observe sigma_y in [+1 -1] manifold
                        Sig_sz(i)=abs(trace(Dens4*(z1)));  % if we observe Sz, equvalent to observe sigma_x in [+1 -1] manifold
                    end
                    
           Signal= Sig +  Signal;    
           Signal_sz= Sig_sz +  Signal_sz; % finally this is going to be our observable
            
end

%%%%% signal processing and plotting only on Sz1D, and that is not our observable, so this part is not much relavant
%%%%%% But, I kept this part,  as this gives only osciilation, no expo decay, so dirrectly I can do FFFT
%%%%%% signal procecceing on the actual obserbavle, Sz, will be done by the other program



            Signal_sz=Signal_sz-(sum(Signal_sz)/length(Signal_sz));
%             figure(200)
%             hold on
%             plot(time,Signal);
%             
%             xlabel('Tau(microsec)') 
%             ylabel('Relative intensity(Arb Unit)') 
%             
            
            figure(100)
            hold on
            plot(time,Signal_sz);
           
            xlabel('Tau(microsec)') 
            ylabel('Relative intensity(Arb Unit)') 
            

            %%%%% forming the exponentoial fucntion for apodisation
            EPO=100/STEP; %we need to tune this number according to choice step, provided step is 1 ns
            for i=1:length(time)
            fun(i)=exp(-i*EPO);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Signal_apo=Signal_sz.*exp(-time*EPO); % Apodisaiton
           %figure(200)
           %plot(time,Signal_apo);
            
            
            %%% FFT of the signal
            npoints=STEP*10+1;
            S=abs (time(1)-time(2));
            spectra1=abs(fftshift(fft(Signal_apo, npoints))); %%%magnitude 
            freq =(1/S).*(-npoints/2:npoints/2-1)./npoints; % creating frequency axis
            figure(dummy)
            hold on
            plot ( freq, abs(spectra1));
            xlim([0 50])
            %ylim([0 3.5])
            xlabel('Frequency(MHz)') 
            ylabel('Relative Intensity(Arb Unit)') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% just an easy way to copy for Origin plot
dummy=dummy+1;
% data1(:,1)=time; data1(:,dummy)=Signal;
data11(:,1)=time; data11(:,dummy)=Signal_sz;

data2(:,1)=freq((npoints+1)/2:1:npoints); data2(:,dummy)=spectra1((npoints+1)/2:1:npoints);

clearvars -except data1 data11 data2 dummy delB colour wn A dd D E STEP jump hp dip offset Signal_sz Signal

toc
end

tt=1;


% 
% END=15;  % in us
% [val, index]= min( abs( data11(:,1) - END));
% for i=2
%   x=data11(1:index,1);
%   y=data11(1:index,i);
%   %plot(x,y)
% 
% myfittype = fittype('a*exp(-(x/T)^b)+c', 'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c', 'T'});
% [f] = fit(x,y,myfittype,'StartPoint', [.024, 3, .005, 8.6]);
% disp (f);
% figure (1000)
% plot(f,x,y,'*')
% end
% npoints=10000;
% S= x(2)-x(1);   
% 
% expo=f(x); 
% y=y-expo;
% y=y-(sum(y(end-4:end))/5);
% 
% EPO=.00; 
% for i=1:length(x)
% fun(i)=exp(-x(i)*EPO);
% end
% y=y.*fun; % Apodisaiton
% 
% figure (1001)
% plot (x,y)
% 
% 
% spectra_final=abs(fftshift(fft(y, npoints))); %%%magnitude 
% freq =(1/S).*(-npoints/2:npoints/2-1)./npoints; % creating frequency axis
% figure(1002)
% hold on
% plot ( freq, abs(spectra_final));
% xlim([0 5])
% %ylim([0 3.5])
% xlabel('Frequency(MHz)') 
% ylabel('Relative Intensity(Arb Unit)') 


