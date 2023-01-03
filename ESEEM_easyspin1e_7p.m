
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dummy=1;  % this variable will be used for different magfield offset
colour='rgbkcmy';
offset=[-3.5,-2.3,-1.2,0,1.2,2.3,3.5]; % in mT
%offset=[ 1];



%%% creation of random numbers %%%%
Nn=7;
realisation=10; %the final result is outcome of 10 different proton configuration
num_dd=Nn*(Nn-1)/2;
A=[1 -1 1 -1 1 -1 1]; B=[-1 1 -1 1 -1 1 -1]; % A random choice of +1 and -1 of my choice
plusmin=[A;B;A;B;A;B;A;B;A;B];
plusmin=plusmin';

for i=1:realisation
    XX1(:,i) = rand(Nn,1).*plusmin(:,i);
    XX2(:,i) = pi*rand(num_dd,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STEP=100;
for delB=1:length(offset)
  
Signal=zeros(1,STEP);
Signal_sz=zeros(1,STEP);

tic

for RR=1:1:realisation  %%%%% this loop is not in use




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temperature=5; %in K
kb=1.3806503e-23; %%% Boltzmann constant: m^2 kg / s^2 K
hplank=6.626068e-34;  %%% Planck constant m^2 kg / s  
beta=hplank/(kb*temperature)*10^6; % usec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining the system
Sys1=[1,  1/2,1/2,1/2,1/2,1/2,1/2,1/2];

%% using easyspin definig the operators
Identity=sop(Sys1,'e1');
x1 = sop(Sys1,'x1'); y1 = sop(Sys1,'y1'); z1 = sop(Sys1,'z1');  % electron operator

x2 = sop(Sys1,'x2'); y2 = sop(Sys1,'y2'); z2 = sop(Sys1,'z2');  % nucleus operators
x3 = sop(Sys1,'x3'); y3 = sop(Sys1,'y3'); z3 = sop(Sys1,'z3');
x4 = sop(Sys1,'x4'); y4 = sop(Sys1,'y4'); z4 = sop(Sys1,'z4');
x5 = sop(Sys1,'x5'); y5 = sop(Sys1,'y5'); z5 = sop(Sys1,'z5');
x6 = sop(Sys1,'x6'); y6 = sop(Sys1,'y6'); z6 = sop(Sys1,'z6');
x7 = sop(Sys1,'x7'); y7 = sop(Sys1,'y7'); z7 = sop(Sys1,'z7');
x8 = sop(Sys1,'x8'); y8 = sop(Sys1,'y8'); z8 = sop(Sys1,'z8');



%%
Ne=1;
Nn=7;
%%%% Spin operators
N=Ne+Nn;
Iz=zeros(length(x1),length(x1),N);
Ix=Iz; Iy=Iz;
%reassinging the operators as in future we have write the Hamiltonian in a loop
    Ix(:,:,1)=x1;Iy(:,:,1)=y1;Iz(:,:,1)=z1;     % electron operator
    
    Ix(:,:,2)=x2;Iy(:,:,2)=y2;Iz(:,:,2)=z2;     % nucleus operators
    Ix(:,:,3)=x3;Iy(:,:,3)=y3;Iz(:,:,3)=z3;
    Ix(:,:,4)=x4;Iy(:,:,4)=y4;Iz(:,:,4)=z4;
    Ix(:,:,5)=x5;Iy(:,:,5)=y5;Iz(:,:,5)=z5;
    Ix(:,:,6)=x6;Iy(:,:,6)=y6;Iz(:,:,6)=z6;
    Ix(:,:,7)=x7;Iy(:,:,7)=y7;Iz(:,:,7)=z7;
    Ix(:,:,8)=x8;Iy(:,:,8)=y8;Iz(:,:,8)=z8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Mag=[0:0.01:200]; 
shift=25; % to move the CT away from zero, first peak of Ho is also at 23.5 mT
detuning=offset(delB); % detuning from CT, mT % 1.2, 2.3,3.5mT is used in experiment
Magfield=shift+detuning; % in mT;


                
                we=(Magfield-shift)*28.024 ; % electron Larmour , MHz/mT
                wn=Magfield*0.04258; % proton larmour ,  MHz/mT, 
                

                S=1; D=-45000; E= 1*4500; % MHz, Anisoropy parameters, makes CT 9 GHz

                H0=  we*Iz(:,:,1); %  Electron zeeman
                H1=   D*(Iz(:,:,1)*Iz(:,:,1) - Identity*S*(S+1)/2) + E*(Ix(:,:,1)*Ix(:,:,1)-Iy(:,:,1)*Iy(:,:,1)); %electron anisotropy
             
                % e-n hyperfine intearction
                hp=4;  % in MHz, experimentally also we have seen similar value
               
                for i=1:Nn
                   %A(i)=(1+ (i-(Nn+1)/2)/24)*hp;  % this way I can make hyperfine linearly space from 7 to 9 MHz 
                   A(i)=hp+XX1(i, RR);
                end
                
                 %%%in general, psedusecular ternm is maller than secular tern, I kept it 1:2 ratio 
           
                Hhp=zeros(length(H0), length(H0));
                for i=1:Nn
                    Hhp= Hhp -  wn*Iz(:,:,i+1); 
                end 
                
                for i=1: Nn 
                Hhp= Hhp + A(i)*Iz(:,:,1)*Iz(:,:,i+1) + (A(i)/2)*Iz(:,:,1)*Ix(:,:,i+1);
                end
 
                %%%%%%%%%% dipolar hamiltonian  %%% %%%%%%%%%%%%%%%%
                Hj=zeros(length(H0), length(H0));
                test=1;
                dip=0.010;% dipolar coupling, MHz, 10 KHz proton-proton interaction is resonable 
                num_dd=Nn*(Nn-1)/2;
                for i=1:num_dd
                     dd(i)=-dip * (3*cos(XX2(i,RR))^2 -1);
                end
                
                for k=1:Nn-1
                    for l=k+1:Nn                
                    
                       Hj= Hj + dd(test)*(2*Iz(:,:,k+1)*Iz(:,:,l+1)) ;
                       Hj= Hj + dd(test)*(-Ix(:,:,k+1)*Ix(:,:,l+1)-Iy(:,:,k+1)*Iy(:,:,l+1));
                       
                       test=test+1;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Total Hamiltonian
               H=H0+H1+Hhp+Hj;
             
            
              %%%%%%%%%%%%%%%%%%%
              % diagonalisation for fun, not required for ESEEM
             % [V, Hd] = eig(H);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%

              %%%%%%% calculating equlibrium density matrix %%%%%%%%%%
              rho_eq=expm(-H*beta)./trace(expm(-H*beta));
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
              % Forming the time axis
                  % STEP=1000; % shoul be back to 1000, defined before
                   jump=.1; % increament of tau ,Step of 100ns, as maximum frequency we want to capture is less than 10 Mhz
                 % jump*STEP will decide the total time evoulution
                 
                    Sz1=0.5*(Ix(:,:,1)*Iy(:,:,1)+Iy(:,:,1)*Ix(:,:,1)); %The pulse propagator in +1 -1 subspace
                    Sz1D=Sz1;
              % some propagator that need NOT to be calculated for each tau
              
                    Pulse = expm(-1i*(pi/2)*Sz1D); % 90 degree pulse
                    Dens1 = Pulse*(rho_eq)*Pulse';  % forming the density matrix after 90 degree pulse   
                    Pulse2 = expm(-1i*pi*Sz1D);  % 180 degree pulse
                    freevol = expm(-2i*pi*jump*H); % the smallest unit of free evolution
              % Time evolution  % parallel computing
%               Signal=zeros(1,STEP);
%               Signal_sz=zeros(1,STEP);
              
                    parfor i=1:STEP
                        
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
           Signal_sz= Sig_sz +  Signal_sz; % finally this is going to be ESEEM decay incase of Clock transition
            
end

%%%%% signal processing and plotting only on Sz1D, and that is not our observable, so this part is not much relavant
%%%%%% But, I kept this part,  as this gives only osciilation, no expo decay, so dirrectly I can do FFFT
%%%%%% signal procecceing on the actual obserbavle, Sz, will be done by the other program



            %Signal=Signal-(sum(Signal)/length(Signal));
            figure(1)
            hold on
            plot(time,Signal,colour(dummy));
            
            xlabel('Tau(microsec)') 
            ylabel('Relative intensity(Arb Unit)') 
            
            
            figure(2)  %%% ESEEM signal in te paper
            hold on
            plot(time,Signal_sz,colour(dummy));
           
            xlabel('Tau(microsec)') 
            ylabel('Relative intensity(Arb Unit)') 
            

            
            
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% just an easy way to copy for Origin plot
dummy=dummy+1;
data1(:,1)=time; data1(:,dummy)=Signal;
data11(:,1)=time; data11(:,dummy)=Signal_sz;  %%% ESEEM signal in the paper


clearvars -except data1 data11 dummy delB colour wn A dd D E STEP jump hp dip offset Signal_sz Signal XX1   XX2   realisation 

toc
end
save 'ESEEM_results.mat'


%%%%%%%%%%%%%%%%%%%%%%%%  END of Program  %%%%%%%%%%%%%%%%%%%%%%%%%%%
