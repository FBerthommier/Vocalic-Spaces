% simulation of Fant's model
% Makes the Figure 4b
% F. Berthommier 11/24/2021

function simuFANT

load simuFANT.mat

no_vibration=0; % Lossy
n=200; % high number of tubelets to avoid rounding effects
Lg=17.5; % length in cm

% Parameters of Table III
Lc=n*0.9; Lip=n-Lc; l=n*0.3;  
A=1; delta=Lc-l;
Psi1=[0.3*delta 2*A  -A 1.5]; 
Psi2=[5*pi/3 pi pi/3 pi/3]; 
Omega=[Lc/2 -1.5*A A/2 Lg]; 

L1=round(Omega(1)-(l/2));L2=Lc-l-L1;
P=[A.*ones(1,L1) Omega(2).*ones(1,l) A.*ones(1,L2) Omega(3).*ones(1,Lip)];
% Load of your own transmission line model with [length,area]
% [f1b, f2b]=TLM([ones(n,1).*Omega(4)./n,expif(P)'], no_vibration);

pastheta=pi/30; 
theta=0:pastheta:2*pi+pastheta;
rho=1;
indthetab=[1  11  16  21  31  41 46   51]; % Characteristic vowels

for k=1:length(theta)
disp(k)

Pval= Omega + rho.*Psi1.*cos(Psi2-theta(k)); 
L1=round(Pval(1)-(l/2));L2=Lc-l-L1;
% Makes a large tube
P=[A.*ones(1,L1) Pval(2).*ones(1,l) A.*ones(1,L2) Pval(3).*ones(1,Lip)];
% Load of your own transmission line model with [length,area]
% [f1(k), f2(k)]=TLM([ones(n,1).*Pval(4)./n,expif(P)'], no_vibration);

end;

figure(1)
axis([600 2300 225 675]); axis ij
xlabel('F2 (Hz)')
ylabel('F1 (Hz)')
hold on
plot(f2,f1, 'b'); plot(f2(indthetab),f1(indthetab), 'bO'); % C2

plot(f2b,f1b,'bO');
text(f2b-150,f1b-50,'(f1n,f2n)');

function vecout=expif(vecin)
vecout=(vecin < 1).*exp(vecin-1)+(vecin >=1).*vecin;

