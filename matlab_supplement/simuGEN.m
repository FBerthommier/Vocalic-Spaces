% simulation of the generic function
% Makes Figure 1 
% F. Berthommier 11/24/2021

function simuGEN

load simuGEN.mat

no_vibration=1; % Lossless
n=100; % Number of tubelets
pas=1/n; 
x=pas/2:pas:1-(pas/2); 
v1=cos(pi*x); 
v2=cos(3*pi*x); 
L=17.5; % Length in cm
rho=1; 

% Load of your own transmission line model with [length,area]
% [f1b, f2b]=TLM([ones(n,1).*L./n,ones(n,1)], no_vibration); % neutral tube
disp(f1b); disp(f2b);

pastheta=pi/60; 
theta=0:pastheta:2*pi+pastheta;
indthetab=[ 1 21 31 41 61 81 91  101]; % Characteristic vowels

for k=1:length(theta)
disp(k)
% Three-phase mixing function
P= 1+2*rho*cos(theta(k)).*v1+(4/3)*rho*sin(pi/3)*sin(theta(k)).*v2;
   % Load of your own transmission line model with [length,area]
   % [f1(k), f2(k)]=TLM([ones(n,1).*L./n,expif(P)'], no_vibration);
end;

figure(1)
axis([600 2400 150 1000]); axis ij
xlabel('F2 (Hz)')
ylabel('F1 (Hz)')
hold on
plot(f2,f1, 'b'); plot(f2(indthetab),f1(indthetab), 'bO'); % C2

plot(f2b,f1b,'bO');
text(f2b-150,f1b-50,'(f1n,f2n)');

function vecout=expif(vecin)
vecout=(vecin < 1).*exp(vecin-1)+(vecin >=1).*vecin;

    






