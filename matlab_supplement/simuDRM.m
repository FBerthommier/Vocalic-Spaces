% simulation of the DRM
% Makes the Figure 3
% F. Berthommier 11/24/2021

function simuDRM

load simuDRM.mat

no_vibration=1; % Lossless
n=120; % adapted to the DRM cut
pas=1/n; 
x=pas/2:pas:1-(pas/2); 
v1=cos(pi*x); 
v2=cos(3*pi*x); 
L=17.5; % Length in cm
rho=1;
% Parameters of Table II
if 0
Omega=[1 1 1 1]; 
Psi2=[atan(2*sin(pi/3)/3) atan(-4*sin(pi/3)/3) atan(-4*sin(pi/3)/3) atan(2*sin(pi/3)/3)]; 
Psi1=[2 1 -1 -2]./cos(Psi2);
else % Let's try the transformation function
    ijk=[3 -1 1
        0.5 0 2.5
        2-[0.5 0 2.5]
        2-[3 -1 1]];
    [Omega, Psi1, Psi2]=ijk2psi(ijk); 
end;

% Load of your own transmission line model with [length,area]
% [f1b, f2b]=TLM([ones(n,1).*L./n,ones(n,1)], no_vibration); % neutral tube

pastheta=pi/60; 
theta=0:pastheta:2*pi+pastheta;
indthetab=[1 21 31 41 61 81 91  101]; % Characteristic vowels

for k=1:length(theta)
disp(k)
% Coordination function
Pval = Omega + rho.*Psi1.*cos(Psi2-theta(k));
% Makes a large ntube
P=[Pval(1).*ones(1,n/6) Pval(2).*ones(1,n/3) Pval(3).*ones(1,n/3) Pval(4).*ones(1,n/6)];

a1(k)=2*sum(expif(P).*v1)./n;
a2(k)=2*sum(expif(P).*v2)./n;

a10(k)=2*sum(P.*v1)./n;
a20(k)=2*sum(P.*v2)./n;

    f1est(k)=f1b - f1b*((a1(k)/2) + (a2(k)/2)^2);  
    f2est(k)=f2b - f2b*a2(k)/2;
    f1est0(k)=f1b - f1b*a1(k)/2;
    f2est0(k)=f2b - f2b*a2(k)/2;
    f1est1(k)=f1b - f1b*a10(k)/2;
    f2est1(k)=f2b - f2b*a20(k)/2;
    % Load of your own transmission line model with [length,area]
    % [f1(k), f2(k)]=TLM([ones(n,1).*L./n,expif(P)'], no_vibration);
end;

figure(1)
axis([250 2650 0 1000]);  axis ij
xlabel('F2 (Hz)')
ylabel('F1 (Hz)')
hold on
plot(f2,f1, 'b'); plot(f2(indthetab),f1(indthetab), 'bO'); % C2
plot(f2est,f1est, 'c'); plot(f2est(indthetab),f1est(indthetab), 'cO'); % "est"
plot(f2est0,f1est0, 'g'); plot(f2est0(indthetab),f1est0(indthetab), 'gO'); % SE1
plot(f2est1,f1est1, 'g'); plot(f2est1(indthetab),f1est1(indthetab), 'gO'); % SE2

plot(f2b,f1b,'bO');
text(f2b-150,f1b-50,'(f1n,f2n)');
return

%%%% Transformation function used optionally 
function [om, psi1, psi2]=ijk2psi(ijk)
[a, b]=size(ijk);
co=zeros(a,b);
phase=[pi/3 pi 5*pi/3];
for k=1:a 
    om(k)=mean(ijk(k,:));
    P1=(ijk(k,1)-ijk(k,3))*sin(pi/3); P2=0.5*(ijk(k,1)+ijk(k,3))-ijk(k,2);
    psi2(k)=atan(P1/P2);
    ind=find(ijk(k,:)-om(k)); % first non zero element within i,j,k
    psi1(k)=(ijk(k,ind(1))-om(k))/cos(psi2(k)-phase(ind(1)));
end;

function vecout=expif(vecin)
vecout=(vecin < 1).*exp(vecin-1)+(vecin >=1).*vecin;

