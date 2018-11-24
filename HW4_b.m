q = 1.6e-19 ;
Ks = 11.9 ; 
Kox = 3.95 ; 
eps0 = 8.85e-12 ;
ni = 1.45e16 ; 
W = 1.0e-6 ; 
L = 14.0e-9 ; 
tox = 1.0e-9 ; 
meff = 100.0e-4 ; 
Na = 6.0e24 ; 
sb = 00 ; 
Vfb = -1.0 ;
T = 600 ;
k = 1.38e-23 ;
delta = k*T/q ;
Vg = 0.6 ;
%Vg = 0.8 ;
Vd_init = 0;
Vd_final = 0.8;
Vd_step = 0.025;
Eg = 1.12;
Cox = eps0*Kox/tox;
gamma = sqrt(2*eps0*Ks*k*T*Na/Cox);
i = 1;

for Vd = Vd_init:Vd_step:Vd_final
     x(i) = Vd;
    lhs = Vg;
    dVds = Vd/20;
    int = 0;
    for Vd = 0:dVds:Vd
        for y = -2*Eg:0.00001:2*Eg
       
%           A = (2*k*T*Na)/eps0/Ks;
%           B = exp(-q*y/k/T)+(q*y/k/T)-1;
%           C = (ni/Na)^2;
%           D = (exp(-q*Vd/k/T))*(exp(q*y/k/T)-1);
%           F = -(q*y/k/T);
%           E = sqrt(A*(B+C*(D+F)));
            if Vg > Vfb
                rhs = Vfb+y+gamma*sqrt((exp(-q*y/k/T)+(q*y/k/T)-1)+(ni*ni/Na/Na)*(exp(q*(y-Vd)/k/T)-(q*y/k/T)-1));
                if (lhs/rhs >= 0.99 && lhs/rhs <= 1.01)
                    Fis = y;    
                end
            elseif Vg < Vfb
                rhs = Vfb+y-gamma*sqrt((exp(-q*y/k/T)+(q*y/k/T)-1)+(ni*ni/Na/Na)*(exp(q*(y-Vd)/k/T)-(q*y/k/T)-1));
                if (lhs/rhs >= 0.99 && lhs/rhs <= 1.01)
                    Fis = y;
                end
            end
        end  
        dFi = abs(Fis-delta)/20;
        inty = 0;
        for Fi = (delta/10):dFi:Fis
%       A = (2*k*T*Na)/(eps0*Ks);
%       B = exp(-q*Fi/k/T)+(q*Fi/k/T)-1;
%       C = (ni/Na)^2;
%       D = (exp(-q*Vd/k/T))*(exp(q*Fi/k/T)-1);
%       F = -(q*Fi/k/T);
%       E = sqrt(A*(B+C*(D+F)));
            numerator = (ni*ni/Na)*exp(q*(Fi-Vd)/k/T);
            denominator = sqrt((2*k*T*Na/eps0/Ks)*((exp(q*(-Fi)/k/T)+q*Fi/k/T-1)+(ni*ni/Na/Na)*(exp(q*(-Vd)/k/T)*(exp(q*Fi/k/T)-1)-q*Fi/k/T)));
            func = numerator/denominator;
            inty = inty + func*dFi;
        end
        int = int+inty*dVds;
    end
   
    Id = (q*meff*W/L)*int;
    rania(i) = Id*1.0e06;
    i = i+1; 
end
%plotting
figure('Color','White');
h1 = plot(x,rania);
%h2 = plot(x,rania);
% Line options
set(h1,'linestyle','-');
set(h1,'linewidth',2)
set(h1,'markersize',8)
set(h1,'marker','o','markeredgecolor','r','markerfacecolor','w')
% Figure options
set(gcf,'Colormap',pink);
% Axis options
axis tight;
set(gca,'color',[1 1 1]*1.0);
set(gca,'fontsize',12);
set(gca,'layer','top');
set(gca,'linewidth',1);
set(gca,'TickLength',[0.025 0.03]);
xlabel('DRAIN VOLTAGE, V_G_S [V]','Fontsize',16);
ylabel('(CURRENT), (I_D) [\muA]','Fontsize',16);
title('(Id)-Vd','Fontsize',16)
disp(h1)
grid on