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
sb = 0.0 ; 
Vfb = -1 ;
T = 600 ;
k = 1.38e-23 ;
delta = k.*T./q ;
%Vds = 0.1 ;
Vds = 0.8 ;
Vg_init = -0.1;
Vg_final = 0.8;
Vg_step = 0.01;
i = 1;

for Vg = Vg_init:Vg_step:Vg_final
  dVds = Vds./20;
  int = 0;
  for Vd = 0:dVds:Vds
    for y = 0:0.001:1.2
      lhs = Vg; 
      A = (2.*k.*T.*Na)./(eps0.*Ks);
      B = exp(-q.*y./k./T)+(q.*y./k./T)-1;
      C = (ni./Na).^2;
      D = (exp(-q.*Vd./k./T)).*(exp(q.*y./k./T)-1);
      F = -(q.*y./k./T);
      E = sqrt(A.*(B+C.*(D+F)));
      rhs = Vfb+y+(Ks.*tox./Kox).*E;
      if (Vg == 0)
        Fis = 0;
      elseif (lhs./rhs >= 0.99 && lhs./rhs <= 1.01)
        Fis = y;
      end
    end 
    dFi = (Fis-delta)./20; 
    for Fi = delta:dFi:Fis
      A = (2.*k.*T.*Na)./(eps0.*Ks);
      B = exp(-q.*Fi./k./T)+(q.*Fi./k./T)-1;
      C = (ni./Na).^2;
      D = (exp(-q.*Vd./k./T)).*(exp(q.*Fi./k./T)-1);
      F = -(q.*Fi./k./T);
      E = sqrt(A.*(B+C.*(D+F)));
      if E == 0;
        E = 1.0e-8;
      elseif E != 0 && Fi < 0
        E = -E;  
      end  
      int = int + ((exp(q.*(Fi-Vd)./(k.*T)))./E).*dFi.*dVds;
    end  
  end      
  if(int < 0)
    int = 0;
  end
  x(i) = Vg;
  Id = (q.*meff.*W.*(ni.^2)./L./Na).*int;
  rania(i) = Id*1.0e06;
  i = i+1;
end  

%plotting
figure('Color','White');
h1 = semilogy(x,rania);
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
xlabel('GATE VOLTAGE, V_G_S [V]','Fontsize',16);
ylabel('log(CURRENT), log(I_D) [\muA]','Fontsize',16);
title('log(Id)-Vg','Fontsize',16)
disp(h1)
