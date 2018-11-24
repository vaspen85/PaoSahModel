%d
dIoff = [8.0e3 7.0e3 7.0e-2 1.0e-5 1.0e-7 7.0e-9]
T = [100 200 300 400 500 600]
dS = [ 10 3 0.07 0.09 0.1 0.12]
%a
aIoff = [1.0e-5 2.0e-5]
Vth = [0.44 0.42]
S = [0.07 0.07]
%b
Ion = [320 1200]

subplot(2,1,1)
semilogy(T,dIoff)
xlabel('T','Fontsize',16)
ylabel('log(Ioff)','Fontsize',16)
title('log(Ioff)-T','Fontsize',16)
grid on
subplot(2,1,2)
plot(T,dS)
xlabel('T','Fontsize',16)
ylabel('S','Fontsize',16)
title('S-T','Fontsize',16)