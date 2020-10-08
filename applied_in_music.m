% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saí­da do filtro u e a resposta d
clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control
P=0;% número certo de zeros
Q=1;% número certo de polos

% coeficientes dos multiplicadores no circuito
% u(1:P+1): coeficientes de feed forward
% u(P+2:end): coeficientes de feedback
% u = [b0 .. bP a1 .. aQ]
% y(n) = b0x(n)+..+bPx(n-P)+a1y(n-1)..aQy(n-Q)
u=[0.5 -0.5]
%u=[0.5 sqrt(2)/2 0.5 -1 -2]

[original,fs] = audioread('Rise From The Ashes.ogg');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro
%x=ones(1,62000); if step function is applied, does not converge to u 

pkg load signal % para usar downsample()
downsample_factor = 2;
x = downsample(x,downsample_factor);
% esse algoritmo precisa que xN != 0
% para acelerar a convergência, pode-se aumentar a potência da entrada
min_x=3200;
max_x=30000;
x = x(min_x:max_x);

d=filter(u(1:P+1),[1 -u(P+2:end)],x);% d of desired response, same length as x

Pmax=0;
Qmax=1;
tol=1e-13;
[y,w,filters,err,step,n] = adaptive_filter(x,d,Pmax,Qmax,tol,3);

disp('Número de iterações usadas:')
disp(n)

N=length(w);
hfilter = zeros(1,N);% handle of i-th coefficient line (plot)
for i=1:N
  figure
  hfilter(i)=plot(filters(:,i,1:n));
  hold on
  if (i <= length(u))
    plot(u(i)*ones(1,n))
  else
    plot(zeros(1,n))
  end
  string = sprintf('%iº coeficiente',i);
  title(string)
  grid on
end

disp('Filtro desconhecido');
disp(u(1:P+1))
disp(u(P+2:end))
disp('Filtro calculado');
disp(w(1:Pmax+1))
disp(w(Pmax+2:end))

figure
herr=stem(err(1:n));% handle of error line (plot)
title('Erro')
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
herrdB=plot(20*log10(abs(err(1:n))));% handle of error dB line (plot)
title('Erro em dB')
grid on

% checking these plots, we see if our filter rapidly converges to the
% unknown filter

figure
hold on
plot_step=10^(round(log10(n))-1);
for i=1:plot_step:n
    plot(filters(:,:,i));
end
stem([u zeros(1,N-length(u))])
title('filtros')

figure
hd=plot(d,'-b');% handle of d signal line (plot)
hold on

% handle of d signal line (plot)
hy=plot(filter(w(1:Pmax+1),[1 -w(Pmax+2:end)],x),'-r');

title('Respostas')
time_step = (max_x-min_x)*downsample_factor/(100*fs);
time_step = 10^(round(log10(time_step))-1);%"arrendonda para baixo"
xt = min_x*downsample_factor/fs:time_step:max_x*downsample_factor/fs;
% altera a unidade do eixo x de amostras para segundos
set(gca,'xticklabel',xt);
xlabel('tempo(s)')
legend('filtro desconhecido','filtro adaptativo')

figure
hx=plot(x,'-b');%handle to input line (plot)
hold on
hstep=plot(step,'-r');%handle to step line (plot)
title('step size \mu')
legend('entrada','step \mu');