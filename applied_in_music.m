% FIR filter identification
% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saída do filtro u e a resposta d
clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control

% coeficientes dos multiplicadores no circuito
% u = [b0 .. bP]
% y(n) = b0x(n)+..+bPx(n-P)
u=[0.5 sqrt(2)/2 0.5]%input("Entre com os coeficientes do filtro desconhecido:")

[original,fs] = audioread('Rise From The Ashes.ogg');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro
%x=ones(1,62000); if step function is applied, does not converge! 

pkg load signal % para usar downsample()
downsample_factor = 2;
x = downsample(x,downsample_factor);
min_x=3200;% esse algoritmo precisa que xN != 0
max_x=30000;
x = x(min_x:max_x);

d=filter(u,[1],x);% d of desired response, same length as x. Could be conv(u,x)

% para convergir: erro percentual entre dois filtros consecutivos 
tol = 1e-6;% |h(n)-h(n-1)| / |h(n)|
N=6;%input("Entre com o número de coeficientes que deseja usar:")
L=length(x);

[y,w,err,filter_mat,n] = adaptive_filter(x,d,N=6,tol=1e-5);

disp('Número de iterações usadas:')
disp(n)

w=filter_mat(:,:,n);
disp('Filtro calculado');
disp(w)

hcoeff=zeros(1,N);
for i=1:N
  figure
  hcoeff(i)=plot(filter_mat(:,i,1:n));
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

figure
hd=plot(d,'-b');
hold on
hy=plot(y,'-r');
title('Respostas')
##xt = min_x*downsample_factor/fs:max_x*downsample_factor/fs;
##% altera a unidade do eixo x de amostras para segundos
##set(gca,'xticklabel',xt);
##xlabel('tempo(s)')
grid on
legend('filtro desconhecido','filtro adaptativo')

figure
herr=stem(err(1:n));
title('Erro')
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
herrdB=plot(20*log10(abs(err(1:n))));
title('Erro em dB')
grid on

figure
stem(w)
hold on
stem([u zeros(1,N-length(u))])
title('filtros')
legend('filtro calculado','filtro desconhecido')

figure
norma_erro_filtro=zeros(1,n);
for i=1:n
    norma_erro_filtro(i)=20*log10(norm(filter_mat(:,:,i) - [u zeros(1,N-length(u))]));% erro em dB
end
hnormaerr=plot(norma_erro_filtro);
grid on
title('2-norma do erro entre os filtros em dB')

% checking these plots, we see if our filter rapidly converges to the
% unknown filter

hfiltros=figure;
hold on
for i=1:10:n
    plot(filter_mat(:,:,i));
end
stem([u zeros(1,N-length(u))])
title('filtros')
