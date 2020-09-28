% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saída do filtro u e a resposta d
clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control
P=1;% número certo de zeros
Q=1;% número certo de polos
%u=[0.5 sqrt(2)/2 0.5 -1 -2]%input("Entre com os coeficientes do filtro desconhecido:")
u=[0.5 -0.5]%coeficientes dos multiplicadores

original = load('rise_original','-ascii');
fs = 44100;% original sampling frequency
% original = audioread('Rise From The Ashes.mp3');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro

pkg load signal % para usar downsample()
downsample_factor = 2;
x = downsample(x,downsample_factor);
min_x=3200;% esse algoritmo precisa que xN != 0
max_x=30000;
x = x(min_x:max_x);

% filtrada = load('rise_filtrado','-ascii');
% d = filtrada(:,1);
% d = downsample(d,downsample_factor);
% d = d(400:3000+N-1); % 3k + N-1 para ter o tamanho de conv(w,x)
% d=conv([u zeros(1,N-length(u))],x); % d of desired response
d=filter(u(1:P+1),[1 -u(P+2:end)],x);% d of desired response, same length as x

n=1; % index of iteration being performed/ sample being taken into account

Pmax=1;
Qmax=1;
N=Pmax+Qmax+1;%input("Entre com o número de coeficientes que deseja usar:")
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last Pmax+1 inputs AND last Qmax outputs
filter_mat=zeros(1,N,L);
alfa=zeros(Pmax+1,L);
beta=zeros(Qmax,L);

% max_iter=5000; % máximo de iterações

% para convergir: erro percentual entre dois filtros consecutivos 
tol = 1e-9;% |h(n)-h(n-1)| / |h(n)|

% itera sobre as amostras
for n=1:L % cálculo de filtro em n+1 usando filtro em n
  
    %xN contém as últimas Pmax+1 entradas e Qmax saídas
    xN = [x(n) xN(1:Pmax) d(n) xN(Pmax+2:N-1)];
	  err(n)=d(n) - filter_mat(:,:,n)*xN.';% aproximação do erro: xN é aproximação do sinal completo

    if n>1
      for i=1:Pmax+1 
        alfa (i,n) = filter([1],[1 -filter_mat(:,Pmax+2:end,n)],xN(i));
      end
      for j=1:Qmax
        beta (j,n) = filter([1],[1 -filter_mat(:,Pmax+2:end,n)],xN(Pmax+1+j));
      end
    else
      alfa(:,1) = xN(1:Pmax+1);
      beta(:,1) = xN(Pmax+2:end);
    end
    
    delta_filter = -2*err(n)*[alfa(:,n)' beta(:,n)'];
    filter_mat(:,:,n+1) = filter_mat(:,:,n) + delta_filter;
    % testa se já convergiu
    if(norm(delta_filter)/norm(filter_mat(:,:,n)) < tol)
        break;
    end
    n=n+1;
end
n=n-1;
disp('Número de iterações usadas:')
disp(n)

##disp('Resultados intermediários:')
##disp(filter_mat(:,:,1:n))
for i=1:N
  figure
  plot(filter_mat(:,i,1:n))
  hold on
  if (i <= length(u))
    plot(u(i)*ones(1,n))
  else
    plot(zeros(1,n))
  end
  string = sprintf('%iº coeficiente',i);
  title(string)
end

w=filter_mat(:,:,n);
disp('Filtro desconhecido');
disp(u(1:P+1))
disp(u(P+2:end))
disp('Filtro calculado');
disp(w(1:Pmax+1))
disp(w(Pmax+2:end))

figure
stem(err(1:n))
title('Erro')
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
plot(20*log10(abs(err(1:n))))
title('Erro em dB')
grid on

% checking these plots, we see that our filter rapidly converges to the
% unknown filter

 figure
 hold on
 for i=1:n
     plot(filter_mat(:,:,i));
 end
 stem([u zeros(1,N-length(u))])
 title('filtros')

figure
plot(d,'-b')
hold on
plot(filter(w(1:Pmax+1),[1 -w(Pmax+2:end)],x),'-r')
title('Respostas')
xt = min_x*downsample_factor/fs:max_x*downsample_factor/fs;
% altera a unidade do eixo x de amostras para segundos
set(gca,'xticklabel',xt);
xlabel('tempo(s)')
legend('filtro desconhecido','filtro adaptativo')