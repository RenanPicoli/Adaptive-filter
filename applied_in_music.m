% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saída do filtro u e a resposta d
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

n=1; % index of iteration being performed/ sample being taken into account

Pmax=0;
Qmax=1;
N=Pmax+Qmax+1;%input("Entre com o número de coeficientes que deseja usar:")
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last Pmax+1 inputs AND last Qmax outputs
filter_mat=zeros(1,N,L);
alfa=zeros(Pmax+1,L);
beta=zeros(Qmax,L);
step=zeros(1,L);% this parameter is adjusted to accelerate convergence

% para convergir: erro percentual entre dois filtros consecutivos 
tol = 1e-13;% |h(n)-h(n-1)| / |h(n)|

% itera sobre as amostras
for n=1:L % cálculo de filtro em n+1 usando filtro em n
    %xN contém as últimas Pmax+1 entradas e Qmax saídas
    % atualiza com última entrada
    xN = [x(n) xN(1:Pmax) xN(Pmax+2:N)];
    step(n)=min(1/(2*xN*xN.'),10000);% this parameter is adjusted to accelerate convergence
    
    yn = filter_mat(:,:,n)*xN.';% saída atual
    % aproximação do erro: xN é aproximação do sinal completo
	  err(n) = d(n) - yn;

    if n>Qmax
      for i=1:Pmax+1 
        alfa(i,n) = xN(i) + filter(filter_mat(:,Pmax+2:end,n),[1],alfa(i,n-Qmax:n-1))(Qmax);
      end
      for j=1:Qmax
        beta(j,n) = xN(Pmax+1+j) + filter(filter_mat(:,Pmax+2:end,n),[1],beta(j,n-Qmax:n-1)(Qmax));
      end
    else
      alfa(:,n) = xN(1:Pmax+1);% 0;
      beta(:,n) = xN(Pmax+2:end);% 0;
    end
    
    delta_filter = 2*step(n)*err(n)*[alfa(:,n)' beta(:,n)'];
    filter_mat(:,:,n+1) = filter_mat(:,:,n) + delta_filter;
    %xN contém as últimas Pmax+1 entradas e Qmax saídas
    % atualiza com a saída atual (será o y(n-1) da próxima iteração)
    xN = [xN(1:Pmax+1) yn xN(Pmax+2:N-1)];
    % testa se já convergiu
    if(norm(delta_filter)/norm(filter_mat(:,:,n)) < tol)
        break;
    end
    n=n+1;
end
n=n-1;
disp('Número de iterações usadas:')
disp(n)

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
  grid on
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

% checking these plots, we see if our filter rapidly converges to the
% unknown filter

figure
hold on
plot_step=10^(round(log10(n))-1);
for i=1:plot_step:n
    plot(filter_mat(:,:,i));
end
stem([u zeros(1,N-length(u))])
title('filtros')

figure
plot(d,'-b')
hold on
plot(filter(w(1:Pmax+1),[1 -w(Pmax+2:end)],x),'-r')
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
legend('entrada','step \mu')
