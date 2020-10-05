% FIR filter identification
% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saída do filtro u e a resposta d
clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control
u=[0.5 sqrt(2)/2 0.5]%input("Entre com os coeficientes do filtro desconhecido:")

% para o método do gradiente descendente
step_grad = 1;

original = load('rise_original','-ascii');
fs = 44100;% original sampling frequency

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
d=filter(u,[1],x);% d of desired response, same length as x. Could be conv(u,x)

n=1; % index of iteration being performed/ sample being taken into account

N=6;%input("Entre com o número de coeficientes que deseja usar:")
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last N inputs
filter_mat=zeros(1,N,L);

% max_iter=5000; % máximo de iterações

% para convergir: erro percentual entre dois filtros consecutivos 
tol = 1e-6;% |h(n)-h(n-1)| / |h(n)|

% itera sobre as amostras
for n=1:L % cálculo de filtro em n+1 usando filtro em n
  
    %xN contém as últimas N entradas
    xN = [x(n) xN(1:N-1)];
	  err(n)=d(n) - filter_mat(:,:,n)*xN.';% aproximação do erro: xN é aproximação do sinal completo
    
    if (xN*xN.' ~= 0)
        delta_filter = step_grad*err(n)*xN/(xN*xN.');
        filter_mat(:,:,n+1) = filter_mat(:,:,n) + delta_filter;
        % testa se já convergiu
        if(norm(delta_filter)/norm(filter_mat(:,:,n)) < tol)
            break;
        end
    end
    n=n+1;	
end

n=n-1;
w=filter_mat(:,:,n);
disp('Filtro calculado');
disp(w)

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

figure
plot(d,'-b')
hold on
plot(conv(w,x),'-r')
title('Respostas')
xt = min_x*downsample_factor/fs:max_x*downsample_factor/fs;
% altera a unidade do eixo x de amostras para segundos
set(gca,'xticklabel',xt);
xlabel('tempo(s)')
legend('filtro desconhecido','filtro adaptativo')

figure
stem(err(1:n))
title('Erro')
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
plot(20*log10(abs(err(1:n))))
title('Erro em dB')
grid on

disp('Número de iterações usadas:')
disp(n)

figure
stem(w)
hold on
stem([u zeros(1,N-length(u))])
title('filtros')
legend('filtro calculado','filtro desconhecido')

% checking these plots, we see that our filter rapidly converges to the
% unknown filter

% figure
% hold on
% for i=1:100:n
%     plot(filter_mat(:,:,i));
% end
% stem([u zeros(1,N-length(u))])
% title('filtros')
figure
norma_erro_filtro=zeros(1,n);
for i=1:n
    norma_erro_filtro(i)=20*log10(norm(filter_mat(:,:,i) - [u zeros(1,N-length(u))]));% erro em dB
end
plot(norma_erro_filtro);
title('2-norma do erro entre os filtros em dB')
