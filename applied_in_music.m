clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control
u=[1 2 3];%input("Entre com os coeficientes do filtro desconhecido:")

N=10;%input("Entre com o número de coeficientes que deseja usar:")

step = 1e-4%input("Entre com o tamanho do passo:") (para cálculo do grad: não influencia o nr de iterações)

tol=1e-2;%input("Entre com o valor da tolerância no erro entre y e d (norma euclidiana):")

% para o método do gradiente descendente
step_grad = 1e-3;
max_iter=700; % máximo de iterações 

filter=ones(1,N,max_iter);

original = load('rise_original','-ascii');
fs = 44100;% original sampling frequency
% original = audioread('Rise From The Ashes.mp3');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro

% pkg load signal % para usar downsample()
downsample_factor = 10;
x = downsample(x,downsample_factor);
x = x(1:30000);
%x=[zeros(1,N) sin(50*(1:1000*N))];
% x=[zeros(1,N) rand(1, 1000*N)];

filtrada = load('rise_filtrado','-ascii');
d = filtrada(:,1);
d = downsample(d,downsample_factor);
d = d(1:30000+N-1); % 30k + N-1 para ter o tamanho de conv(w,x)
%d=conv(u,x); % d of desired response

% do{}while() equivalent for adaptive filtering
condition=true;
n=1; % index of iteration being performed

delta=step*eye(N);
err=zeros(1,N);

% método do gradiente descendente
while condition % cálculo de filtro em n+1 usando filtro em n
	err(n)=norm(conv(filter(:,:,n),x) - d);

    % cálculo do gradiente
	for j=1:N
		grad(j) = (norm(conv(filter(:,:,n) + delta(j,:),x) - d)-err(n))/step;
	end

	filter(:,:,n+1) = filter(:,:,n) - step_grad*grad;

	err(n+1)=norm(conv(filter(:,:,n+1),x) - d);
	if(abs(err(n+1))<tol)
		disp('Término por atingir precisão desejada.')
		condition=false;
	end
	if(n+1>=max_iter)
		disp('Término por número de iterações.')
		condition=false;
	end

	n=n+1;	
end

figure

plot(d,'-b')
hold on
w=filter(:,:,n);
plot(conv(w,x),'-r')

title('Respostas')
legend('filtro desconhecido','filtro adaptativo')

figure
stem(err)
title('Erro (2-norma da diferença)')
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
plot(20*log10(err))
title('Erro em dB')
grid on

disp('Número de iterações usadas:')
disp(n)

disp('Valor final do erro entre d e sequência gerada pelo filtro w:')
disp(err(n))

figure
stem(w)
hold on
stem([u zeros(1,N-length(u))])
title('filtros')
legend('filtro calculado','filtro desconhecido')