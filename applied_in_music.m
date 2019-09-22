clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control
u=[0.5 sqrt(2)/2 0.5]%input("Entre com os coeficientes do filtro desconhecido:")

N=6;%input("Entre com o número de coeficientes que deseja usar:")

step = 1e-4;%input("Entre com o tamanho do passo:") (para cálculo do grad: não influencia o nr de iterações)

tol=1e-2;%input("Entre com o valor da tolerância no erro entre y e d (norma euclidiana):")

% para o método do gradiente descendente
step_grad = 5e-6;%valor INICIAL da constante para o gradiente descendente
max_iter=5000; % máximo de iterações 

filter=ones(1,N,max_iter);
% filter(1,:,1)=[1 2 3];

original = load('rise_original','-ascii');
fs = 44100;% original sampling frequency
% original = audioread('Rise From The Ashes.mp3');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro

pkg load signal % para usar downsample()
downsample_factor = 2;
x = downsample(x,downsample_factor);
min_x=700;
max_x=360000;
x = x(min_x:max_x);
%x=[zeros(1,N) sin(50*(1:1000*N))];
% x=[zeros(1,N) rand(1, 1000*N)];

% filtrada = load('rise_filtrado','-ascii');
% d = filtrada(:,1);
% d = downsample(d,downsample_factor);
% d = d(400:3000+N-1); % 3k + N-1 para ter o tamanho de conv(w,x)
d=conv([u zeros(1,N-length(u))],x); % d of desired response

% do{}while() equivalent for adaptive filtering
condition=true;
n=1; % index of iteration being performed

delta=step*eye(N);
err=zeros(1,N);
grad = zeros(max_iter,N);

% método do gradiente descendente
while condition % cálculo de filtro em n+1 usando filtro em n
	err(n)=norm(conv(filter(:,:,n),x) - d); % sempre positivo, tem mínimo global onde se anula

    % cálculo do gradiente
	for j=1:N
		grad(n,j) = (norm(conv(filter(:,:,n) + delta(j,:),x) - d)-err(n))/step;
	end

  if (norm(grad(n,:))!=0)
    step_grad = err(n)/(norm(grad(n,:)))^2;% análogo ao método de Newton-Raphson
  end
	filter(:,:,n+1) = filter(:,:,n) - step_grad*grad(n,:);

	err(n+1)=norm(conv(filter(:,:,n+1),x) - d);
	if(abs(err(n+1))<tol)
		disp('Término por atingir precisão desejada.')
		condition=false;
    elseif(n+1>=max_iter)
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
xt = min_x*downsample_factor/fs:max_x*downsample_factor/fs;
% altera a unidade do eixo x de amostras para segundos
set(gca,'xticklabel',xt);
xlabel('tempo(s)')
legend('filtro desconhecido','filtro adaptativo')

figure
stem(err(1:n))
title('Erro (2-norma da diferença)')
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
plot(20*log10(err(1:n)))
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

figure
hold on
for i=1:n
    plot(filter(:,:,i));
end
stem([u zeros(1,N-length(u))])
title('filtros')

figure
hold on
for i=1:n-1
    plot(grad(i,:));
end
title('grad')

% fir=zeros(3,2000);
% for i=1:2000
% fir(1,i)=filter(1,1,i);
% fir(2,i)=filter(1,2,i);
% fir(3,i)=filter(1,3,i);
% end

figure
hold on
norma_grad=zeros(1,n);
for i=1:n
norma_grad(i)=norm(grad(i,:));
end
title('evolução da norma do grad')
plot(norma_grad)