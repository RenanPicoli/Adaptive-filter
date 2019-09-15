function [y,w,n] = adaptive_filter(x=[zeros(1,6) sin(50*(1:1000*6))],u=[1 2 3],N=6,step=5e-6,tol=1e-2,max_iter=1e5)

close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)

% pkg load control
%u: input: Entre com os coeficientes do filtro desconhecido

%N: input: Entre com o número de coeficientes que deseja usar

%step: input: Entre com o tamanho do passo

%tol: input: Entre com o valor da tolerãncia no erro entre y e d (norma euclidiana)

%max_iter: input: máximo de iterações 

%y: output: saída do FA para entrada x

%w: output: filtro obtido

%n: output: número de iteracoes usadas

disp(u)
disp(N)

if length(u)<N
	u=[u zeros(1,N-length(u))];
else
	N=length(u);
	warning("N foi aumentado para %d para que os vetores d e conv(w,x) tenham o mesmo tamanho!",length(u))
endif

filter=zeros(1,N,max_iter);

% x=[zeros(1,N) sin(50*(1:1000*N))];
% x=[zeros(1,N) rand(1, 1000*N)];

d=conv(u,x); % d of desired response


% do{}while() equivalent for adaptive filtering
condition=true;
n=1; % index of iteration being performed

delta=step*eye(N);
err=zeros(1,N);

while condition % cálculo de filtro em n+1 usando filtro em n
	err(n)=norm(conv(filter(:,:,n),x) - d);

	for j=1:N
		grad(j) =  (norm(conv(filter(:,:,n) + delta(j,:),x) - d)-err(n))/step;
	end

	filter(:,:,n+1) = filter(:,:,n) - step*grad;

	err(n+1)=norm(conv(filter(:,:,n+1),x) - d);
	if(abs(err(n+1))<tol)
		disp("Término por atingir precisão desejada.")
		condition=false;
	end
	if(n+1>=max_iter)
		disp("Término por número de iterações.")
		condition=false;
	end

	n=n+1;	
end

figure

plot(d,'-b')
hold on
w=filter(:,:,n)
plot(conv(w,x),'-r')

title("Respostas")
legend("filtro desconhecido","filtro adaptativo")

figure
plot(err)
title("Erro (2-norma da diferença)")
grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
figure
plot(20*log10(err))
title("Erro em dB")
grid on

disp("Número de iterações usadas:")
disp(n)

disp("Valor final do erro entre d e sequência gerada pelo filtro w:")
disp(err(n))

endfunction
