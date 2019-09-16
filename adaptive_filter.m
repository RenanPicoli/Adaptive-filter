% uncoment if you're using octave:
% function [y,w,n] = adaptive_filter(x=[zeros(1,6) sin(50*(1:1000*6))],u=[1 2 3],N=6,step=5e-6,tol=1e-2,max_iter=1e5)
function [y,w,n] = adaptive_filter(x,u,N,step,tol,max_iter)

close all

% Objetivo: fazer um filtro cuja sa�da y para a excita��o x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)

% pkg load control
%u: input: Entre com os coeficientes do filtro desconhecido

%N: input: Entre com o n�mero de coeficientes que deseja usar

%step: input: Entre com o tamanho do passo

%tol: input: Entre com o valor da toler�ncia no erro entre y e d (norma euclidiana)

%max_iter: input: m�ximo de itera��es 

%y: output: sa�da do FA para entrada x

%w: output: filtro obtido

%n: output: n�mero de iteracoes usadas

disp(u)
disp(N)

if length(u)<N
	u=[u zeros(1,N-length(u))];
else
	N=length(u);
	warning('N foi aumentado para %d para que os vetores d e conv(w,x) tenham o mesmo tamanho!',length(u))
end

filter=zeros(1,N,max_iter);

% x=[zeros(1,N) sin(50*(1:1000*N))];
% x=[zeros(1,N) rand(1, 1000*N)];

d=conv(u,x); % d of desired response


% do{}while() equivalent for adaptive filtering
condition=true;
n=1; % index of iteration being performed

delta=step*eye(N);
err=zeros(1,N);

while condition % c�lculo de filtro em n+1 usando filtro em n
	err(n)=norm(conv(filter(:,:,n),x) - d);

	for j=1:N
		grad(j) =  (norm(conv(filter(:,:,n) + delta(j,:),x) - d)-err(n))/step;
	end

	filter(:,:,n+1) = filter(:,:,n) - step*grad;

	err(n+1)=norm(conv(filter(:,:,n+1),x) - d);
	if(abs(err(n+1))<tol)
		disp('T�rmino por atingir precis�o desejada.')
		condition=false;
	end
	if(n+1>=max_iter)
		disp('T�rmino por n�mero de itera��es.')
		condition=false;
	end

	n=n+1;	
end

figure

plot(d,'-b')
hold on
w=filter(:,:,n)
plot(conv(w,x),'-r')

title('Respostas')
legend('filtro desconhecido','filtro adaptativo')

figure
stem(err)
title('Erro (2-norma da diferen�a)')
grid on

% algumas vezes pode ser mais f�cil observar em escala logaritmica
figure
plot(20*log10(err))
title('Erro em dB')
grid on

disp('N�mero de itera��es usadas:')
disp(n)

disp('Valor final do erro entre d e sequ�ncia gerada pelo filtro w:')
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

% uncoment if you are using octave
%endfunction
end
