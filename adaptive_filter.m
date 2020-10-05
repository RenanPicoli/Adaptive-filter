%Objetivo: fazer um filtro cuja sa�da y para a excita��o x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)
%x: input: input signal
%u: input: Entre com os coeficientes do filtro desconhecido
%N: input: Entre com o n�mero de coeficientes que deseja usar
%tol: input: Entre com o valor da toler�ncia no erro entre y e d (norma euclidiana)
%y: output: sa�da do FA para entrada x
%w: output: filtro obtido
%n: output: n�mero de iteracoes usadas

function [y,w,n] = adaptive_filter(x,d,N=6,tol=1e-5)
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last N samples
filter=zeros(1,N,L);
step_grad=1;

% itera sobre as amostras
for n=1:L % c�lculo de filtro em n+1 usando filtro em n
    xN = [x(n) xN(1:N-1)];
	err(n)=d(n) - filter(:,:,n)*xN.';% aproxima��o do erro: xN � aproxima��o do sinal completo
    
    if (xN*xN.' ~= 0)
        delta_filter = step_grad*err(n)*xN/(xN*xN.');
        filter(:,:,n+1) = filter(:,:,n) + delta_filter;
        % testa se j� convergiu
        if(norm(delta_filter)/norm(filter(:,:,n)) < tol)
            break;
        end
    end
    n=n+1;	
end

n=n-1;
w=filter(:,:,n);
y=conv(w,x);
disp('Filtro calculado');
disp(w)
disp('N�mero de itera��es usadas:')
disp(n)
end