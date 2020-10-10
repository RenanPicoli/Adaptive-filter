%Objetivo: fazer um filtro cuja sa�da y para a excita��o x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)
%x: input: input signal
%d: input: Resposta do filtro desconhecido � entrada x
%N: input: Nn�mero m�ximo de coeficientes do filtro
%tol: input: Toler�ncia no erro entre y e d (norma euclidiana)
%y: output: sa�da do FA para entrada x
%w: output: filtro obtido
%err: output: erros em cada itera��o (y-d)
%filter_mat: output: filtros intermedi�rios obtidos
%n: output: n�mero de iteracoes usadas

function [y,w,err,filter_mat,n] = adaptive_filter(x,d,N=6,tol=1e-5)
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last N inputs
filter_mat=zeros(1,N,L);
% para o m�todo do gradiente descendente
step_grad = 1;

% itera sobre as amostras
for n=1:L % c�lculo de filtro em n+1 usando filtro em n
  
    %xN cont�m as �ltimas N entradas
    xN = [x(n) xN(1:N-1)];
	  err(n)=d(n) - filter_mat(:,:,n)*xN.';% aproxima��o do erro: xN � aproxima��o do sinal completo
    
    if (xN*xN.' ~= 0)
        delta_filter = step_grad*err(n)*xN/(xN*xN.');
        filter_mat(:,:,n+1) = filter_mat(:,:,n) + delta_filter;
        % testa se j� convergiu
        if(norm(delta_filter)/norm(filter_mat(:,:,n)) < tol)
            break;
        end
    end
    n=n+1;
end
n=n-1;
w=filter_mat(:,:,n);
y=conv(w,x);
end