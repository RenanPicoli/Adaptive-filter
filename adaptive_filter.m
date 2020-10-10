%Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)
%x: input: input signal
%d: input: Resposta do filtro desconhecido à entrada x
%N: input: Nnúmero máximo de coeficientes do filtro
%tol: input: Tolerância no erro entre y e d (norma euclidiana)
%y: output: saída do FA para entrada x
%w: output: filtro obtido
%err: output: erros em cada iteração (y-d)
%filter_mat: output: filtros intermediários obtidos
%n: output: número de iteracoes usadas

function [y,w,err,filter_mat,n] = adaptive_filter(x,d,N=6,tol=1e-5)
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last N inputs
filter_mat=zeros(1,N,L);
% para o método do gradiente descendente
step_grad = 1;

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
y=conv(w,x);
end