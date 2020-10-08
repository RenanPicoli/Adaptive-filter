%Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)
%x: input: input signal
%u: input: Entre com os coeficientes do filtro desconhecido
%N: input: Entre com o número de coeficientes que deseja usar
%tol: input: Entre com o valor da tolerância no erro entre y e d (norma euclidiana)
%y: output: saída do FA para entrada x
%w: output: filtro obtido
%n: output: número de iteracoes usadas

function [y,w,n] = adaptive_filter(x,d,Pmax,Qmax,tol=1e-5)
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

w=filter(:,:,n);
y=filter(w(1:Pmax+1),[1 -w(Pmax+2:end)],x);

end
