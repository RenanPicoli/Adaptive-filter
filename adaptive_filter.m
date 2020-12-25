%Objetivo: fazer um filtro cuja sa�da y para a excita��o x seja igual a resposta desejada d (ou seja, descobrir o filtro u usado)
%x: input: input signal
%d: input: Resposta do filtro desconhecido � entrada x
%Pmax: input: n�mero m�ximo -1 de coeficientes na por��o feed forward
%Qmax: input: n�mero m�ximo de coeficientes na por��o feedback
%tol:input: erro percentual entre dois filtros consecutivos |h(n+1)-h(n)| / |h(n)|
%varargin: input: argumento opcional: valor fixo de step size
%y: output: sa�da do FA para entrada x
%w: output: filtro obtido
%filters: output: filtros intermedi�rios calculados
%err: output: erros em cada amostra
%step: output: step size em cada amostra
%n: output: n�mero de itera��es usadas

% �ltimo argumento pode ser um valor fixo de step
function [y,w,filters,err,step,n] = adaptive_filter(x,d,Pmax,Qmax,tol,varargin)
N=Pmax+Qmax+1;
L=length(x);
err=zeros(1,L);
xN=zeros(1,N);% vector with last Pmax+1 inputs AND last Qmax outputs
filter_mat=zeros(1,N,L);
step=zeros(1,L);% this parameter is adjusted to accelerate convergence
if (nargin > 5)
  string = sprintf('Using input step=%d\n',varargin{1});
  disp(string);
end
    
% itera sobre as amostras
for n=1:L % c�lculo de filtro em n+1 usando filtro em n
    %xN cont�m as �ltimas Pmax+1 entradas e Qmax sa�das
    % atualiza com �ltima entrada
    xN = [x(n) xN(1:Pmax) xN(Pmax+2:N)];
    if (nargin == 5)
      step(n)=min(1/(2*xN*xN.'),1e4);% this parameter is adjusted to accelerate convergence
    else
      step(n)=varargin{1};
    end
  
    yn = filter_mat(:,:,n)*xN.';% sa�da atual
    % aproxima��o do erro: xN � aproxima��o do sinal completo
	  err(n) = d(n) - yn;
    
    delta_filter = 2*step(n)*err(n)*xN;% [alfa beta] approx xN: Feintuch's approximation
    filter_mat(:,:,n+1) = filter_mat(:,:,n) + delta_filter;
    %xN cont�m as �ltimas Pmax+1 entradas e Qmax sa�das
    % atualiza com a sa�da atual (ser� o y(n-1) da pr�xima itera��o)
    xN = [xN(1:Pmax+1) yn xN(Pmax+2:N-1)];
    % testa se j� convergiu
    if(norm(delta_filter)/norm(filter_mat(:,:,n)) < tol)
      if (norm(filter_mat(:,:,n))==0)
        disp('Divis�o por zero na itera��o:');
        n
       end;
      break;
    end
    n=n+1;
end
n=n-1;

w=filter_mat(:,:,n);
y=filter(w(1:Pmax+1),[1 -w(Pmax+2:end)],x);
filters = filter_mat;
end
