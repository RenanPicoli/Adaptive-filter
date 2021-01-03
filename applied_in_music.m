% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saí­da do filtro u e a resposta d
clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

%P: número certo de zeros
%Q: número certo de polos
fs=22050;
ts=1/fs;

pkg load control% para usar tf, bode, etc
##z=tf('z',ts)
##N=1;% ordem do filtro
##% Dados do passa-banda
##fc1=5000;% freq de corte 1
##fc2=10000;% freq de corte 2
##pkg load signal % para usar butter, downsample()
##[b, a]=ellip(N,5,40,[fc1 fc2]*2/fs)
##[z,p,g]=ellip(N,5,40,[fc1 fc2]*2/fs)
##
##% Compensation for unitary gain
##b=(1/g)*b

b=[1 0 -2 1]
a=[1.000000   0.590110   0.582896   0.302579   0.076053]

% diplay P, Q, direct form 1 coeffs
b_direct_form_1 = b
a_direct_form_1 = -a(2:end)
P=length(b_direct_form_1)-1
Q=length(a_direct_form_1)

% coeficientes dos multiplicadores no circuito
% u(1:P+1): coeficientes de feed forward
% u(P+2:end): coeficientes de feedback
% u = [b0 .. bP a1 .. aQ]
% y(n) = b0x(n)+..+bPx(n-P)+a1y(n-1)..aQy(n-Q)
u=[b_direct_form_1 a_direct_form_1]

[original,fs2] = audioread('Rise From The Ashes.ogg');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro
%x=ones(1,62000); if step function is applied, does not converge to u 

downsample_factor = 2;
x = downsample(x,downsample_factor);
% esse algoritmo precisa que xN != 0
% para acelerar a convergência, pode-se aumentar a potência da entrada
min_x=3200;
max_x=300_000;
x = x(min_x:max_x);

d=filter(u(1:P+1),[1 -u(P+2:end)],x);% d of desired response, same length as x

Pmax=3;
Qmax=4;
tol=1e-13;

%%% Adaptive Filter with adapted step size %%%%%%%%

[y2,w2,filters2,err2,step2,n2] = adaptive_filter(x,d,Pmax,Qmax,tol);

disp('Número de iterações usadas:')
disp(n2)

N=length(w2);

hfilter2 = zeros(1,N);% handle of i-th coefficient line (plot)
for i=1:N
  % subplot para comparação entre os resultados obtidos pela variação do step
  subplot(2,6,i);
  hfilter2(i)=plot(filters2(:,i,1:n2));
  hold on
  if (i <= length(u))
    plot(u(i)*ones(1,n2))
  else
    plot(zeros(1,n2))
  end
  string = sprintf('%iº coeficiente',i);
  title(string)
  grid on
end

disp('Filtro desconhecido');
disp(u(1:P+1))
disp(u(P+2:end))
disp('Filtro calculado');
disp(w2(1:Pmax+1))
disp(w2(Pmax+2:end))

##herr2=stem(err(1:n2));% handle of error line (plot)
##title('Erro')
##grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
subplot(2,6,[9 12])
herrdB2=plot(20*log10(abs(err2(1:n2))));% handle of error dB line (plot)
title('Erro em dB')
grid on

% checking these plots, we see if our filter rapidly converges to the
% unknown filter

figure
hold on
plot_step=10^(round(log10(n2))-1);
for i=1:plot_step:n2
    plot(filters2(:,:,i));
end
stem([u zeros(1,N-length(u))])
title('filtros: \mu variável')

figure
hd=plot(d,'-b');% handle of d signal line (plot)
hold on

% handle of d signal line (plot)
hy2=plot(filter(w2(1:Pmax+1),[1 -w2(Pmax+2:end)],x),'-r');

title('Respostas: \mu variável')
time_step = (max_x-min_x)*downsample_factor/(100*fs);
time_step = 10^(round(log10(time_step))-1);%"arrendonda para baixo"
xt = min_x*downsample_factor/fs:time_step:max_x*downsample_factor/fs;
% altera a unidade do eixo x de amostras para segundos
set(gca,'xticklabel',xt);
xlabel('tempo(s)')
legend('filtro desconhecido','filtro adaptativo')

figure
hx=plot(x,'-b');%handle to input line (plot)
hold on
hstep2=plot(step2(1:n2),'-r');%handle to step line (plot)
title('step size \mu')
legend('entrada','step \mu');