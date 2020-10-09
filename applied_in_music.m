% ver: https://en.wikipedia.org/wiki/Least_mean_squares_filter#Normalized_least_mean_squares_filter_(NLMS)
% assumirei ausência de interferencia entre a saí­da do filtro u e a resposta d
clear;
close all

% Objetivo: fazer um filtro cuja saída y para a excitação x seja igual a resposta desejada d 

% pkg load control
P=0;% número certo de zeros
Q=1;% número certo de polos

% coeficientes dos multiplicadores no circuito
% u(1:P+1): coeficientes de feed forward
% u(P+2:end): coeficientes de feedback
% u = [b0 .. bP a1 .. aQ]
% y(n) = b0x(n)+..+bPx(n-P)+a1y(n-1)..aQy(n-Q)
u=[0.5 -0.5]
%u=[0.5 sqrt(2)/2 0.5 -1 -2]

[original,fs] = audioread('Rise From The Ashes.ogg');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro
%x=ones(1,62000); if step function is applied, does not converge to u 

pkg load signal % para usar downsample()
downsample_factor = 2;
x = downsample(x,downsample_factor);
% esse algoritmo precisa que xN != 0
% para acelerar a convergência, pode-se aumentar a potência da entrada
min_x=3200;
max_x=30000;
x = x(min_x:max_x);

d=filter(u(1:P+1),[1 -u(P+2:end)],x);% d of desired response, same length as x

Pmax=0;
Qmax=1;
tol=1e-13;

%%% Adaptive Filter with constant step size %%%%%%%%

[y1,w1,filters1,err1,step1,n1] = adaptive_filter(x,d,Pmax,Qmax,tol,100);

disp('Número de iterações usadas:')
disp(n1)

N=length(w1);
hfilter1 = zeros(1,N);% handle of i-th coefficient line (plot)
for i=1:N
  % subplot para comparação entre os resultados obtidos pela variação do step
  subplot(2,4,i);
  hfilter1(i)=plot(filters1(:,i,1:n1));
  hold on
  if (i <= length(u))
    plot(u(i)*ones(1,n1))
  else
    plot(zeros(1,n1))
  end
  string = sprintf('%iº coeficiente',i);
  title(string)
  grid on
end
subplot_fig = gcf;

disp('Filtro desconhecido');
disp(u(1:P+1))
disp(u(P+2:end))
disp('Filtro calculado');
disp(w1(1:Pmax+1))
disp(w1(Pmax+2:end))

##herr1=stem(err(1:n1));% handle of error line (plot)
##title('Erro')
##grid on

% algumas vezes pode ser mais fácil observar em escala logaritmica
subplot(2,4,N+1)
herrdB1=plot(20*log10(abs(err1(1:n1))));% handle of error dB line (plot)
title('Erro em dB')
grid on

% checking these plots, we see if our filter rapidly converges to the
% unknown filter

figure
hold on
plot_step=10^(round(log10(n1))-1);
for i=1:plot_step:n1
    plot(filters1(:,:,i));
end
stem([u zeros(1,N-length(u))])
title('filtros: \mu constante')

figure
hd=plot(d,'-b');% handle of d signal line (plot)
hold on

% handle of d signal line (plot)
hy1=plot(filter(w1(1:Pmax+1),[1 -w1(Pmax+2:end)],x),'-r');

title('Respostas: \mu constante')
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
hstep1=plot(step1(1:n1),'-r');%handle to step line (plot)
title('step size \mu')
legend('entrada','step \mu');

set(groot,"currentfigure",subplot_fig);
subplot(2,4,N+2)
hstepsubplot1=copyobj(hstep1,gca);
grid on;

%%% Adaptive Filter with adapted step size %%%%%%%%

[y2,w2,filters2,err2,step2,n2] = adaptive_filter(x,d,Pmax,Qmax,tol);

disp('Número de iterações usadas:')
disp(n2)

N=length(w2);
set(groot,"currentfigure",subplot_fig);

hfilter2 = zeros(1,N);% handle of i-th coefficient line (plot)
for i=1:N
  % subplot para comparação entre os resultados obtidos pela variação do step
  subplot(2,4,i+4);
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
subplot(2,4,N+1+4)
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

set(groot,"currentfigure",subplot_fig);
subplot(2,4,N+2+4)
hstepsubplot2=copyobj(hstep2,gca);
grid on;

%%%%%% plots na mesma figura, para comparação
hcoeff=zeros(1,N);% handle da figura
for i=1:N
  % figure para comparação entre os resultados obtidos pela variação do step
  hcoeff(i)=figure;
  plot(filters1(:,i,1:n1));
  hold on
  plot(filters2(:,i,1:n2));
  if (i <= length(u))
    plot(u(i)*ones(1,max(n1,n2)))
  else
    plot(zeros(1,max(n1,n2)))
  end
  string = sprintf('%iº coeficiente',i);
  title(string)
  grid on
  legend('\mu constante','\mu variável');
end

% figure para comparação entre os resultados obtidos pela variação do step
herrdB=figure;% handle da figura
ax=axes;
herrdB1_new=copyobj(herrdB1,ax);
hold on;
herrdB2_new=copyobj(herrdB2,ax);
set(herrdB2_new,"color",[1 0 0]);
title('Erro em dB')
grid on
legend('\mu constante','\mu variável');