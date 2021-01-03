%clear;
%close all

fs=22050;
ts=1/fs;
pkg load signal % para usar butter, downsample()
z=tf('z',ts)
N=2;% ordem do filtro
% Dados do passa-banda
fc1=7000;% freq de corte 1
fc2=10000;% freqde corte 2
[b, a]=butter(N,[fc1 fc2]*2/fs)
[z,p,g]=butter(N,[fc1 fc2]*2/fs)

% Compensation for unitary gain
b=(1/g)*b

% equivalente a freq(b,a)
##figure;
##G=tf(b,a,ts)
##% bode (G)
##[mag,pha,w]=bode(G);
##subplot(2,1,1)
##plot(w/max(w),mag2db(mag));
##grid on
##subplot(2,1,2)
##plot(w/max(w),pha);
##grid on

figure;
freqz(b,a)

% applying filter on sampled signal
[original,fs] = audioread('Rise From The Ashes.ogg');

x = original(:,1);% x foi gravada com 2 canais, vamos pegar apenas o primeiro 

pkg load signal % para usar downsample()
downsample_factor = 2;
x = downsample(x,downsample_factor);
L=length(x)
fs = fs/downsample_factor;% update fs
% esse algoritmo precisa que xN != 0
% para acelerar a convergência, pode-se aumentar a potência da entrada
min_x=3200;
max_x=300_000;
x = x(min_x:max_x);

L=length(x);
X=fft(x);
f=fs*linspace(-0.5,0.5,L);
figure;
plot(f,abs(X));
grid on
title('fft(|X|): original')

y=filter(b,a,x);
Y=fft(y);
figure;
plot(f,abs(Y));
grid on
title('fft(|Y|): filtrada')

% diplay P, Q, direct form 1 coeffs
b_direct_form_1 = b
a_direct_form_1 = -a(2:end)
P=length(b_direct_form_1)-1
Q=length(a_direct_form_1)