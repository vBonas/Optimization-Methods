clc
clear all
close all

##=====================

%Início;
tol = 0.0001;
h = 0.001;

%Passo K;
dk = 0.001;
kmax = 10-dk;
Tpontos = round(kmax/dk);
i = 0;
w0 = -2;

%Armazenamento;
W = zeros (Tpontos, 1);
J = zeros (Tpontos, 1);

for k = 0:dk:kmax;

   w = w0+k;
  w1 = w+h;

  j = (w^2)-(2*w)+2;
 j1 = (w1^2)-(2*w1)+2;

  delta = j-j1;

  i = i+1;

  if abs(delta)<tol;
    w = w1;
##    for k2 = k:Tpontos;
##     W(k2) = w;
##     J(k2) = w^2-2*w+2;
##    endfor;
    break;
  endif;

  if delta<0;
    h = -h/2;
  endif;

  W(i) = w;
  J(i) = j;

endfor;

plot(W,J);

