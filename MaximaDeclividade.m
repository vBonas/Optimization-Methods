clc
clear all
close all

##=========================

function j = J(w)
  w0 = w(1);
  z0 = w(2);
  j = (w0^2)+(z0^2)+(2*w0)+z0+4;
end
%Início;
w0 = -2;
z0 = -2;

w = [w0,z0];

%Passo K;

dk = 1;
kmax = 100;

dk1 = 1;
kmax1 = 10;

Tpontos = round(kmax/dk);
i = 0;
v = 0;

lambda0 = 3;
hx = 0.01;
tol = 0.0001;

%Armazenamento;
## X = zeros (Tpontos, 1);
## Y = zeros (Tpontos, 1);
## Z = zeros (Tpontos, 1);

        for k1 = 0:dk1:kmax1;

            j = J(w);
           Gj = [2*w(1)+2, 2*w(2)+1];
            h = Gj;

            i = i+1;
            v = v+1;

           for kk = 0:dk:kmax;

              lambda = lambda0+hx;

               jl0 = J(w+h*lambda0);
               jl1 = J(w+h*lambda);

                delta = jl0-jl1;

                if abs(delta)<tol;
                  lambda_otimo = lambda;
                  break;
                endif;

                if delta<0;
                  hx = -hx/2;
                endif;

                lambda0 = lambda;

            endfor;
            w = w+lambda*h;
            xx = w;
            XX(v) = w(1);
            YY(v) = w(2);
            JJ(v) = J(w);


endfor;
w
J(w)

[X,Y]=meshgrid(-1500:150:1500,-1500:150:1500);
Z=(X.^2)+(Y.^2)+(2*X)+Y.+4;

figure
surf (X,Y,Z)
hold on
plot3(XX,YY,JJ,'r')

figure
contour (X,Y,Z)
hold on
plot(XX,YY)
