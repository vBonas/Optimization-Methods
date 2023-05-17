clc
clear all
close all

##=========================

function j = J(g)
  w0 = g(1);
  z0 = g(2);
  j = (w0^2)+(z0^2)+(2*w0)+z0+4;
end

%Início;
w0 = -20;
z0 = -20;

dx = [10, 10];
i = 1;

 h(1) = 1;
 h(2) = 1;

w = [w0,z0];
dw = [0,0];

%Passo K;
dk = 0.001;
kmax = 1-dk;
Tpontos = round(kmax/dk);
v = 0;

lambda0 = 3;
hx = 0.001;
tol = 0.000001;
n = 2;

for k = 0:dk:kmax;

       j = J(w);
       i = 1;
       v = v+1;

      while (i <= n)

        dw(i) = w(i) + dx(i);

        dj = J(dw);

           if dj <= j;
             h(i) = dx(i);
             else
             dx(i) = -dx(i);
             dw(i) = w(i) + dx(i);

             dj = J(dw);
             if dj <= j;
               h(i) = dx(i);
             else
               dx(i) = dx(i)/2;
               h(i) = dx(i);
             endif
           endif

           if abs(dx(i)) <= tol;
                h(i) = 0;
           endif
                  i = i+1;

       endwhile

               if (sqrt(h(1)^2+h(2)^2) < tol) #norma
                 break;
               endif

               for kk = 0:dk:kmax

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

%Gráfico;
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
