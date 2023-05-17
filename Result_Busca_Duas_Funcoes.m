
##=========================

function j = J(g = zeros(6,1), x = zeros(5,1), r=zeros(6,1))

   vq0 = g(1);
   vd0 = g(2);
   vq1 = g(3);
   vd1 = g(4);
   vq2 = g(5);
   vd2 = g(6);

   deltau = [vq0; vd0; vq1; vd1; vq2; vd2];

   G = [0 , 0.0250, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0;
        0, 0.0499, 0, 0.0250, 0, 0;
        0, 0, 0, 0, 0, 0;
        0, 0.0749, 0, 0.0499, 0, 0.0250;
        0, 0, 0, 0, 0, 0];

   Phi = [0, 0.9988, 0, 1, 0;
          0, 0, 1, 0, 1;
          0, 1.9963, 0, 1, 0;
          0, 0, 2, 0, 1;
          0, 2.9925, 0, 1, 0;
          0, 0, 3, 0, 1];

   F = Phi*x;
   R = r;

   Qq0N = [200, 0, 0, 0, 0, 0;
            0, 1000, 0, 0, 0, 0;
            0, 0, 200, 0, 0, 0;
            0, 0, 0, 1000, 0, 0;
            0, 0, 0, 0, 200, 0;
            0, 0, 0, 0, 0, 1000];

   GamapiN = [0, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0;
              0, 0, 0, 1, 0, 0;
              0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 1];

   j = (G*deltau + F - R)'*Qq0N*(G*deltau + F - R) + deltau'*GamapiN*deltau;

endfunction

function busca = result_opt(x = zeros(5,1), r = zeros(6,1))

 #Passo K;
  dk = 0.001;
  kmax = 10-dk;
  Tpontos = round(kmax/dk);
  lambda0 = 3;
  hx = 10;
  tol = 0.00001;
  n = 6;

  #Início;
  vq0 = -2;
  vq1 = -2;
  vq2 = -2;
  vd0 = -2;
  vd1 = -2;
  vd2 = -2;
  w = [vq0; vd0; vq1; vd1; vq2; vd2];
  dx = [0.1; 10; 10; 10; 10; 10];
  i = 1;

  h = [1; 1; 1; 1; 1; 1];

  dw = [0; 0; 0; 0; 0; 0];

  for k = 0:dk:kmax

       j = J(w, x, r);
       i = 1;

      while (i <= n)

        dw(i) = w(i) + dx(i);

        dj = J(dw, x, r);

           if dj <= j
             h(i) = dx(i);
             else
               dx(i) = -dx(i);
               dw(i) = w(i) + dx(i);

               dj = J(dw, x, r);

             if dj <= j
               h(i) = dx(i);
             else
               dx(i) = dx(i)/2;
               h(i) = dx(i);
             endif
           endif

           if abs(dx(i)) <= tol
                h(i) = 0;
           endif
                  i = i+1;

       endwhile

               if (sqrt(h(1)^2+h(2)^2+h(3)^2+h(4)^2+h(5)^2+h(6)^2) < tol) #norma
                 break
               endif

               for kk = 0:dk:kmax

                  lambda = lambda0+hx;

                  jl0 = J(w+h*lambda0);
                  jl1 = J(w+h*lambda);

                    delta = jl0-jl1;

                    if abs(delta)<tol
                      lambda_otimo = lambda;
                      break
                    endif

                    if delta<0
                      hx = -hx/2;
                    endif

                    lambda0 = lambda;

              endfor
                w = w+lambda*h;

  endfor

  busca = w;

endfunction

##delta_otimo = result_opt
w = zeros(6,1);

x = [0.5; 0; 0; 0; 1.8667];

Phi = [0, 0.9988, 0, 1, 0;
       0, 0, 1, 0, 1;
       0, 1.9963, 0, 1, 0;
       0, 0, 2, 0, 1;
       0, 2.9925, 0, 1, 0;
       0, 0, 3, 0, 1];

Kmpc = [0 , 0.0386, 0, 0.1542, 0, 0.3469;
        2.3942 , 0, 3.3897, 0, 3.8974, 0];

R = [0; 1.2444; 0; 1.244; 0; 1.2444];

ua = Kmpc*(R-Phi*x)
ub = result_opt(w,x,R)
ub1 = ub(1)
ub2 = ub(2)
