function dxdt = tfexeu(t,x)
  
    A = -diag([0.2 1]);
    
    dxdt = A*x;
