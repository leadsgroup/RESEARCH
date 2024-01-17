naca4  = '2412';
alpha  = 2;
npanel = 20;


[cl,cd,cm,x,y,cp] = hess_smith(naca4,alpha,npanel);
cl
cd
cm