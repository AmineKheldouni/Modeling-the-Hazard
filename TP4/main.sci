funcprot(0);
K = 1;
lambda = 1;
mu = 3;
rho = lambda/mu;
// rho < 1 => Traffic non bouché, rho > 1 : Saturation de la file d'attente.

function [res, t]=evol_markov(i)
    res = 0;
    if (i==0) then
        t = grand(1,1,'exp',1/(lambda));
    else
        t = grand(1,1,'exp',1/(lambda+mu));
    end
    if (rand() <= lambda/(lambda+mu) | i==0) then
        res = i + 1
    else
        res = i - 1
    end
endfunction

function [X,T]=simul_markov(N, xini)
    T = [0];
    X = [xini];
    xetat = xini;
    for i=1:N do
        [res,t] = evol_markov(xetat);
        xetat = res;
        T = [T, t];
        X = [X, xetat];
    end
    //plot2d2(cumsum(T),X);
endfunction

function [X,T]=simul_markov_ergo(xini, Tf)
    T = [0];
    t=0;
    X = [xini];
    xetat = xini;
    while t<Tf
        [res,t_int] = evol_markov(xetat);
        xetat = res;
        T = [T, t_int];
        X = [X, xetat];
        t=t+t_int
    end
    //plot2d2(cumsum(T),X);
endfunction
//Num() = 70;
//[X,T]=simul_markov(N,0);

// Vérification de l'espérance de X_t et sa variance :
// Espérance :
function [m]=int_ergodique(Tf)
    [X,T]=simul_markov_ergo(0,Tf);
    m=T(1:$-1)*X(1:$-1)';
    m=m/Tf;
endfunction

Tf = 10000;
E = int_ergodique(Tf);
E - (rho/(1-rho))

Varxt = variance(X);
Varxt - (rho/(1-rho).^2)

// Question 3


