funcprot(0);
K = 1;
lambda = 1;
mu = 3;
rho = lambda/mu;
// rho < 1 => Traffic non bouché, rho > 1 : Saturation de la file d'attente.

//Question 1

//fonction qui simule l'évolution de la chaine de Markov
//elle part d'un état i, puis donne le temps d'attente et l'évolution suivante 
//de la chaine
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

//fonction qui simule l'évolution de la chaine de Markov avec un nombre d'états différents
//à donner
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
    //tracer de l'évolution de la chaine
    //plot2d2(cumsum(T),X);
endfunction

//fonction de simulation selon les valeurs de rho
numwin=0;
//rho<1
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov pour rho<1","temps","Etats de la chaine");

//rho=1
lambda=1;
mu=1;
numwin = numwin + 1;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov pour rho=1","temps","Etats de la chaine");

//rho>1
lambda=5;
mu=1;
numwin = numwin + 1;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov pour rho>1","temps","Etats de la chaine");


    
    //Question 2
    
    //on remet des valeurs de lambda et mu pour avoir rho<1
    lambda=1;
    mu=3;
    //fonction qui simule la chaine de Markov et garde les temps de saut et les etats
    //pour faire l'intégrale discrète comme un produit scalaire
    
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
        //on effectue l'intégrale sur les bons indices de temps et d'espace
        m=T(2:$)*X(1:$-1)';
        //on divise par la bonne durée
        m=m/sum(T(2:$));
    endfunction
    
    //variance
    function [Var]=var_ergodique(Tf)
        [X,T]=simul_markov_ergo(0,Tf);
        m=T(2:$)*X(1:$-1)';
        m=m/sum(T(2:$));
        //calcul de la variance
        Var=T(2:$)*((X(1:$-1)-m).^2)';
        Var=Var/sum(T(2:$));
        
    endfunction
    
    
    //Test des valeurs numériques
    Tf = 50000;
    E = int_ergodique(Tf);
    disp(E);
    disp(abs(E - (rho/(1-rho))));
    
    Varxt = var_ergodique(Tf);
    disp(Varxt - (rho/(1-rho).^2));
    //Evolution de la précision en fonction de rho

//on fait varier mu de 1.1 à 10 pour avoir l'évolution de la précision en fonction de rho
function [PrecisionE,PrecisionVar,rho]=tracer_precision(Nombre_rho)
    PrecisionE=zeros(1,Nombre_rho);
    PrecisionVar=zeros(1,Nombre_rho);
    vecteur_mu=ones(1,Nombre_rho);
    vecteur_mu=cumsum(vecteur_mu);
    vecteur_mu(1,1)=1.1;
    vecteur_rho=1./vecteur_mu;
    for i=1:Nombre_rho do
        lambda=1;
        mu=vecteur_mu(1,i);
        rho_test=vecteur_rho(1,i);
        //calcul par integrale de E
        E = int_ergodique(Tf);
        PrecisionE(1,i)=abs(E - (rho_test/(1-rho_test)));
        //calcul de la variance
        Varxt = var_ergodique(Tf);
        PrecisionVar(1,i)=abs(Varxt - (rho_test/(1-rho_test) .^2));
    end
    rho=vecteur_rho;
endfunction


//tracer des precision numeriques
[PrecisionE,PrecisionVar,rho]=tracer_precision(10);
numwin=numwin+1;
//Precision de E
xset("window",numwin);
clf(numwin);
plot2d2(rho,PrecisionE);
xtitle("evolution de la precision du calcul de l esperance en fonction de rho","valeurs de rho","Valeur de la precision sur E");

//Tracer de la precision de la variance
numwin=numwin+1;
//Precision de Var
xset("window",numwin);
clf(numwin);
plot2d2(rho,PrecisionVar);
xtitle("evolution de la precision du calcul de la variance en fonction de rho","valeurs de rho","Valeur de la precision sur Var");


// Question 3

//pour tracer la distribution, on utilise la fonction indicatrice de l'evenement
lambda=1;
mu=3;
rho = lambda/mu;
function [distribution]=distrib_ergodique(Tf)
    [X,T]=simul_markov_ergo(0,Tf);
    max_distrib=max(X);
    valeur_distrib=zeros(1,max_distrib+1);
    for i=0:max_distrib do
        //valeur du produit des indicatrices par les temps
        esp=T(2:$)*((X(1:$-1)==i)*i)';
        //moyenne temporelle
        esp=esp/sum(T(2:$));
        valeur_distrib(1,i+1)=esp;
    end
    disp(esp);
    distribution=valeur_distrib;
endfunction

distribution=distrib_ergodique(Tf);
//tracer de la distribution
numwin=numwin+1;
xset("window",numwin);
clf(numwin);
plot(distribution);
xtitle("distribution de X en regime stationnaire","valeurs des etats","Probabilité des états");




//Question subsidiaire 1

//on modélise la file d'attente avec des lois uniformes et K=1
//Les lois uniformes sont: 2/lambda pour l'arrivée d'un client
//2/mu pour les temps de services
//on va simuler les deux lois uniformes pour les temps de changement de la chaine
//et la variable realisant le minimum determinera si on augmente ou diminue le nombre
//de clients
lambda=1;
mu=5;
function [res, t]=evol_markov_Uniforme(i)
    res = 0;
    t1=0;
    t2=0;
    if (i==0) then
        t = grand(1,1,"unf",0,2/lambda);
    else
        t1 = grand(1,1,"unf",0,2/lambda);
        t2=grand(1,1,"unf",0,2/mu);
        t=min(t1,t2);
    end
    //si le minimum du temps d'attente est réalisé en t1, (loi d arrivée)
    //on augmente le nombre de client, sinon on le diminue
    if (t1<t2| (i==0)) then
        res = i + 1
    else
        res = i - 1
    end
endfunction



//fonction qui simule l'évolution de la chaine de Markov avec un nombre d'états différents
//à donner
function [X,T]=simul_markov_Uniforme(N, xini)
    T = [0];
    X = [xini];
    xetat = xini;
    for i=1:N do
        [res,t] = evol_markov_Uniforme(xetat);
        xetat = res;
        T = [T, t];
        X = [X, xetat];
    end
    //tracer de l'évolution de la chaine
    //plot2d2(cumsum(T),X);
endfunction


//fonction de simulation selon les valeurs de rho
numwin=numwin+1;
//rho<1
lambda=1;
mu=5;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov_Uniforme(100,0);
disp("X vaut");
disp(X);
disp("T vaut");
disp(T);
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov (Lois Uniformes temps d arrivée service) pour rho=1/5","temps","Etats de la chaine");

//rho=1
lambda=1;
mu=1;
numwin = numwin + 1;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov_Uniforme(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov (Lois Uniformes temps d arrivée service) pour rho=1","temps","Etats de la chaine");

//rho>1
lambda=5;
mu=1;
numwin = numwin + 1;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov_Uniforme(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov (Lois Uniformes temps d arrivée service) pour rho=5","temps","Etats de la chaine");





//Question subsidiaire 2 
//on suppose que tous les guichets on le même mu
//ICI rho vaut lambda/K*mu comme dans le livre de cours 
//dans le modèle à K guichets, tant que la file a 0 clients on simule l'arrivée d'un client
//avec un temps d'attente de Loi exponientielle de parametre lamdba
//et dès que l on plus de 1 clients (i clients) on a un client qui arrive avec probabilité lammbda/(lambda+(i^K)*mu)
//et un client en moins avec proba (i^K)*mu/lambda+K*mu
K=5;
function [res, t]=evol_markov_K(i)
    res = 0;
    if (i==0) then
        t = grand(1,1,'exp',1/(lambda));
    else
        t = grand(1,1,'exp',1/(lambda+min(K,i)*mu));
    end
    if (rand() <= lambda/(lambda+min(K,i)*mu) | (i==0)) then
        res = i + 1
    else
        res = i - 1
    end
endfunction



//fonction qui simule l'évolution de la chaine de Markov avec un nombre d'états différents
//à donner
function [X,T]=simul_markov_K(N, xini)
    T = [0];
    X = [xini];
    xetat = xini;
    for i=1:N do
        [res,t] = evol_markov_K(xetat);
        xetat = res;
        T = [T, t];
        X = [X, xetat];
    end
    //tracer de l'évolution de la chaine
    //plot2d2(cumsum(T),X);
endfunction


//fonction de simulation selon les valeurs de rho
numwin=numwin+1;
//rho<1
lambda=1;
mu=5;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov_K(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov pour rho=1/25, K=5","temps","Etats de la chaine");

//rho=1
lambda=5;
mu=1;
numwin = numwin + 1;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov_K(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov pour rho=1, K=5","temps","Etats de la chaine");

//rho>1
lambda=10;
mu=1;
numwin = numwin + 1;
xset("window",numwin);
clf(numwin);
[X,T]=simul_markov_K(100,0)
plot2d2(cumsum(T),X);
xtitle("evolution de la chaine de Markov pour rho=2, K=5","temps","Etats de la chaine");