funcprot(0);
n=10; // Nombre de pages
// alpha in 0.8 0.9
alpha = 0.8;

function show_adj(Adj,diameters)
    [lhs,rhs]=argn(0);
    if rhs < 2 then diameters = 30*ones(1,n);end
    graph = mat_2_graph(sparse(Adj),1,'node-node');
    graph('node_x')=300*cos(2*%pi*(1:n)/(n+1));
    graph('node_y')=300*sin(2*%pi*(1:n)/(n+1));
    graph('node_name')=string([1:n]);
    graph('node_diam')= diameters;
    graph('node_color')= 1:n;
    show_graph(graph);
    rep=[1 1 1 1 2 2 2 2 2 2 2 2 2];
    plot_graph(graph,rep);
endfunction

Adj=grand(n,n,'bin',1,0.2);
//show_adj(Adj);

// Construction de la matrice de transition P
// associ´ee `a une matrice d’adjacence.
// Pss: transition d’origine,
// P: matrice de google
// z: vecteur de teleportation
// d: vecteur vaut 1 si le degré vaut zero et 0 sinon

function [P,Pss,P1,d,z,alpha]=google(Adj)
    Pss = Adj;
    alpha = 0.8;
    d = ones(n,1);
    z = ones(1,n)/n;
    e = ones(n,1);

    for i=1:n do
        if sum(Adj(i,:)) ~= 0 then
            Pss(i,:) = Adj(i,:) / sum(Adj(i,:));
            d(i,1) = 0;
        end
    end

    P1 = Pss;
    for i=1:n do
        if sum(Adj(i,:)) == 0 then
            P1(i,:) = z;
        end
    end
    disp(size(alpha*P1))
    disp(size((1-alpha)*e*z))
    P = alpha*P1 + (1-alpha) * e * z;
endfunction

[P,Pss,Pprim,d,z,alpha]=google(Adj);

// verification que P est stochastique
sum(P,'c');

e = ones(n,1);
x= rand(n,1)
y1= P'*x;
y2= alpha*Pss'*x + (alpha*d*z)'*x+ ((1-alpha)*e*z)'*x;
disp(y1)
disp(y2)
disp(y1 - y2)

[evals,X] =spec(P');
disp(evals)
disp(X)

pi = abs(evals(:,1)/sum(evals(:,1)));

disp(pi)
disp(sum(pi))

clf();
//show_adj(Adj,int(300*pi'));

function [pi]=pi_iterative()
    p=ones(n,1);
    k = 1;
    while k < 100000
        pn = P'*p;
        k = k + 1;
        if norm(pn-p,%inf) < 10*%eps then
             break;
        end
        p = pn;
    end
    pi= pn/sum(pn);
endfunction

pi = pi_iterative();
clean(P'*pi - pi);
disp(pi)
disp(sum(pi))
disp(P'*pi - pi)


function [pi]=pi_iterative_sparse()
    p=ones(n,1);
    k = 1;
    while k < 100000
        pn = alpha*Pss'*p + (alpha*d*z)'*p+ ((1-alpha)*e*z)'*p;
        k = k + 1;
        if norm(pn-p,%inf) < 10*%eps then
             break;
        end
        p = pn;
    end
    pi= abs(p/sum(p));
endfunction

pi=pi_iterative_sparse();
clean(P'*pi - pi);
disp(pi)
disp(P'*pi- pi)

//Question 7
function []=maximizePageRank(p,m, Adj)
    Adj_copy = Adj;
    k = 1;
    PR = pi_iterative_sparse();
    PR = sum(PR(1,m));
    while k < 100000
        for i=m+1:n do
    
        end
    end
endfunction


//Question 8
function y=r(x)
    y=x.^2
endfunction

n=4;
P=rand(n,n)
pr=sum(P,'c');
P = P ./ (pr*ones(1,n));

function [cerg]=ergodique_markov_T(T,P)
    //on prend la loi initiale u uniforme de X0
    //on rappelle que la loi de Xt est (P^t)'u
    vecteur=[1:n];
    vecteur=vecteur';
    loiInit=ones(n,1)/n;
    Matrice=eye(n,n);
    Esperance=0;
    for i=0:(T-1) do
        Esperance=Esperance+((Matrice'*loiInit)')*(r(vecteur));
        Matrice=Matrice*P;
    end
    cerg=Esperance/T;
endfunction


function [cerg,pi]=ergodique_markov(P)
    p=ones(n,1);
    vecteur=[1:n];
    vecteur=vecteur';
    k = 1;
    while k < 100000
        pn = P'*p;
        k = k + 1;
        if norm(pn-p,%inf) < 10*%eps then
             break;
        end
        p = pn;
    end
    pi= pn/sum(pn);
    cerg=(pi')*r(vecteur);
endfunction

disp(ergodique_markov_T(10,P));

// test
T=100000;
CT=ergodique_markov_T(T,P);
[c,pi]=ergodique_markov(P);
disp("Test");
disp(c-CT);

// Le noyau de P-I est engendr´e par ones(n,1)
[x0,K]=linsolve(P- eye(n,n),zeros(n,1));
disp("x0");
disp(x0);
disp("K");
disp(K);


//Question 9
// le projecteur spectral sur Espace propre associ´e a 1
pi=pi';
Pr = ones(n,1)*pi; // [pi;pi;pi;....]
A = P-eye(n,n); // A -Id
S = Pr - inv(Pr-A) // Pr-A est inversible
// v´erifier que S*Pr et Pr*S sont nuls
disp("s*Pr");
clean(S*Pr);
disp(S*Pr);
disp("Pr*S");
clean(Pr*S);
disp(Pr*S);
// A*w + R - c= 0
// A*c = 0
R = r([1:n]');
// v´erifions que w=-S*R et c=Pr*R sont solution du systeme linaire
w= -S*R;
c= Pr*R;
disp("A*w+R-c");
disp(A*w + R -c);
disp("A*c");
disp(A*c);
// Noter que w n’est pas unique, on peut rajouter `a w les elts du noyau de A

// Montrons inversement que c doit ^etre egal `a Pr*R
// Pr*A est nul
disp("Pr*A");
disp(Pr*A);
// on doit donc avoir
// Pr*R - Pr*c = 0 et A*c =0
// en sommant
// Pr*R = (Pr-A)*c
// c = (Pr-A)^-1 *Pr*R
// c = (Pr-S)*Pr*R = Pr*Pr*R -S*Pr*R = Pr*R
// car Pr est un projecteur Pr^2 = Pr et S*Pr = 0
disp("Pr.^2-Pr");
clean(Pr.^2-Pr);
disp(Pr.^2-Pr);
disp("S*Pr");
clean(S*Pr);
disp(S*Pr);
// conclusion c doit valoir Pr*R
// on le v´erifie avec linsolve
[x0,K]=linsolve([A,-eye(n,n);zeros(n,n),A],[R;zeros(n,1)]);
disp("x0 et K");
disp(x0);
// on v´erifie bien que e = Pr*RK);
disp("Pr*r");
disp(Pr*R);



P1=rand(n,n);
pr=sum(P1,'c');
P1 = P1 ./ (pr*ones(1,n));
z=grand(1,n,'unf',0,1);
z=z/sum(z);
alpha = 0.8;
P = alpha*P1 + (1-alpha)*ones(n,1)*z;

// les couts Rm(i,j)
Rm = grand(n,n,'unf',0,1);

//Question 10

// On le v´erifie numeriquement
// trouver la solution de
// w = alpha*P1*w + sum(P.*Rm,’c’)
[x0,K]=linsolve(alpha*P1- eye(n,n),sum(P.*Rm,'c'));

w = x0;
disp("w");
disp(w-alpha*P1*w-sum(P.*Rm,'c'));
// calcul de c
c = (1-alpha)*z*w
// (w,c) solution du pb ergodique ?
disp("verification de la solution (w,c) trouvée ");
disp(size(P));
disp(size(R));
disp("test de w");
disp(w + c - (P*w + sum(P.*Rm,'c')));
// Maintenant on peut utiliser une m´ethode iterative


//Question 11

function [w]=iterative_c(tol)
        res1=ones(n,1);
        res2=alpha*P1*res1+sum(P.*Rm,'c');
        while(((res2-res1).^(2)/n)>tol)
            res1=res2;
            res2=alpha*P1*res2+sum(P.*Rm,'c');
        end
        w=res2;
endfunction

w=iterative_c(10*%eps);
disp("w valeur");
disp(w);
disp("test w");
disp(alpha*P1*w+sum(P.*Rm,'c')-w);
// calcul de c
c = (1-alpha)*z*w
// (w,c) solution du pb ergodique ?
disp(w + c - (P*w + sum(P.*Rm,'c')));


//Question 12
function [w]=algo_iter(tol, max_iter)
    w0 = zeros()
    w1 = ones()
    k = 0
    while (abs(w1-w0) > tol & k<max_iter)
        w0 = w1;
        k = k+1;
        // expr = expression de droite
        // Calcul de nu_k optimal
        expr = -%inf
        nu_k = 0
        for x=1:n do
            expr_tmp = ...
            if (expr_tmp > expr)
                expr = expr_tmp;
                nu_k = ...
            end
        end
        // On a nu_k optimal 
        // Résolution du système pour trouver w_{k+1} :
        w1 = ... ;        
    end
endfunction
