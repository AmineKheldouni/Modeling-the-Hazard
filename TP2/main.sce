function res=search_pattern(M,Text)
 indice_partial_match=0;
    res=0;
    for i = 1:length(Text)
    if part(M,indice_partial_match+1:indice_partial_match+1) == part(Text,i) then
    indice_partial_match=indice_partial_match+1;
    if indice_partial_match == length(M) then
    // printf("position=%d, indice=%d\n",i,indice_partial_match);
    res=i; return;
    Pour
    end
    else
    proposition=part(M,1:indice_partial_match)+part(Text,i);
    // on calcule le nouvel indice de matching partiel
    indice_partial_match=length(overlap(M,proposition));
    end
    //printf("position=%d, indice=%d\n",i,indice_partial_match);
    end
endfunction


function res=overlap(x,y)
    // Calcule un prefixe de longueur maximale de x qui est un suffixe de y
    for i=length(x):-1:1
       //QUESTION: prefixe = le prefixe de x de longueur i
        prefixe = part(x,1:i);
        //QUESTION: suffixe = le suffixe de j de longueur i
        suffixe = part(y,max(1,length(y)-i+1):length(y));
        if (prefixe == suffixe) then
          res=prefixe;
          return;
        end
    end
    res="";
endfunction


function res=inc(n)
    // sert à décaler de 1 les indices des matrices commencant en 0
    res=n+1;
endfunction

function Matrice=markov_chain(M, E, q)
    // Calcule la matrice de transtion de la chaine de Markov
    //
    // M : la chaine cherchée,
    //
    // E: l’aphabet,
    //
    // q: la probabité de chaque lettre supposée i.i.d. selon cette loi
    card_alphabet=length(E);
    // Construction de la matrice de transition
    Matrice=zeros(length(M)+1,length(M)+1);
    for i=0:length(M)-1
        for j=1:card_alphabet
           // On rajoute la lettre E(j) à l’état courant
            proposition=part(M,1:i)+part(E,j);
            // on calcule le nouvel etat grâce à la fonction overlap
            to=overlap(M,proposition);
            indice_etat=length(to);// l’indice de l’etat, c’est sa longueur
            // rajout de la probabilité q(j) à la transition de i -> "to"
            //QUESTION: Matrice(inc(i), inc(indice_etat)) = <À COMPLÉTER>
            Matrice(inc(i), inc(indice_etat)) = Matrice(inc(i), inc(indice_etat)) + q(j);
        end
    end
    // Lorsque l’on a atteint l’etat "length(M)", c’est gagné, on s’arrête.
    Matrice(inc(length(M)),inc(length(M))) = 1;
endfunction

M = "PP"
E = "PF"
q = [1/2;1/2]
markov_chain(M,E,q);

function [moyenne,V]=moyenne_variance_tau(M,E,q)
    // Calcule la moyenne et la variance de τ
    P = markov_chain(M,E,q);
    N = max(size(P))-1;
    Q = P(1:N,1:N);
    M=eye(N,N) - Q;
    M_moins_1=M^(-1);// calcule l’inverse la matrice M
    // Calcul de la moyenne E(τ )
    moyennes = M_moins_1* ones(N,1);
    moyenne=moyennes(1,1);
    // Calcul de la variance
    VV = M_moins_1 * 2 * Q * moyennes; // Calcul de E(τ (τ − 1))
    V=VV(1,1)+moyenne-moyenne^2;// On en déduit la variance.
endfunction;

// Tirage à Pile ou Face
alphabet="PF";
proba_alphabet=[1/2,1/2];
W="P";
// loi géomètrique on doit trouver 2 et 2
[m,v]=moyenne_variance_tau(W,alphabet, proba_alphabet)

W='FP';
// on doit trouver 4 et 4
[m,v]=moyenne_variance_tau(W,alphabet, proba_alphabet)


W='PP';
// on doit trouver 6 et 22
[m,v]=moyenne_variance_tau(W,alphabet, proba_alphabet)
// Cette loi est différente de la précedente


// Le génome, 4 bases CATG
alphabet="CATG";
proba_alphabet=[1/4,1/4,1/4,1/4];
W='CCTAAGGA'
[m,v]=moyenne_variance_tau(W,alphabet, proba_alphabet)
// Alphabet de 27 lettres
alphabet="abcdefghijklmnopqrstuvwwyz"+" ";
N=length(alphabet);
proba_alphabet=ones(1,N)/N;
W='bonjour';
[m,v]=moyenne_variance_tau(W,alphabet, proba_alphabet)
