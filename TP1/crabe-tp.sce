clear

// Loi normale
function [x]=normale(y,m,s2)
  x=%e^(-(y-m).^2/2/s2)/sqrt(2*%pi*s2)
endfunction;

// Ouvrir le fichier de données (nombre de crabes par intervalle)
x=fscanfMat('crabe.txt');
//x=x;


// intervalles
y=.580+.002+.004*[0:28];
yM=y+.002;
ym=y-.002;
Max=25;

// Dessiner la loi normale correspondante
norm_approx = normale(y, 0.647, 0.0003638);
plot2d(y, norm_approx);

// Tracer l'histogramme
tmp = 0
for i=1:29 do
    for j=1:int(x(i)) do
        tmp = tmp +1;
        data_crabe(tmp) = y(i) + 0.002;
    end
end

histplot(y, data_crabe);

// Données
pi0=[0.3;0.6;0.1];
pi=pi0;
mu=[.6; .6; .6];
s2=[1 ;1; 1]/100;

// Algorithme EM pour les crabes
//------------------------------

N=2000;
pop = 3;
R=zeros(2*pop,1);
R(:,1)=[mu(1);mu(2);mu(3);s2(1);s2(2);s2(3)];

rho = zeros(pop, 1000);
for k=1:1000 do
    for i=1:pop do
        rho(i,:) = pi(i)*normale(data_crabe', R(i), R(i+pop));
    end;
    for i=1:pop do
        for k=1:1000 do
        rho(i,k) = rho(i,k)./sum(rho(:,k));
        end;
    end;

    for i=1:pop do
        R(i) = (rho(i,:)*data_crabe)/sum(rho(i,:));
        R(pop+i) = ((((data_crabe'-R(i)).*(data_crabe'-R(i)))*rho(i,:)')) / sum(rho(i,:));
        pi(i) = mean(rho(i,:));
    end;
    pi(pop) = 1 - sum(pi(1:pop-1));
end;

disp(size(R));
disp(pi);

// Affichages
histplot(y, data_crabe);

norm_approxEM = zeros(29,1);

for j=1:pop do
    norm_approxEM(:,1) = norm_approxEM(:,1) + pi(j)*normale(, R(j), R(pop+j))';
end;

plot(y',norm_approxEM, color='red');
