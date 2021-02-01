% TS217- TP Egalisation
% Pascal Vallet (IPB)
% 2017

clear all;
close all;
clc;

%% Paramètres

% Longueur de la séquence binaire transmise
N=2000; % à fixer en fonction de l'ordre de grandeur de la probabilité d'erreur à estimer
f=-1/2:1/N:1/2-1/N;

% SNR en dB
SNR=[0:20];
sigma=10.^(-SNR/20); % ecart-type du bruit correspondant

% Quelques canaux à tester ...

h=[1;0.5]; % CANAL TEST 1
%h=[1;0.1;0.9]; % CANAL TEST 2
%h=[0.5;0.8;0.5]/0.5; % CANAL TEST 3



K=length(h); % longueur du canal
P=K+2 ; % ordre des filtres ZF/MMSE/DFE (filtre direct), A COMPLETER
%Q= ; % ordre du filtre de retour DFE, A COMPLETER

%% Simulation des signaux

% Bits
bits=rand(N,1) > 0.5;
s=2*(bits-0.5); % réalisation d'une séquence de N symboles i.i.d. BPSK

% Quelques variables ...
zZF=zeros(N,length(SNR)); % sortie de l'égaliseur ZF
bitsZF=zeros(N,length(SNR)); % bits estimés après ZF
berZF=zeros(1,length(SNR)); % BER après ZF
zMMSE=zeros(N,length(SNR)); % sortie de l'égaliseur MMSE
bitsMMSE=zeros(N,length(SNR)); % bits estimés après MMSE
berMMSE=zeros(1,length(SNR)); % BER après MMSE
zDFE=zeros(N,length(SNR)); % sortie de l'égaliseur DFE
bitsDFE=zeros(N,length(SNR)); % bits estimés après DFE
berDFE=zeros(1,length(SNR)); % BER après DFE

for i=1:length(SNR) 
    
    %w=sigma(i)*(randn(N,1)+1j*randn(N,1));
    y= conv(s,h,'same'); % Observations en sortie du canal, A COMPLETER
    
    % Filtre ZF
    dZF= ; % retard optimal du filtre ZF, A COMPLETER
    fZF= ; % vecteur colonne des coefficients du filtre ZF, A COMPLETER
    
    % Filtre MMSE
    dMMSE= ; % retard optimal du filtre MMSE, A COMPLETER
    fMMSE= ;  % vecteur colonne des coefficients du filtre MMSE, A COMPLETER
    
    % Filtre DFE
    dDFE= ; % retard optimal du filtre DFE, A COMPLETER
    fDFE= ; % coefficients du filtre DFE, A COMPLETER
    
    % Egalisation ZF/MMSE
    for n=P+K-1:N 
        zZF(n,i)=real(fZF'*y(n:-1:n-P+1));
        zMMSE(n,i)=real(fMMSE'*y(n:-1:n-P+1));  
    end
    bitsZF(:,i)=real(zZF(:,i)) > 0; 
    berZF(:,i)=sum(abs(bits(P+K-1-(dZF-1):N-(dZF-1))-bitsZF(P+K-1:N,i)),1)/(N-(P+K-2)); 
    bitsMMSE(:,i)=real(zMMSE(:,i)) > 0; 
    berMMSE(i)=sum(abs(bits(P+K-1-(dMMSE-1):N-(dMMSE-1))-bitsMMSE(P+K-1:N,i)),1)/(N-(P+K-2)); 
    
    % Egalisation DFE
    M=max(P+K-1,Q+dDFE+1);
    sh=zeros(N,1);
    sh(1:M)=s(1:M).'; % initialisation des symboles estimés précédents
    for n=M:N
        zDFE(n,i)=real(fDFE'*[y(n:-1:n-P+1);sh(n-1-(dDFE-1):-1:n-(dDFE-1)-Q)]);
        sh(n-(dDFE-1))=2*(real(zDFE(n,i)) > 0)-1; % Estimation du symbole BPSK émis  
        %sh(n-(dDFE(i)-1))=s(n-(dDFE(i)-1)); % si estimation parfaite
    end 
    bitsDFE(:,i)=real(zDFE(:,i)) > 0; 
    berDFE(i)= sum(abs(bits(M-(dDFE-1):N-(dDFE-1))-bitsDFE(M:N,i)),1)/(N-M); 
end

%% Graphes

% Fonction de transfert canal/ZF

figure;
freqz(h,1,pi*f);
% Constellation des échantillons reçus y(n)
figure;
plot(real(y),imag(y),'o');
hold on;
plot([1,-1],[0,0],'x');
xlabel('partie réelle');
ylabel('partie imaginaire');
xlim([-3,3]);
legend('Symboles reçus','Symboles Bpsk');
title('Constellation des symboles reçus sans bruit additif : Canal 2');

% Courbes de BER
% Probabilité d'erreur du canal AWGN = Q(sqrt(2 Eb/N0)), où Eb=1 et N0 = sigma^2
figure;
%semilogy(SNR,berZF,'-b^',SNR,berMMSE,'-rd',SNR,berDFE,'-gv',SNR,1-normcdf(sqrt(2./sigma.^2),0,1),'k-','LineWidth',3);
semilogy(SNR,berZF,'-b^',SNR,berMMSE,'-rd',SNR,berDFE,'-gv',SNR,0.5*(1-erf(sqrt(1./sigma.^2))),'k-','LineWidth',3);
grid on;
xlabel('SNR');
ylabel('BER')
legend('ZF','MMSE','DFE','AWGN','Location','SouthWest');




