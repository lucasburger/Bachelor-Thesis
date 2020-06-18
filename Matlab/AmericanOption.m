function result = AmericanPut(K,r,D0,sigma,T)
% Funktion für den Wert amerikanischer Call- und Put-Optionen auf eine
% dividendenzahlende Aktie S.


a = 5; N = 50; M = 5;       
T0 = sigma^2*T/2;           % transformierte Endzeit
h = 2*a/N;                  % Ortsschrittweite
s = T0/M;                   % Zeitschrittweite

% Abkürzungen
alpha = s/h^2; k = 2*r/sigma^2; k0 = 2*(r-D0)/sigma^2; x = [-a:h:a]; S = K*exp(x);

% Definition von f(i,j)
f = nan(M+1,N+1);
for j = 1:M+1
    f(:,j) = exp(0.5*(k0-1)*x + 0.25*(k0-1)^2*(j-1)*s*ones(1,N+1) + k*(j-1)*s).*max(0,1-exp(x));
end
u = f(:,1)';

A = diag(-0.5*alpha,alpha + 1,-0.5*alpha);
B = diag(0.5*alpha,-alpha+1,0.5*alpha);
% Berechnen der Cholesky-Zerlegung
G = chol(A,'lower');

% Schleife für die put-Option
for j = 2:M+1
    d = zeros(N+1,1);
    d(1)  = 0.5*alpha*f(1,j-1)+0.5*alpha*f(1,j);
    d(N+1) = 0.5*alpha*f(N+1,j-1)+0.5*alpha*f(N+1,j);
    b = B*u + d;
    Gprime=G';
    % Rückwärtsiteration zum Lösen von Gprime*w=b
    w = zeros(N+1,1);
    for i = n:-1:1                                 
        for k = i+1:n
            b(i) = b(i) - Gprime(i,k) * w(k) ;
        end 
        w(i) = b(i) / Gprime(i,i) ;
    end 
    
    % Vorwärtsiteration zum Lösen von G*u=w
    u_new = zeros(N+1,1);
    for i = 1:n                                   
        for k = 1:i-1                              
            w(i) = w(i) - G(i,k) * u_new(k) ;
        end 
        u_new(i) = w(i) / G(i,i) ;
    end 
    u = u_new;
end

% Rücktransformation
P = k*

result = 
end

