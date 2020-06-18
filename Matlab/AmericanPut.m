function result = AmericanPut(K,r,D0,sigma,T)
% Berechnung des Wertes einer Amerikanischen Put-Option. result = [S,P],
% wobei S der Aktienkurs und P der entsprechende Optionswert ist.

a = 5; N = 100; M = 5;       
T0 = sigma^2*T/2;           % transformierte Endzeit
h = 2*a/N; s = T0/M;        % Schrittweiten
alpha = s/h^2; k = 2*r/sigma^2; k0 = 2*(r-D0)/sigma^2; % Abkuerzungen
x = [-a:h:a]; S = K*exp(x);

% Definition von f(i,j)
f = nan(N+1,M+1);
for j = 1:M+1
    f(:,j) = exp(0.5*(k0-1)*x + 0.25*(k0-1)^2*(j-1)*s*ones(1,N+1)... 
        + k*(j-1)*s*ones(1,N+1))'.*max(0,1-exp(x))';
end
u = f(:,1);

A = (alpha + 1)*diag(ones(N+1,1)) + -0.5*alpha*diag(ones(N,1),1)...
    + -0.5*alpha*diag(ones(N,1),-1);
B = (-alpha + 1)*diag(ones(N+1,1)) + 0.5*alpha*diag(ones(N,1),1)...
    + 0.5*alpha*diag(ones(N,1),-1);

% Berechnen der Cholesky-Zerlegung
G = chol(A,'lower');

% Schleife zur Iteration bis tau = sigma^2*T/2
for j = 2:M+1
    d = zeros(N+1,1);
    d(1)  = 0.5*alpha*f(1,j-1)+0.5*alpha*f(1,j);
    d(N+1) = 0.5*alpha*f(N+1,j-1)+0.5*alpha*f(N+1,j);
    b = B*u + d; Gprime=G';
    % Rueckwaertsiteration zum Loesen von Gprime*v=b
    v = zeros(N+1,1);
    for i = N+1:-1:1                                 
        for l = i+1:N+1
            b(i) = b(i) - Gprime(i,l) * v(l) ;
        end 
        v(i) = b(i) / Gprime(i,i) ;
    end
    
    % Vorwaertsiteration zum Loesen von G*u=w
    u_new = zeros(N+1,1);
    for i = 1:N+1                                   
        for l = 1:i-1                              
            v(i) = v(i) - G(i,l) * u_new(l) ;
        end
        % Kontrolle der Nebenbedingung
        u_new(i) = max(v(i) / G(i,i),f(i,j));
    end
    u = u_new;
end

% Ruecktransformation
P = K*exp(-0.5*(k0-1)*x - 0.25*(k0-1)^2*T0 - k*T0).*u';
result = [S;P];
end