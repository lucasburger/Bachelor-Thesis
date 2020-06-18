function result = AmericanCall(K,r,D0,sigma,T)
% Berechnung des Wertes einer Amerikanischen Put-Option. result = [S,P],
% wobei S der Aktienkurs und P der entsprechende Optionswert ist.

a = 5; N = 10000; M = 5;       
T0 = sigma^2*T/2;           % transformierte Endzeit
h = 1.2*a/N; s = T0/M;        % Schrittweiten
alpha = s/h^2; k = 2*r/sigma^2; k0 = 2*(r-D0)/sigma^2; % Abkuerzungen
x = [-a:h:0.2*a]; S = K*exp(x);

% Definition von f(i,j)
f = nan(N+1,M+1);
for j = 1:M+1
    f(:,j) = exp(0.5*(k0-1)*x + 0.25*(k0-1)^2*(j-1)*s*ones(1,N+1)... 
        + k*(j-1)*s*ones(1,N+1))'.*max(0,exp(x)-1)';
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
    % Vorwaertsiteration zum Loesen von G*w=b
    w = zeros(N+1,1);
    for i = 1:N+1                                 
        for l = 1:i-1
            b(i) = b(i) - G(i,l) * w(l) ;
        end 
        w(i) = b(i) / G(i,i) ;
    end
    
    % Rueckwaertsiteration zum Loesen von Gprime*u=w
    u_new = zeros(N+1,1);
    for i = N+1:-1:1                                   
        for l = i+1:N+1                              
            w(i) = w(i) - Gprime(i,l) * u_new(l) ;
        end
        % Kontrolle der Nebenbedingung
        u_new(i) = max(w(i) / Gprime(i,i),f(i,j));
    end
    u = u_new;
end

for i = 1:N-1
    if u(i) == f(i,M+1)
        disp(i)
        disp(K*exp(-a+i*h))
        break
    end
end


% Ruecktransformation
C = K*exp(-0.5*(k0-1)*x - 0.25*(k0-1)^2*T0 - k*T0).*u';
result = [S;C];

end