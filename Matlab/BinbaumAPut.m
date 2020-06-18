function result = BinbaumAPut(S0,K,r,sigma,T,N)
%Funktion zur Bewertung einer amerikanischen Put-Option zu gegebenen
%Parametern. result = BinbaumAPut(S0,K,r,sigma,T,N)

% Setzen der Parameter
dt = T/N;
u = exp(sigma*sqrt(dt));
d = 1/u;
p = (sqrt(dt)/(2*sigma))*(r-0.5*sigma^2) + 0.5;

S = nan(N+1,N+1); S(1,1) = S0;
% Initialisieren des Binomialbaums
for i = 2:N+1
    S(:,i) = S(:,i-1)*u;
    S(i,i) = S(i-1,i-1)*d;
end

V = max(K-S,0);

% "present value Matrix"
pvV = exp(-r*dt)*spdiags([p*ones(N+1,1) (1-p)*ones(N+1,1)], [1 0], N+1, N+1);
% Rueckwaertsiteration durch den Baum
for i = N:-1:1
    S(:,i) = max( pvV(1:i,1:i+1)*S(:,i+1), V(:,i) );
end

result = V(1,1);

end

