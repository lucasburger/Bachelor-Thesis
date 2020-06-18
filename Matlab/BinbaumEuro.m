function result = BinbaumEuro(S0,K,r,sigma,T,N)
% Berechnet mithilfe der Binomialmethode einen approximativen
% Wert einer europaeischen Call- und Put-Option.

% Setzen der Parameter
dt = T/N ;
u = exp(sigma*sqrt(dt));
d = 1/u ;
p = (sqrt(dt)/(2*sigma))*(r-0.5*sigma^2) + 0.5;

% Bestimmen der Endwerte des Aktienkurses S(T)
S = S0*(u.^(N:-1:0)'.*d.^(0:N)');

% Bestimmen der Optionswerte im Endzeitpunkt T
Call = max(S-K,0);  Put = max(K-S,0);

% Eine present-value-Matrix. p auf der Haupt-, (1-p) auf der 1. Nebendiagonalen
pvV = exp(-r*dt)*spdiags( [p*ones(N+1,1) (1-p)*ones(N+1,1)], [0 1], N+1, N+1 );

for i = N:-1:1 % Rueckwaertsiteration zur Bestimmung von Put/Call(0,0)
    Call = pvV(1:i,1:i+1)*Call;
    Put = pvV(1:i,1:i+1)*Put;
end

result = [Call, Put];
end




