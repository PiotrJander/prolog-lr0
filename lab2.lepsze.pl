zrobliste(0,[]):- !.
zrobliste(A,[A|L]) :- B is A-1, zrobliste(B,L).

% Zdefiniować predykaty:

% a) suma(L, S) wtw, gdy S = suma elementów listy L

suma0([], 0).
suma0([E|L], W+E) :- suma0(L, W).

suma([], 0).
suma([E|L], W2) :- suma(L, W1), W2 is W1+E. 

sumaa(L, X) :- sumaa(L, 0, X).
sumaa([], A, A).
sumaa([E|L], A, B) :- C is A+E, sumaa(L, C, B).

%zrobliste(1000000,L), sumaa(L,X).








% b) dlugosc(L, K) wtw, gdy K = liczba elementów listy L (length/2)

dlugosc([],0).
dlugosc([_|T],X) :- dlugosc(T,Y), X is Y+1.


dlugosca(L, X) :- dlugosca(L, 0, X).	
dlugosca([], A, A).			      
dlugosca([E|L], A, B) :- C is A+1, dlugosca(L, C, B).














% c) min(L, M) wtw, gdy M jest minimalnym elementem L (L = lista
%  np. liczb całkowitych) 

minw([E],E).
minw([E|L], M) :- minw(L,M), M=<E.
minw([E|L], E) :- minw(L,M), E<M.
% złożoność??!!!!!! !











max(A,B,A) :- A >= B, !.
max(A,B,B) :- A < B.


max2(A,B,C) :- A < B -> C=B ; C=A.


min([E|L], M) :- szukaj(L, E, M).
szukaj([], M, M).
szukaj([E1|L], E2 , M) :- E1>=E2, !, szukaj(L, E2, M). 
szukaj([E1|L], E2 , M) :- E1<E2, szukaj(L, E1, M).  
%szukaj([E1|L], E2 , M) :- szukaj(L, E1, M).  









% d) odwroc(L, R) wtw, gdy R jest odwróconą listą L
%  (np. odwroc([1,2,3,4], [4,3,2,1]) - sukces) 

% naive reverse (koszt czasowy i pamieciowy: n^2), definicja odwracania

nrev([], []).
nrev([E|L], R) :-  nrev(L, LR), append(LR, [E], R).

% odwroc z akumulatorem, koszt (czas, pamiec): 1n

odwroc(L, R) :-  odwroc(L, [], R).

odwroc([], R, R).
odwroc([E|L], A, R) :-  odwroc(L, [E|A], R).

% L <- R ?
odwroc2(L,R) :- rownaDlugosc(L,R), odwroc(L, [], R).

rownaDlugosc([],[]).
rownaDlugosc([_|L1], [_|L2]) :- rownaDlugosc(L1,L2).

% palindrom

palindrom(L) :-  odwroc(L, L).   % banalne

pal(L) :-  pal(L, []).

pal(L, L).                  % parzysta dlugosc
pal([_|L], L).              % nieparzysta, np. kajak, czyli [k,a,j,a,k]
pal([E|L], A) :-  pal(L, [E|A]).


% ----------------------------------------------------------------------------

% slowo(Slowo) == Slowo= a^n b^n, n>0 i/lub n >= 0

% definicja i korzystanie:
%   slowo([k,a,j,a,k])  lub ewentualnie:
%   slowo("kajak") -- lista kodow znakow (czyli [0'k,0'a,0'j,0'a,0'k])

slowo(S) :-  slowo(S, []).

slowo([a|S], A) :-  slowo(S, [b|A]).    % wstawiamy b, nie a (sic!)
slowo(S, S).                            % wersja: n >= 0
% slowo([b|S], [b|S]).                  % wersja: n >= 1


% slowoR( S, Reszta ) wtw, gdy S = a^n b^n  * Reszta, n>=1

% wersja z akumulatorem (z append) -- nienajlepsza, tylko na poczatek

slowoR(S, R) :-  slowoR(S, [], R).

slowoR([ a | S ], A, R) :-  slowoR(S, [ b | A ], R).
slowoR([ b | S ], [ b | A ], R) :-  append(A, R, S).


% wersja bez akumulatora (bez append), tj. R jest akumulatorem (sic!)

slowoS([ a | S ], R ) :-  slowoS( S, [ b | R ] ).
slowoS(S, S).                         % wersja: n >= 0

% slowoS( [ b | S ], [ b | S ] ).     % wersja: n >= 1

% Niestety zle, 
%   sukces np. dla slowoS([b,c],[b,c]),
%       czyli akceptuje puste slowa,
%       zarowno dla zapytan ustalonych, jak i zmiennych:
%           ?- slowoS([b,c],X).
%           X = [b,c]

% Poprawka, wersja dla n >= 1:

slowoS2([a | S], R) :-  slowoS2(S, [b | R]). 
slowoS2([a, b | S], S).



% ----------------------------------------------------------------------------

% Flaga Polska - rozne warianty (rozny koszt, rozne sposobu zapisu)
%    dane - lista zlozona z b oraz c

% fp1 - koszt n^2, + pamieciowy

fp1([], []).
fp1([b|L], [b|F]) :-  fp1(L, F).
fp1([c|L], F)     :-  fp1(L, FL), append(FL, [c], F).



% fp2 - dwa akumulatory (na stale b oraz na stale c),
%        czyli klasyka: koszyki Dijkstry (flaga k kolorowa)
%      
% akumulatory: parametry wEjsciowe (dane)
%       koszt liniowy (2n)

fp2(L, F) :- fp2(L, [], [], F).

fp2([], B, C, F) :- append(B, C, F).
fp2([b|L], B, C, F) :-  fp2(L, [b|B], C, F).
fp2([c|L], B, C, F) :-  fp2(L, B, [c|C], F).



% fp3 - dwa akumulatory (na stale b oraz na stale c)
%       akumulatory: parametry wYjsciowe (wyniki)
%       koszt liniowy (2n)

fp3(L, F) :- fp3(L, B, C), append(B, C, F).

fp3([], [], []).
fp3([b|L], [b|B], C) :-  fp3(L, B, C).
fp3([c|L], B, [c|C]) :-  fp3(L, B, C).



% fp4 - jeden akumulator (na stale c)     <== rozwiazanie najlepsze
%       koszt liniowy (1n)

fp4(L, F) :- fp4(L, [], F).

fp4([], F, F).
fp4([b|L], C, [b|F]) :- fp4(L, C, F).
fp4([c|L], C, F)     :- fp4(L, [c|C], F).

% w Haskellu to by było mniej więcej coś takiego:
% pf (0:xs) ak = 0 : pf xs ak
% pf (1:xs) ak = pf xs (1:ak) 
% pf [] ak = ak
% tylko, że tu idziemy i wracamy, a w Prologu tylko idziemy + unifikacja... :)


% fp5 - jeden akumulator (na stale c), wyjsciowy (dosyc sztuczne)
%       koszt liniowy (1n)

fp5(L, F) :- fp5(L, C, F, C).

fp5([], [], F, F).
fp5([b|L], C, [b|F], R) :- fp5(L, C, F, R).
fp5([c|L], [c|C], F, R) :- fp5(L, C, F, R).


% flaga holenderska - 3 listy otwarte (Zera, Jedynki, Dwójki)

%hl(L,X) :- d(L, X,ZK, ZK,JK, JK,[]).
hl(L,X) :- d(L, ZP,ZK, JP,JK, DP,DK), ZP=X, ZK=JP, JK=DP, DK=[].
d([0|L], [0|ZP],ZK, JP,JK, DP,DK) :- d(L, ZP,ZK, JP,JK, DP,DK).
d([1|L], ZP,ZK, [1|JP],JK, DP,DK) :- d(L, ZP,ZK, JP,JK, DP,DK).
d([2|L], ZP,ZK, JP,JK, [2|DP],DK) :- d(L, ZP,ZK, JP,JK, DP,DK).
d([], Z,Z, J,J, D,D).

% złączenie list odbywa się na końcu (wersja z równaniami), albo na początku


% ----------------------------------------------------------------------------

% QuickSort (3 wersje)

% (1) definicja algorytmu, czyli z (kosztownym) appendem

qsort([], []).
qsort([X | L], S) :-
      partition(L, X, M, W),
      qsort(M, SM),
      qsort(W, SW),
      append(SM, [X|SW], S).   % a co gdyby bylo: append(SM, SW, S) ?

% partition(L, E, M, W) == podzial L wzgledem E na dwie listy:
%                             M - mniejszych lub rownych E
%                             W - wiekszych od E
% uwaga: pierwszym parametrem powinna byc lista L (indeksacja)

partition([], _, [], []).
partition([X | L ], E, [X | Lm], Lw) :-  X =< E, !, partition(L, E, Lm, Lw).
partition([X | L ], E, Lm, [X | Lw]) :-  X > E, partition(L, E, Lm, Lw).



% (2) wersja z akumulatorem (czyli bez append)

qsortA(L, S) :-  qsortA(L, [], S).

% qsort(L, A, S) ==  S = sort(L) * A  (konkatenacja)

qsortA([], S, S).
qsortA([X | L], Akum, Sort) :-
      partition(L, X, L1, L2),
      qsortA(L2, Akum, S2),
      qsortA(L1, [X | S2], Sort).



% (3) ewentualnie (jak starczy czasu) wersja jeszcze ulepszona
% tj. logarytmiczna glebokosc stosu rekurencyjnych wywolan qsort

qsortLog(L, S) :- qsortLog(L, [], S).

% qsortLog(L, A, S) S = sort(L) * A  (konkatenacja)

qsortLog([], S, S).
qsortLog([H | T], A, S) :-
     partition(T, H, M, W, 0, R),
     qsortLogPom(H, A, M, W, R, S).

% predykat pomocniczy, aby nie wykonywac wielokrotnie partition

qsortLogPom(H, A, Lm, Lw, R, S) :-
    R >= 0, !, qsortLog(Lw, A, SLw), qsortLog(Lm, [H | SLw], S).
qsortLogPom(H, A, Lm, Lw, R, S) :-
    R < 0, qsortLog(Lm, [H | SLw], S), qsortLog(Lw, A, SLw).

% partition, wersja z akumulatorem, ktory jest roznica miedzy
% liczba elementow w listach wynikowych

partition([], _, [], [], S, S).
partition([A | L], X, [A | M], D, S, W) :-
     A =< X, !, NS is S+1,
     partition(L, X, M, D, NS, W).
partition([A | L], X, M, [A | D], S, W):-
     A > X, NS is S-1,
     partition(L, X, M, D, NS, W).


% rozplaszcz(L,W)

rozplaszcz(L,W) :- rozplaszcz(L,[],[],X), odwroc(X,W).
rozplaszcz([],[],A,A) :- ! .
rozplaszcz([],[E|R],A,W) :- !, rozplaszcz(E,R,A,W).
rozplaszcz([E|L],R,A,W) :- !, rozplaszcz(E,[L|R],A,W).
rozplaszcz(X,R,A,W) :- rozplaszcz([],R,[X|A],W).

% rozplaszcz(Biezaca, Reszta, Akumulator, Wynik)
% Reszta to stos zawierający te części list, których nie możemy na raz
% obrobić na pozycji Bieżącej
% wykrzykniki (czerwone!) zapobiegają punktom powrotu
% zamiast nich można by zrobić predykat nielist:

nielist([]) :- !, fail.
nielist([_|_]) :- !, fail.
nielist(_).

%albo

nielist2(X) :- functor(X,N,K), (N,K)\=([],0), (N,K)\=('.',2).

% i użyć go w ostatniej klauzuli

list([]).
list([_|_]).

nielist3(X) :- \+(list(X)).
