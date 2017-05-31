% Zadania na trzeci tydzien zajec z Prologu.

% Sposoby reprezentacji struktur danych
% -------------------------------------
%   (a) termowa (stale, funkcje),
%   (b) klauzulowa (klauzule unarne, i nie tylko)


% drzewa binarne, drzewa BST


% drzewo(D) == D jest drzewem binarnym

drzewo(nil).
drzewo(tree(L, _, P)) :-  drzewo(L), drzewo(P).


% insertBST( DrzewoBST, Elem, NoweDrzewoBST )

insert(nil, Elem, tree(nil, Elem, nil)).
% insert(tree(L, E, P), E, tree(L, E, P)).   % nie wstawia wielokrotnie
insert(tree(L, W, P), E, tree(NL, W, P)) :-
       E =< W, !,
       insert(L, E, NL).
insert(tree(L, W, P), E, tree(L, W, NP) ) :-
       E > W,
       insert(P, E, NP).


% stworzDrzewo( ListaElem, DrzewoBST ) - wersja z akumulatorem

stworzDrzewo(L, D) :-  stworzDrzewo(L, nil, D).

stworzDrzewo([], D, D).
stworzDrzewo([ E | L ], T, D) :-
             insert(T, E, NT),
             stworzDrzewo(L, NT, D).

% liscie(DrzewoBinarne, ListaLisci - kolejnosc od lewej), z akumulatorem

liscie(D, LL) :-  liscie(D, LL, []).

liscie(nil, L, L).
liscie(tree(nil, E, nil), [E | L], L) :-  !.
liscie(tree(L, _, P), LL, TL) :-
%	(L \= nil; P \= nil),   % zle! (wiele odpowiedzi/sukcesow)
	\+ (L=nil, P=nil),	% lub  (L,P) \= (nil,nil)
        liscie(L, LL, L1),
        liscie(P, L1, TL).

% warunek \+ (...) konieczny, wpp. sukces dla _podlisty_ lisci
% warunek postaci: (L \= nil; P \= nil) zly (wiele odpowiedzi/sukcesow)

% wersja tylko z :- ! w I klauzuli - czerwone odcięcie - zuo ;)


% wszerz( DrzewoBinarne, ListaWierzcholkowWszerz )
%   wersja prymitywna: kolejka - zwykla lista, czyli append
%   warto wspomniec: sa listy roznicowe

wszerz(D, LW) :-  wszerz_pom([D], LW).

wszerz_pom([], []).
wszerz_pom([ nil | KD ], LW) :-
           wszerz_pom( KD, LW ).
wszerz_pom([ tree( L, W, P ) | KD ], [ W | LW ]) :-
           append(KD, [ L, P ], NKD),
           wszerz_pom(NKD, LW).




% ------------------------- Grafy --------------------------------------------



% connectT(G, A, B) - DAG, reprezentacja termowa: listy krawedzi

connectT(G, A, B) :-  member(kr(A, B), G).
connectT(G, A, C) :-  member(kr(A, B), G), connectT(G, B, C).

% connectK(A, B) - DAG, reprezentacja klauzulowa

connectK(A, B) :-  edge(A, B).
connectK(A, C) :-  edge(A, B), connectK(B, C).


% path(A, B, Path) - DAG, reprezentacja klauzulowa

path(A, B, [A, B] )  :-  edge(A, B).
path(A, C, [A | P] ) :-  edge(A, B), path(B, C, P).


% pathC(A, B, P) - dowolny graf (cykliczny)
%    inicjacja listy odwiedzonych wierzcholkow na []

pathC(X, Y, P) :-  pathC(X, Y, [], P).

pathC(X, Y, _Visited, [X, Y]) :-  edge(X, Y).
pathC(X, Y,  Visited, [X | P]) :-
             edge(X, Z),
             nonmember(Z, Visited),
             pathC(Z, Y, [X | Visited], P).


% inicjacja listy odwiedzonych wierzcholkow na [A]

path2(X, Y, P) :-  path2(X, Y, [X], P).

path2(X, Y, _Visited, [X, Y]) :-  edge(X, Y).
path2(X, Y,  Visited, [X | P]) :-
             edge( X, Z ),
             nonmember(Z, Visited),
             path2(Z, Y, [Z | Visited], P).

% inicjacja listy odwiedzonych wierzcholkow na [A,B]

path3(X, Y, P) :-  path3(X, Y, [X,Y], P).

path3(X, Y, _Visited, [X, Y]) :-  edge(X, Y).
path3(X, Y,  Visited, [X | P]) :-
             edge( X, Z ),
             nonmember(Z, Visited),
             path3(Z, Y, [Z | Visited], P).

edge(a,a).
edge(a,b).
edge(a,c).
edge(b,b).
edge(b,d).
edge(c,d).
edge(d,a).



path4(X, Y, P) :-  path4(X, Y, [X], P), P=[_,_|_].

path4(X, X, _Visited, [X]).
path4(X, Y,  Visited, [X | P]) :-
             edge( X, Z ),
             nonmember(Z, Visited),
             path4(Z, Y, [Z | Visited], P).


/*
  Inicjacja listy wierzcholkow odwiedzonych.
  ------------------------------------------

Podstawowe sposoby: [], [A], [B], [A,B]

  * B -- jesli B zmienna, to nomember(X, V) poniesie porazke
         (member odniesie sukces uzgadniajac B z X),
         czyli brak odpowiedzi (nie dziala dla zapytan ze zmienna)

  * [] czy [A] ?

    jesli w grafie nie ma cykli z A do A (np. edge(a,a)), to nie ma
    wiekszej roznicy, w szczegolnosci liczba wszystkich znalezionych sciezek
    (zapytanie z samymi zmiennymi) jest taka sama w obu przypadkach

    jesli w grafie jest cykl z A do A, to wersja z [] znajduje sporo
    wiecej sciezek, zawierajacych podcykl A-A;

*/


% euler(Graf, Euler's Path)

:- op(500, xfx, <->).

euler([ W1 <-> W2 ], [ W1, W2 ]).
euler([ W1 <-> W2 ], [ W2, W1 ]).
euler(G, [ W1, W2 | Dalej ]) :-
      pobierz(G, W1 <-> W2, G1),
      euler(G1, [ W2 | Dalej ]).
      
pobierz([ W1 <-> W2 | G ], W1 <-> W2, G).
pobierz([ W1 <-> W2 | G ], W2 <-> W1, G).
pobierz([ E | G ], Kr, [ E | NG ]) :-  pobierz(G, Kr, NG).

/*
| ?- euler([a<->b, b<->c, c<->d], P).
P = [a,b,c,d] ? ;
P = [d,c,b,a] ? ;
no
*/



/* kolejka ********************************************/

% TODO difference list and queues explained

% init(Kolejka) - inicjacja kolejki (na pustą)
init(P-P).

% get(Elem, Kolejka, NowaKolejka) - pobranie
%get(Elem, [Elem|P]-K, P-K).
get(Elem, P-K, P1-K) :- \+ empty(P-K), P=[Elem|P1].

% put(Elem, Kolejka, NowaKolejka) - wstawienie
%put(Elem, P-K, P-K1) :- K=[Elem|K1].
put(Elem, P-[Elem|K1], P-K1).

% empty(Kolejka) - czy kolejka pusta
empty(P-K) :- P==K.

% init(Q), put(5,Q,Q1), put(6,Q1,Q2), get(X,Q2,Q3), write(X), put(7,Q3,Q4), put(8,Q4,Q5), get(Y,Q5,Q6), write(Y), get(Z,Q6,Q7), write(Z). 




/* wszerz ********************************************/

% wszerz( DrzewoBinarne, ListaWierzcholkowWszerz )
%   wszerz na listach różnicowych

wszerz2(D, LW) :- init(K), put(D,K,K1), wszerz2_pom(K1, LW).

wszerz2_pom(K, []) :-
	empty(K), ! .

wszerz2_pom(K, LW) :-
	\+ empty(K),
	get(D,K,K1),
	D=nil,
	!,
	wszerz2_pom(K1, LW).

wszerz2_pom(K, NLW) :-
	\+ empty(K),
	get(D,K,K1),
	D=tree(L,W,P),
	!,
	put(L,K1,K2), 
	put(P,K2,K3),
	NLW=[W|LW],
	wszerz2_pom(K3,LW).



/* wciągnięcie równości ***********/



wszerz3(D, LW) :- init(K), put(D,K,K1), wszerz3_pom(K1, LW).

wszerz3_pom(K, []) :-
	empty(K), ! .

wszerz3_pom(K, LW) :-
%	\+ empty(K),
	get(nil,K,K1),
	!,
	wszerz3_pom(K1, LW).

wszerz3_pom(K, [W|LW]) :-
%	\+ empty(K),
	get(tree(L,W,P),K,K1),
	!,
	put(L,K1,K2), 
	put(P,K2,K3),
	wszerz3_pom(K3,LW).

% odcięcia są zielone - bez nich też działa, ale tworzą się niepotrzebne punkty powrotu.
