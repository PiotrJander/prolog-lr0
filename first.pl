% 0. Środowisko Prologu
%    (dwa tryby pracy:  tryb zapytań i tryb definiowania,   ładowanie programu).

%    Program w Prologu = zbiór definicji relacji (relacyj wg WMT);
%    definicja pojedynczej relacji = zbiór faktów i reguł opisujących
%                       (wszystkie) przypadki, w których relacja zachodzi.
 
% 1. Dana relacja bycia dzieckiem, czyli relacja 3-argumentowa zachodząca
%    między dzieckiem, jego matką i jego ojcem.

%    Na przykład:

% dziecko(Dziecko, Matka, Ojciec)
dziecko(jasio, ewa, jan).
dziecko(stasio, ewa, jan).
dziecko(basia, anna, piotr).
dziecko(jan, ela, jakub).
dziecko(ela, hermenegilda, horacy).

% 1.a) Sformułować przeróżne zapytania dotyczące relacji bycia dzieckiem,
%      od najprostszych typu: czy Jasio jest dzieckiem Ewy i Jana, po
%      bardziej złożone (dzieci pary rodziców albo jednego rodzica itp.).

% 1.b) Zdefiniować predykaty:

ojciec(Dziecko, Ojciec) :- 
    dziecko(Dziecko, _, Ojciec).
matka(Dziecko, Matka) :- 
    dziecko(Dziecko, Matka, _).
rodzic(Dziecko, Rodzic) :- 
    ojciec(Dziecko, Rodzic); 
    % OR
    matka(Dziecko, Rodzic).
babcia(Dziecko, Babcia) :- 
    rodzic(Dziecko, X), 
    matka(X, Babcia).
wnuk(Wnuk, Dziadek) :- 
    rodzic(Wnuk, X), 
    rodzic(X, Dziadek).
przodek(Przodek, Potomek) :- 
    rodzic(Potomek, Przodek); 
    % OR                            
    rodzic(Potomek, X), 
    przodek(Przodek, X).

przodkowie(Y) :- przodek(X, Y), format('~w ~s grandparent~n', [X, "is the"]), nl, fail.

% 3. Operacje na liczbach naturalnych
%    (reprezentacja liczb: stała 0, symbol funkcyjny s/1)

%     Zdefiniować predykaty:
%      a) nat(x) wtw, gdy x jest liczbą naturalną
%      b) plus(x, y, z) wtw, gdy x + y = z
%      c) minus(x, y, z) wtw, gdy x - y = z
%      d) fib(k, n) wtw, gdy n = k-ta liczba Fibonacciego

nat(0).
nat(s(X)) :- nat(X).

plus(0, X, X).
plus(s(X), Y, Z) :- plus(X, s(Y), Z).

minus(X, Y, Z) :- plus(Y, Z, X).

fib(s(0), 0).
fib(s(s(0)), s(0)).
fib(N, s(s(K))) :- fib(X, K), fib(Y, s(K)), plus(X, Y, N).

% 4. Proste relacje na listach.
%    Zdefiniować predykaty:

% a) lista(L) wtw, gdy L jest (prologową) listą
lista([]).
lista([_|_]).

%     b) pierwszy(E, L) wtw, gdy E jest pierwszym elementem L
pierwszy([E|_], E).

%     c) ostatni(E, L) wtw, gdy E jest ostatnim elementem L

ostatni([X], X).
ostatni([_|L], X) :- ostatni(L, X).

%     d) element(E, L) wtw, gdy E jest (dowolnym) elementem L
%        (czyli member/2)

element(E, [E|_]).
element(E, [_|L]) :- element(E, L).

%     e) scal(L1, L2, L3) wtw, gdy L3 = konkatenacja listy L1 z L2
%        (czyli append/3);
%        porównać z definicją funkcji ++ w Haskellu,
%        podać (wiele) zastosowań procedury scal/3

scal([], W, W).
scal([E|L], R, [E|M]) :- scal(L, R, M).
% scal([E|L], R, W) :- scal(L, R, M), W = [E|M].

%     e') intersect(Z1,Z2) wtw, gdy zbiory (listy) Z1 i Z2 mają niepuste przecięcie

intersect([], _) :- false.
intersect(_, []) :- false.
intersect([E|X], Y) :-
    element(E, Y)
    ;
    intersect(X, Y).

%     f) podziel(Lista, NieParz, Parz) == podział danej listy na dwie
%        podlisty zawierające kolejne elementy (odpowiednio) z parzystych
%        (nieparzystych) pozycji
%        (np. podziel([1,3,5,7,9], [1,5,9], [3,7]) - sukces)

podziel([], [], []).
podziel([E], [E], []).
podziel([E,F|X], [E|Y], [F|Z]) :- podziel(X, Y, Z). 

%     g) podlista(P, L) wtw, gdy P jest spójną podlistą L

% podlista(P, L) :- append(_, P, X), append(X, _, L).

prefix([], _).
prefix([E|P], [E|L]) :- prefix(P, L).

podlista([], _).
podlista([E|P], [E|L]) :- prefix(P, L); podlista([E|P], L).
podlista([E|P], [_|L]) :- podlista([E|P], L).

%     h) podciag(P, L)  wtw, gdy P jest podciągiem L
%        (czyli niekoniecznie spójną podlistą)
%        (preferowane rozwiązanie: każdy podciąg wygenerowany jeden raz)

podciag([], _).
podciag([E|P], [E|L]) :- podciag(P, L); podciag([E|P], L).
podciag([E|P], [_|L]) :- podciag([E|P], L).

%     i) wypisz(L) == czytelne wypisanie elementów listy L, z zaznaczeniem
%        jeśli lista pusta (np. elementy oddzielane przecinkami, po
%        ostatnim elemencie kropka)

wypisz([]) :- write('[]'), nl.
wypisz([E]) :- format('~w.~n', [E]).
wypisz([E|X]) :- format('~w, ', [E]), wypisz(X).

%     j) sortowanie przez wstawianie:
%          insertionSort(Lista, Posortowana),
%          insert(Lista, Elem, NowaLista)

% insert([], Elem, [Elem]). 
% insert([E|Lista], Elem, NowaLista) :- 
%     Elem < E,
%     NowaLista = [Elem,E|Lista]
%     ;
%     Elem >= E,
%     insert(Lista, Elem, R), 
%     NowaLista = [E|R].

insert([], Elem, [Elem]).
insert([E|List], Elem, [Elem, E | List]) :- 
    Elem =< E, !.
insert([E|List], Elem, [E|Result]) :-
    Elem > E,
    insert(List, Elem, Result).

% ( 
%     Elem < E -> 
%     NowaLista = [Elem,E|Lista]
%     ; 
%     insert(Lista, Elem, R), 
%     NowaLista = [E|R]
% ).

% TODO czy wersja z akumulatorem jest lepsza?

insertionSort([], []).
insertionSort([E|Lista], Posortowana) :- insertionSort(Lista, L), insert(L, E, Posortowana).

insertionSortAcc(List, Sorted) :- insertionSortAcc(List, [], Sorted).

insertionSortAcc([], Sorted, Sorted).
insertionSortAcc([E|List], Acc, Sorted) :-
    insert(Acc, E, AccNew),
    insertionSortAcc(List, AccNew, Sorted).

l([100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,
71,70,69,68,67,66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,
41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,
11,10,9,8,7,6,5,4,3,2,1]).

% insert from the end
% sort 

%     k) zadanie domowe:
%           srodek(E, L) wtw, gdy E jest środkowym elementem L
%           (lista nieparzystej długości; np. srodek(3,[1,2,3,4,5]))
%        Uwagi:
%          - w tym zadaniu nie używamy jeszcze arytmetyki (nie trzeba)
%          - możliwe rozwiązania zarówno deklaratywne, jak i imperatywne (tylko jako
%            ćwiczenie) o dużym koszcie
%          - poszukiwane rozwiązanie o koszcie liniowym.

sameLength([], []).
sameLength([_|X], [_|Y]) :- sameLength(X, Y).

% srodek(3, X) only yields a singleton list
srodek(E, [E]).
srodek(E, L) :-
    sameLength(X, Y),
    append(X, [E|Y], Z),
    sameLength(Z, L),
    !,
    Z = L.

%%%%%%%%%%%
% Lab 2
%%%%%%%%%%%

% Zdefiniować predykaty:

% a) suma(L, S) wtw, gdy S = suma elementów listy L

suma(L, S) :- suma(L, 0, S).

suma([], A, A).
suma([E|L], A, S) :- B is A + E, suma(L, B, S).

% make list

% range(X, X, []).
% range(X, Y, L) :-
%     Z is X + 1,
%     range(Z, Y, K),
%     L = [X|K].

% % TODO does this even make sense?
% series(X, [X|L]) :-
%     Y is X + 1,
%     series(Y, L).

% take(0, [], _).
% take(N, [E|S], [E|L]) :-
%     M is N - 1,
%     take(M, S, L).

% TODO should we cut for all pattern matching?
zrobliste(0,[]):- !.
zrobliste(A,[A|L]) :- B is A-1, zrobliste(B,L).

% b) dlugosc(L, K) wtw, gdy K = liczba elementów listy L (length/2)

dlugosc(L, K) :- dlugosc(L, 0, K). 

dlugosc([], K, K).
dlugosc([_|L], A, K) :-
    B is A + 1,
    dlugosc(L, B, K).

% c) min(L, M) wtw, gdy M jest minimalnym elementem L (L = lista np. liczb całkowitych)

min([E|L], M) :- min(L, E, M).

% TODO for 100000000 length list -> stack overflow
min([], M, M).
min([E|L], N, M) :-
    E < N, !,
    min(L, E, M)
    ;
    E >= N,
    min(L, N, M).

% d) odwroc(L, R) wtw, gdy R jest odwróconą listą L (np. odwroc([1,2,3,4], [4,3,2,1]) - sukces)

odwroc(L, R) :- odwroc(L, [], R).

odwroc([], R, R).
odwroc([E|L], H, R) :- odwroc(L, [E|H], R).

% e) palindrom(Slowo) wtw, gdy (lista) Slowo jest palindromem (np. palindrom([k,a,j,a,k]), palindrom([1,b,2,b,1]) - sukcesy)

palindrom(L) :- odwroc(L, L).

pal(L) :-  pal(L, []).

pal(L, L).                  % parzysta dlugosc
pal([_|L], L).              % nieparzysta, np. kajak, czyli [k,a,j,a,k]
pal([E|L], A) :-  pal(L, [E|A]).

% f) slowo(Slowo) == Slowo= anbnanbn (Uwaga: bez arytmetyki!) dwa warianty: (*) n > 0 (**) n >= 0 (np. slowo([a,a,b,b]) - sukces)

% TODO why get
% ?- replicate(a, X, [a, a, a]).
% ERROR: Arguments are not sufficiently instantiated

replicate(_, 0, []).
replicate(X, N, [X|L]) :-
    M is N - 1,
    replicate(X, M, L).

ab([]).
ab(L) :-
    ab(K),
    append([a|K], [b], L).

slowo(Slowo) :-
    % replicate(a, N, A),
    % replicate(b, N, B),
    % append(A, B, Half),
    ab(Half),
    append(Half, Half, Slowo).

% g) slowo(Zdanie, Reszta) == Zdanie = Slowo * Reszta, Slowo - jw. (np. slowo([a,a,b,b,c,d], [c,d]) - sukces)

% TODO dlaczego slowo jest zawsze puste?
kstar([]).
kstar(Slowo) :- slowo(X), kstar(Y), append(X, Y, Slowo).

% TODO how to bind variable in SWIPL REPL?

% h) flagaPolska(Lista, Flaga) wtw, gdy Flaga jest posortowaną listą Lista, złożoną ze stałych b,c
% (np. flagaPolska([b,c,b,c], [b,b,c,c]) - sukces)

% idea: keep two stacks for bs and cs
flagaPolska(R, S) :- flagaPolska(R, [], [], S).

flagaPolska([], B, C, S) :- append(B, C, S).
flagaPolska([b|R], B, C, S) :- flagaPolska(R, [b|B], C, S).
flagaPolska([c|R], B, C, S) :- flagaPolska(R, B, [c|C], S).

% i) ew. flagaHolenderska(ListaRWB, RWB) (flaga: red-white-blue)

flagaHolenderska(ListaRWB, RWB) :- flagaHolenderska(ListaRWB, [], [], [], RWB).

flagaHolenderska([], R, W, B, RWB) :- append(R, W, P), append(P, B, RWB).
flagaHolenderska([r|ListaRWB], R, W, B, RWB) :- flagaHolenderska(ListaRWB, [r|R], W, B, RWB).
flagaHolenderska([w|ListaRWB], R, W, B, RWB) :- flagaHolenderska(ListaRWB, R, [w|W], B, RWB).
flagaHolenderska([b|ListaRWB], R, W, B, RWB) :- flagaHolenderska(ListaRWB, R, W, [b|B], RWB).

% j) quickSort(L, S) wtw, gdy S jest wynikiem sortowania L (algorytm QuickSort)

% wersja bez akumulatora

% TODO why doesn't work??
% partition(_Pivot, _Smaller, _Greater, []).
partition(_, [], [], []).
partition(Pivot, Smaller, Greater, [E|List]) :-
    E > Pivot,
    !,
    partition(Pivot, Smaller, [E|Greater], List)
    ;
    E < Pivot,
    !,
    partition(Pivot, [E|Smaller], Greater, List)
    ;
    E = Pivot,
    partition(Pivot, Smaller, Greater, List).

% quickSort([], []).
quickSort([E], [E]).
quickSort([E|L], S) :-
    partition(E, Smaller, Greater, [E|L]),
    quickSort(Smaller, SmallerSorted),
    quickSort(Greater, GreaterSorted),
    append(SmallerSorted, [E|GreaterSorted], S).

% wersja z akumulatorem (czyli bez append)

% k) flatten(L, F) wtw, gdy L jest zagnieżdżoną listą list, których elementami są liczby całkowite, a F jest spłaszczoną listą L (np. flatten([1,[[[[2,[3]]], 4], 5]], [1,2,3,4,5]) - sukces)

flatten([], []).
flatten([E|L], F) :-
    is_list(E),
    !,
    flatten(E, EFlattened),
    flatten(L, LFlattened),
    append(EFlattened, LFlattened, F)
    ;
    flatten(L, LFlattened),
    F = [E|LFlattened].

% flatten(L, F) :- flatten(L, [], F).

% flatten([], A, A).
% flatten([E|L], A, F) :-
%     is_list(E),
%     flatten(E, EFlattened),
%     flatten(L, LFlattened),
%     append(EFlattened, LFlattened, )
%     ;
%     E = [F|Rest],
%     flatten

warm_blooded(penguin).
warm_blooded(human).
 
produce_milk(penguin).
produce_milk(human).
 
have_feathers(penguin).
have_hair(human).
 
mammal(X) :-
  warm_blooded(X),
  produce_milk(X),
  have_hair(X).

%%%%%%%%%%%
% Lab 3
%%%%%%%%%%%

% a) drzewo(D) wtw, gdy D jest drzewem binarnym

drzewo(nil).
drzewo(node(left, _, right)) :- drzewo(left), drzewo(right).

% b) insertBST(DrzewoBST, Elem, NoweDrzewoBST)

insertBST(nil, Elem, node(nil, Elem, nil)).
insertBST(node(Left, Value, Right), Elem, NoweDrzewoBST) :-
    Elem =< Value, !,
    insertBST(Left, Elem, NewLeft),
    NoweDrzewoBST = node(NewLeft, Value, Right)
    ;
    Elem > Value,
    insertBST(Right, Elem, NewRight),
    NoweDrzewoBST = node(Left, Value, NewRight).

drzewko(Drzewko) :-
    insertBST(nil, 5, T1),
    insertBST(T1, 7, T2),
    insertBST(T2, 6, T3),
    insertBST(T3, 2, T4),
    insertBST(T4, 4, Drzewko).

% 3 wersje tej procedury: jeśli element już jest w drzewie, to sukces (ponowne wstawienie elementu lub nie) lub porażka

% b') (ew. jeśli czas pozwoli albo zadanie do domu) deleteBST/3

deleteBST(node(Left, Value, Right), Elem, Tree) :-
    Elem < Value, !,
    deleteBST(Left, Elem, NewLeft),
    Tree = node(NewLeft, Value, Right)
    ;
    Elem > Value, !,
    deleteBST(Right, Elem, NewRight),
    Tree = node(Left, Value, NewRight)
    ;
    Elem = Value,
    (
        (Left, Right) = (nil, nil), !,
        Tree = nil
        ;
        Left = nil, !,
        Tree = Right
        ;
        Right = nil, !,
        Tree = Left
        ;
        minimumBST(Right, Min),
        deleteBST(Right, Min, NewRight),
        Tree = node(Left, Min, NewRight)
    ).

minimumBST(node(Left, _, _), Min) :- minimumBST(Left, Min).
minimumBST(node(nil, Value, _), Value). 

% c) (ew. wypiszBST(D) = wypisanie na ekranie, porządek infiksowy)

wypiszBST(nil).
wypiszBST(node(Left, Value, Right)) :-
    wypiszBST(Left),
    write(Value),
    nl,
    wypiszBST(Right).

% d) wypiszBST(D, L) wtw, gdy L=lista wszystkich wierzchołków D (porządek infiksowy)

listBST(nil, []).
listBST(node(Left, Value, Right), Result) :-
    listBST(Left, LeftList),
    listBST(Right, RightList),
    append(LeftList, [Value|RightList], Result).

% e) stworzBST(L, D) wtw, gdy D jest drzewem BST zawierającym wszystkie elementy listy L (akumulator, ew. bez)

% TODO ogolnie, jaki sens akumulatora? ogonowa rekurencja? wyw rek na koncu?

% stworzBST([], nil).
% stworzBST([E|L], D) :-
%     stworzBST(L, D2),
%     insertBST(D2, E, D).

stworzBST(L, D) :- stworzBST(L, nil, D).

stworzBST([], D, D).
stworzBST([E|L], A, D) :-
    insertBST(A, E, ANew),
    stworzBST(L, ANew, D).

% f) liscie(D, L) wtw, gdy L = lista wszystkich liści, od lewej do prawej

liscie(nil, []).
liscie(node(Left, Value, Right), L) :-
    (Left, Right) = (nil, nil), !,
    L = [Value]
    ;
    % (Left, Right) /= (nil, nil), % TODO should we duplicate condition here?
    liscie(Left, LeftLeaves),
    liscie(Right, RightLeaves),
    append(LeftLeaves, RightLeaves, L).

% g) sortBST(L, S) wtw, gdy S = lista posortowana, przy użyciu drzew BST

sortBST(L, S) :-
    stworzBST(L, D),
    listBST(D, S).

% h) wszerz(D, L) wtw, gdy L = lista wszystkich wierzchołków wszerz

% TODO difference lists
% to potrzeba list roznicowych
% kolejka.

wszerz(Tree, List) :- wszerz([Tree], [], List).

% TODO how come we can have pattern matching on both sides?
% TODO why we get the reverse?

wszerz([], Acc, Acc).
wszerz([nil|Queue], Acc, List) :- wszerz(Queue, Acc, List).
wszerz([node(Left, Value, Right)|Queue], Acc, List) :-
    append(Queue, [Left, Right], NewQueue),
    wszerz(NewQueue, [Value|Acc], List).

% wszerz(D, LW) :-  wszerz_pom([D], LW).

% wszerz_pom([], []).
% wszerz_pom([ nil | KD ], LW) :-
%            wszerz_pom( KD, LW ).
% wszerz_pom([ tree( L, W, P ) | KD ], [ W | LW ]) :-
%            append(KD, [ L, P ], NKD),
%            wszerz_pom(NKD, LW).

% Grafy
% (skierowane, nieskierowane, cykliczne, acykliczne, etykietowane itd.) (reprezentacja termowa oraz reprezentacja klauzulowa)

% Reprezentacja grafów: zbiory krawędzi, czyli

% listy termów postaci np. kr(A, B)
% klauzule unarne typu: edge(A, B).
% Grafy DAG (skierowane, acykliczne)
% Zdefiniować predykaty:

% a) connect(A,B), connect(Graf,A,B) wtw, gdy istnieje ścieżka z A do B.
% Uwaga: ścieżka = niepusty (!) ciąg krawędzi

connect(A, B) :- kr(A, B).
connect(A, B) :- kr(A, C), connect(C, B).

krs(a, b).
krs(b, c).
krs(c, d).
krs(d, e).

kr(a, b).
kr(b, c).
kr(c, a).
kr(c, d).
kr(d, e).

% TODO how to avoid infinite recursion here?
% kr(A, B) :- kr(B, A).

graph([krs(a, b), krs(b, c), krs(c,d), krs(d, e)]).

connect(Graph, A, B) :- member(krs(A, B), Graph).
connect(Graph, A, B) :- 
    member(krs(A, C), Graph),
    connect(Graph, C, B).

% b) path(A,B,P) wtw, gdy P = opis ścieżki z A do B, tzn. P = [A, ..., B]

path(A, B, [krs(A, B)]) :- krs(A, B).
path(A, B, [krs(A, C)|P]) :-
    krs(A, C),
    path(C, B, P).

% Grafy (nie)skierowane, (a)cykliczne
% c) pathC(A,B,P) w dowolnym grafie skierowanym (cyklicznym)

pathC(A, B, P) :- pathC(A, B, [], P).

pathC(A, B, _Visited, [A, B]) :- kr(A, B).
pathC(A, B, Visited, [A|P]) :-
    kr(A, C),
    \+ member(C, Visited),
    pathC(C, B, P).

% d) euler/? - czy dany graf jest grafem Eulera, czyli znalezienie (sprawdzenie) ścieżki Eulera (wprost z definicji):
% ścieżka, która przechodzi przez każdą krawędź grafu dokładnie raz

% Struktury otwarte i różnicowe
% Implementacja kolejki FIFO, czyli:

% a) init(Kolejka) - inicjalizacja kolejki (na pustą)

% b) get(Elem, Kolejka, NowaKolejka) - pobranie

% c) put(Elem, Kolejka, NowaKolejka) - wstawienie

% d) empty(Kolejka) - czy kolejka pusta

% wszerz(DrzewoBinarne, ListaWierzchWszerz)

inc(X, Y) :-
    Y is X + 1.