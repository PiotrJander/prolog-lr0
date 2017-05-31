% Piotr Jander

% TODO
% TODO we'd like debug messages when things dont go as expected
% TODO is the determinism with green cuts everywhere okay?
% TODO need to detect non-LR(0) grammar
% TODO (aut, yes) albo (null, opis)

% shift(State, Terminal, NewState).
% reduce(State, Number, Nonterminal).
% goto(State, Nonterminal, NewState).

%%%%%%%%%%%%%%%
% createLR

createLR(Grammar, Automaton, Result).

% transform the grammar to convenient form
% augment with Z
% start with Z -> start production
% make state-closure
% states are arbitarily numbered
% for each production in the closure, when . Sym, add a transition on symbol

% we can be filling the table on the go
% complete the closure for the state
% for each symbol after dot
% if sym is terminal
    % see if state already exists
    % if not create one with closure
        % if dot at the end, then add a reduce
    % IMPORTANT: state is ident by the first dotted item
    % add a shift
% if sym is nt, then add a goto

% state has its number and list of dotted items, where first one is unique to the state
% first dotted item in state has the somewhere
% all other dotted items have the dot in the middle

% actually, maybe no need to store the closure anywhere
% just store states with the first ident production
% for all other productions, just add transitions
% the state machine is a graph: BFS or DFS?
% in either case, need to mark nodes white/grey/black
% could make states dynamic, or in a dictionary or list
% if dict, then check sicstus, maybe impl own dict

% dotted items: like productions prod, but with dot added
% actually distrubute rhs to prod2(NT, RHS); this can be a dotted item

% start: take start production, put a dot, make closure

% when making the automaton
% retract all relation desc automaton
% dynamically add such relations
% export those relations to a list

%%%%%%%%%%%%%%%%%%%
% example automaton

:- dynamic shift_in/3, reduce_in/3, goto_in/3.

example_automaton([
    shift_in(0, a, 3),
    shift_in(0, b, 4),
    shift_in(1, eof, accept),
    shift_in(2, a, 3),
    shift_in(2, b, 4),
    shift_in(3, a, 3),
    shift_in(3, b, 4),

    reduce_in(4, 1, a),
    reduce_in(5, 2, s),
    reduce_in(6, 2, a),

    goto_in(0, a, 2),
    goto_in(0, s, 1),
    goto_in(2, a, 5),
    goto_in(3, a, 6)
]).

shift(state(From), t(Through), state(To)) :- shift_in(From, Through, To).

reduce(state(S), SuccNumber, nt(Nonterminal)) :-
    reduce_in(S, Number, Nonterminal),
    % Number > 0, !,  % TODO trying to avoid inf loop
    succs(Number, SuccNumber).

goto(state(From), nt(Through), state(To)) :- goto_in(From, Through, To).

%%%%%%%%
% accept

assertall([]).
assertall([Term|Database]) :-
    assert(Term),
    assertall(Database).

retract_automaton :-
    retractall(shift_in(_, _, _)),
    retractall(reduce_in(_, _, _)),
    retractall(goto_in(_, _, _)).

accept(Automaton, Word) :-
    assertall(Automaton),
    automaton(Word, [state(0)]),
    retract_automaton.

automaton([], [state(accept)|_]) :- !.
automaton([Symbol|Word], [State|Stack]) :-  % we'd like the invariant that State is top of the Stack
    shift(State, t(Symbol), NewState),  % if we have a shift action
    !,
    automaton(Word, [NewState, t(Symbol), State | Stack]).
automaton(Word, [State|Stack]) :-
    reduce(State, Number, Nonterminal),  % if we have a reduce action
    !,
    stack_reduce([State|Stack], Number, Nonterminal, NewStack),  % then we reduce (and go to new state)
    automaton(Word, NewStack).

% Reduce the stack for a given production. Number is the number of symbols on the RHS.
% For each symbol on the RHS, a pair of State, Symbol is popped.
% Then the LHS nonterminal is pushed and a transition is applied.
stack_reduce([State|Stack], 0, Nonterminal, [NewState,Nonterminal,State|Stack]) :- 
    !, 
    goto(State, Nonterminal, NewState).
stack_reduce([_,_|Stack], succ(Number), Nonterminal, NewStack) :-
    !,
    stack_reduce(Stack, Number, Nonterminal, NewStack).

%%%%%%%
% tests

test_automaton :-
    example_automaton(A),
    accept(A, [b,b,eof]).

test_example :-
    example_automaton(A),

    accept(A, [b,b, eof]),
    accept(A, [a,b,b, eof]),
    accept(A, [b,a,b, eof]),
    accept(A, [a,b,a,b, eof]),
    accept(A, [a,b,a,a,b, eof]),
    accept(A, [a,a,b,a,a,b, eof]),

    \+ accept(A, [eof]),
    \+ accept(A, [a,a, eof]),
    \+ accept(A, [b, eof]),
    \+ accept(A, [a,b, eof]),
    \+ accept(A, [b,b,a, eof]),
    \+ accept(A, [b,a,b,a, eof]),
    \+ accept(A, [b,b,b, eof]),
    \+ accept(A, [a,b,a,b,a,b, eof]).

%%%%%%%%%
% helpers

% succs(0, _) :-
%     write('Trying to unify 0 with non-zero.'), nl,
%     fail.
% succs(N, _) :-
%     N < 0,
%     write('Trying with a negative number'), nl,
%     fail.
succs(0, 0).
succs(N, succ(S)) :-
    M is N - 1,
    succs(M, S).
