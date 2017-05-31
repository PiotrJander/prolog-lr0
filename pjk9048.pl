% Piotr Jander

% TODO
% TODO we'd like debug messages when things dont go as expected

% shift(State, Terminal, NewState).
% reduce(State, Number, Nonterminal).
% goto(State, Nonterminal, NewState).

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
    goto_in(3, a, 6).
]).

shift_in(0, a, 3).
shift_in(0, b, 4).
shift_in(1, eof, accept).
shift_in(2, a, 3).
shift_in(2, b, 4).
shift_in(3, a, 3).
shift_in(3, b, 4).

shift(state(From), t(Through), state(To)) :- shift_in(From, Through, To).

reduce_in(4, 1, a).
reduce_in(5, 2, s).
reduce_in(6, 2, a).

reduce(state(S), SuccNumber, nt(Nonterminal)) :-
    reduce_in(S, Number, Nonterminal),
    succs(Number, SuccNumber).

goto_in(0, a, 2).
goto_in(0, s, 1).
goto_in(2, a, 5).
goto_in(3, a, 6).

goto(state(From), nt(Through), state(To)) :- goto_in(From, Through, To).

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

automaton([], [state(accept)]).
automaton([Terminal|Word], [State|Stack]) :-  % we'd like the invariant that State is top of the Stack
    shift(State, Terminal, NewState),  % if we have a shift action
    automaton(Word, [NewState, Terminal, State | Stack])
    ;
    reduce(State, Number, Nonterminal),  % if we have a reduce action
    stack_reduce([State|Stack], Number, Nonterminal, NewStack),  % then we reduce (and go to new state)
    automaton([Terminal|Word], NewStack).

% Reduce the stack for a given production. Number is the number of symbols on the RHS.
% For each symbol on the RHS, a pair of State, Symbol is popped.
% Then the LHS nonterminal is pushed and a transition is applied.
stack_reduce([State|Stack], 0, Nonterminal, [NewState,Nonterminal,State|Stack]) :- goto(State, Nonterminal, NewState).
stack_reduce([_,_|Stack], succ(Number), Nonterminal, Stack).

%%%%%%%%%
% helpers

succs(0, 0).
succs(N, succ(S)) :-
    M is N - 1,
    succs(M, S).