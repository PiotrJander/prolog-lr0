% Piotr Jander

% TODO
% TODO we'd like debug messages when things dont go as expected
% TODO is the determinism with green cuts everywhere okay?
% TODO need to detect non-LR(0) grammar
% TODO (aut, yes) albo (null, opis)
% TODO mark params with +/- ?
% TODO special symbols: dot, null

% shift(State, Terminal, NewState).
% reduce(State, Number, Nonterminal).
% goto(State, Nonterminal, NewState).

%%%%%%%%%%%%%%%
% createLR

% createLR(gramatyka(Start, Productions), Automaton, Result).

% % State stores states by number and first dotted rule
% % Assume no state for the dotted rule yet
% make_closure(FirstDottedRule, Productions, Automaton, Result, State) :-
%     % state_add(State, FirstDottedRule, StateNumber, NewState),
%     % now we've added the state
%     % time to recursively complete the closure
%     % look at Dotted Rule
%     % look at Sym after dot
%     % for each such Sym, add a transition, find or add new closure-state
%     % first: dont add transtions, just complete the dotted rules in the closure
%     % then map over this list

% make_dotted_rule(ExpandedDottedRules, NextDottedRule, EnqueuedDottedRules) :-
%     NextDottedRule = dotted_rule(_, RightHandSide),
%     symbol_after_dot(RightHandSide, Symbol),
%     % now we have a symbol

% bar(State, DottedRule, StateNumber, NewState) :-
%     % assume this state doesn't exist yet
%     state_add(State, DottedRule, StateNumber, NewState),

%     % base case
%     dotted_rule(LHS, RHS) = DottedRule,
%     symbol_after_dot(RHS, null),
%     length(RHS, DottedRuleLength),
%     ProductionLength is DottedRuleLength - 1,
%     assert(reduce_in(StateNumber, ProductionLength, RHS))
%     ;
%     % normal case
%     % make transition from the dotted rule
%     % recurse with other dotted rules
%     % passing state around seems pretty unwieldy
%     % use assert as we don't require backtracking

% baz(State, DottedRule, ThisStateNumber, NewState) :-
%     % get the transition
%     dotted_rule(LHS, RHS) = DottedRule,
%     advance_dot(RHS, RHSNew, Symbol),
%     (
%         state_has(State, dotted_rule(LHS, RHSNew), ThatStateNumber)
%         ;
%         % state not in state, call bar to create that new state
%         bar(State, DottedRule, StateNumber, NewState)
%     ),
%     (
%         t(Terminal) = Symbol,
%         assert(shift_in(ThisStateNumber, Terminal, ThatStateNumber))
%         ;
%         nt(Nonterminal) = Symbol,
%         assert(goto_in(ThisStateNumber, Nonterminal, ThatStateNumber))
%     ).

% time for some testing

advance_dot([dot, Symbol | RHS], [Symbol, dot | RHS], Symbol).
advance_dot([First, Second | RHS], [First | RHSNew], Symbol) :-
    advance_dot([Second | RHS], RHSNew, Symbol).

symbol_after_dot([dot], null).
symbol_after_dot([First, Second | RightHandSide], Symbol) :-
    First = dot,
    Second = Symbol
    ;
    symbol_after_dot([Second | RightHandSide], Symbol).

test_symbol_after_dot :-
    \+ symbol_after_dot([], _),
    symbol_after_dot([dot], null),
    symbol_after_dot([dot, t(a)], t(a)),
    symbol_after_dot([t(b), dot, nt(a)], nt(a)).

%%%%%%%%%%%%%%%%
% List of closure states

state_has([ClosureState|State], FirstDottedRule, StateNumber) :-
    ClosureState = closure_state(FirstDottedRule, StateNumber)
    ;
    state_has(State, FirstDottedRule, StateNumber).

state_add(State, FirstDottedRule, StateNumber, NewState) :-
    State = [closure_state(_, PreviousStateNumber) | _],
    StateNumber is PreviousStateNumber + 1,
    NewState = [closure_state(FirstDottedRule, StateNumber) | State]
    ;
    State = [],
    StateNumber = 0,
    NewState = [closure_state(FirstDottedRule, StateNumber)].

test_state :-
    state_add([], foo, FooNum, FooState),
    FooNum = 0,
    state_add(FooState, bar, BarNum, State),
    BarNum = 1,
    state_has(State, bar, BarNum),
    state_has(State, foo, FooNum),
    \+ state_has(State, baz, _).

%%%%%%%%%%%%%%%%%%%
% example automaton

:- dynamic shift_in/3, reduce_in/3, goto_in/3.
:- dynamic closure_state/2.

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
