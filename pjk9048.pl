% Piotr Jander

% TODO
% TODO we'd like debug messages when things dont go as expected
% TODO is the determinism with green cuts everywhere okay?
% TODO need to detect non-LR(0) grammar
% TODO (aut, yes) albo (null, opis)
% TODO maybe S → E eof
% TODO am I sure a state is characterized by just one dotted item?

:- dynamic shift_in/3, reduce_in/3, goto_in/3.
:- dynamic regularized_gramar/2.

%%%%%%%%%%%%%%%
% createLR

% createLR(gramatyka(StartSymbol, Productions), Automaton, Result) :-
%     retract_automaton,
%     retractall(gramatyka/2),
%     regularize_productions(Productions, RegularizedProductions),
%     assert(regularized_gramar(RegularizedProductions)),
%     add_closure_state([], dotted_rule(nt('Z'), [nt(StartSymbol)]), _, _).
%     % at this point we're ready to gather all table relations

%%%%%%%%%%%%%%%%
% construct automaton

add_closure_state(G, DottedRule, StateNumber, NG) :-
    % assume this state doesn't exist yet
    state_add(G, DottedRule, StateNumber, NG),

    % base case
    dotted_rule(_, RHS) = DottedRule,
    (
        % dot at the end of production, reduce
        dot_at_end_of_production(RHS),
        length(RHS, DottedRuleLength),
        ProductionLength is DottedRuleLength - 1,
        assert(reduce_in(StateNumber, ProductionLength, RHS))
        ;
        complete_closure_state(G, DottedRule, StateNumber, NG)
    ).

test_automaton_construction(NG) :-
    add_closure_state([], dotted_rule(nt('Z'), [nt('S')]), _, NG).

complete_closure_state(G, dotted_rule(LHS, RHS), ThisStateNumber, NG) :-
    advance_dot(RHS, RHSNew, Symbol),

    % get or add the next closure state
    (
        state_has(G, dotted_rule(LHS, RHSNew), ThatStateNumber),
        IntermediateG = G
        ;
        add_closure_state(G, dotted_rule(LHS, RHSNew), ThatStateNumber, IntermediateG)
    ),
    (
        t(Terminal) = Symbol,
        assert(shift_in(ThisStateNumber, Terminal, ThatStateNumber)),
        NG = IntermediateG
        ;
        nt(Nonterminal) = Symbol,
        assert(goto_in(ThisStateNumber, Nonterminal, ThatStateNumber)),

        % nonterminal, so keep completing the closure
        productions_for_nonterminal(nt(Nonterminal), NonterminalProductions),
        complete_closure_state1(IntermediateG, nt(Nonterminal), NonterminalProductions, ThisStateNumber, NG)
    ).

complete_closure_state1(G, _, [], _, G).
complete_closure_state1(G, Nonterminal, [NontProd|NonterminalProductions], ThisStateNumber, NG) :-
    complete_closure_state(G, dotted_rule(Nonterminal, NontProd), ThisStateNumber, IntermediateG),
    complete_closure_state1(IntermediateG, Nonterminal, NonterminalProductions, ThisStateNumber, NG).

%%%%%%%%%%%%
% Regularize grammar

gramatyka('S', [
    prod('S', [[nt('A'), nt('A')]]),
    prod('A', [['a', nt('A')], ['b']])
]).

regularized_grammar([
    prod(nt('S'), [[dot, nt('A'), nt('A')]]),
    prod(nt('A'), [[dot, t('a'), nt('A')], [dot, t(b)]])
]).

regularize_productions(Productions, RegularizedProductions) :-
    maplist(regularize_production, Productions, RegularizedProductions).

regularize_production(prod(Nonterminal, RHSs), prod(nt(Nonterminal), RegRHSs)) :-
    maplist(regularize_rhs, RHSs, RegRHSs).

regularize_rhs(RHS, [dot|RegRHS]) :-
    maplist(wrap_terminals, RHS, RegRHS).

wrap_terminals(nt(S), nt(S)) :- !.
wrap_terminals(S, t(S)).

test_regularize_grammar :-
    gramatyka('S', Productions),
    regularized_grammar(RegularizedProductions),
    regularize_productions(Productions, RegularizedProductions).

%%%%%%%%%%%%
% Productions for nonterminal

productions_for_nonterminal(Nonterminal, NonterminalProductions) :-
    regularized_grammar(Productions),
    productions_for_nonterminal1(Nonterminal, Productions, NonterminalProductions).

productions_for_nonterminal1(Nonterminal, [NonProd|Productions], NonterminalProductions) :-
    prod(Nonterminal, NonterminalProductions) = NonProd
    ;
    productions_for_nonterminal1(Nonterminal, Productions, NonterminalProductions).

test_productions_for_nonterminal :-
    productions_for_nonterminal(nt('S'), [[dot, nt('A'), nt('A')]]),
    productions_for_nonterminal(nt('A'), [[dot, t('a'), nt('A')], [dot, t(b)]]).

%%%%%%%%%%%%%%
% Manipulating dotted rules

advance_dot([dot], [dot], null).
advance_dot([dot, Symbol | RHS], [Symbol, dot | RHS], Symbol).
advance_dot([First, Second | RHS], [First | RHSNew], Symbol) :-
    advance_dot([Second | RHS], RHSNew, Symbol).

test_advance_dot :-
    advance_dot([dot], [dot], null),
    advance_dot([dot, t(a)], [t(a), dot], t(a)),
    advance_dot([t(a), dot, t(b)], [t(a), t(b), dot], t(b)),
    \+ advance_dot([t(a), dot], _, _).

dot_at_end_of_production([dot]).
dot_at_end_of_production([dot|_]) :- false.
dot_at_end_of_production([_|RHS]) :- dot_at_end_of_production(RHS).

test_dot_at_end_of_production :-
    dot_at_end_of_production([dot]),
    dot_at_end_of_production([nt(a), t(b), dot]),
    \+ dot_at_end_of_production([nt(a), dot, t(b)]).

%%%%%%%%%%%%%%%%
% List of closure states

state_has([ClosureState|G], FirstDottedRule, StateNumber) :-
    ClosureState = closure_state(FirstDottedRule, StateNumber)
    ;
    state_has(G, FirstDottedRule, StateNumber).

state_add(G, FirstDottedRule, StateNumber, NG) :-
    G = [closure_state(_, PreviousStateNumber) | _],
    StateNumber is PreviousStateNumber + 1,
    NG = [closure_state(FirstDottedRule, StateNumber) | G]
    ;
    G = [],
    StateNumber = 0,
    NG = [closure_state(FirstDottedRule, StateNumber)].

test_state :-
    state_add([], foo, FooNum, FooG),
    FooNum = 0,
    state_add(FooG, bar, BarNum, G),
    BarNum = 1,
    state_has(G, bar, BarNum),
    state_has(G, foo, FooNum),
    \+ state_has(G, baz, _).

%%%%%%%%%%%%%%%%%%%
% example automaton

% shift(State, Terminal, NewState).
% reduce(State, Number, Nonterminal).
% goto(State, Nonterminal, NewState).

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
