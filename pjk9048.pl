% Piotr Jander

% TODO
% TODO we'd like debug messages when things dont go as expected
% TODO is the determinism with green cuts everywhere okay?
% TODO need to detect non-LR(0) grammar
% TODO (aut, yes) albo (null, opis)
% TODO mark params with +/- ?
% TODO special symbols: dot, null
% TODO how to use lib lists in sicstus? maplist

:- dynamic shift_in/3, reduce_in/3, goto_in/3.
:- dynamic regularized_gramar/2.

%%%%%%%%%%%%%%%
% createLR

% createLR(Grammar, Automaton, Result) :-
%     retract_automaton,
%     retractall(gramatyka/2),
%     % invariant: only one grammar at a time    
%     assert(Grammar),
%     % TODO get all automaton relation and put into list, export as automaton

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

gramatyka(’E’,
[prod(’E’, [[nt(’E’), ’+’, nt(’T’)], [nt(’T’)]]),
prod(’T’, [[id], [’(’, nt(’E’), ’)’]]) ]).

gramatyka('S', [
    prod('S', [[nt('A'), nt('A')]]),
    prod('A', [['a', nt('A')], ['b']])
]).

regularized_gramar([
    prod(nt('S'), [[dot, nt('A'), nt('A')]]),
    prod(nt('A'), [[dot, t('a'), nt('A')], [dot, t(b)]])
]).

regularize_grammar(gramatyka(_, Productions), regularized_gramar(RegularizedProductions)) :-
    maplist(regularize_production, Productions, RegularizedProductions).
    % regularize_productions(Productions, RegularizedProductions).

% regularize_productions([], []).
% regularize_productions([Prod | Productions], [RegProd | RegularizedProductions]) :-
%     regularize_production(Prod, RegProd),
%     regularize_productions(Productions, RegularizedProductions).

regularize_production(prod(Nonterminal, RHSs), prod(nt(Nonterminal), RegRHSs)) :-
    maplist(regularize_rhs, RHSs, RegRHSs).
    % regularize_rhsides(RHSs RegRHSs).

% regularize_rhsides([], []).
% regularize_rhsides([RHS | RHSs], [RegRHS | RegRHSs]) :-
%     regularize_rhs(RHS, RegRHS),
%     regularize_rhsides(RHSs, RegRHSs).

regularize_rhs(RHS, [dot|RegRHS]) :-
    maplist(wrap_terminals, RHS, RegRHS).

wrap_terminals(nt(S), nt(S)) :- !.
wrap_terminals(S, t(S)).


% the time when we take productions from the grammar is when we complete our closure
% but then we have dot at start and terminals t(...)

add_closure_state(G, DottedRule, StateNumber, NG) :-
    % assume this state doesn't exist yet
    state_add(G, DottedRule, StateNumber, NG),

    % base case
    dotted_rule(LHS, RHS) = DottedRule,
    symbol_after_dot(RHS, null),
    length(RHS, DottedRuleLength),
    ProductionLength is DottedRuleLength - 1,
    assert(reduce_in(StateNumber, ProductionLength, RHS))
    ;
    % normal case
    % make transition from the dotted rule
    % recurse with other dotted rules
    % passing state around seems pretty unwieldy
    % use assert as we don't require backtracking

complete_closure_state(G, DottedRule, ThisStateNumber, NG) :-
    dotted_rule(LHS, RHS) = DottedRule,
    advance_dot(RHS, RHSNew, Symbol),

    % get or add the next closure state
    (
        state_has(G, dotted_rule(LHS, RHSNew), ThatStateNumber),
        NG = State
        ;
        add_closure_state(G, DottedRule, StateNumber, NG)
    ),
    (
        t(Terminal) = Symbol,
        assert(shift_in(ThisStateNumber, Terminal, ThatStateNumber))
        ;
        nt(Nonterminal) = Symbol,
        assert(goto_in(ThisStateNumber, Nonterminal, ThatStateNumber))
        % notterminal, so complete closure
        % need grammar passed
        % but grammar is static, so maybe simply assert it
    ).

%%%%%%%%%%%%%%
% Manipulating dotted rules

advance_dot([dot, Symbol | RHS], [Symbol, dot | RHS], Symbol).
advance_dot([First, Second | RHS], [First | RHSNew], Symbol) :-
    advance_dot([Second | RHS], RHSNew, Symbol).

test_advance_dot :-
    advance_dot([dot, t(a)], [t(a), dot], t(a)),
    advance_dot([t(a), dot, t(b)], [t(a), t(b), dot], t(b)),
    \+ advance_dot([t(a), dot], _, _).

symbol_after_dot([dot], null).
symbol_after_dot([dot, Symbol | RightHandSide], Symbol).
symbol_after_dot([First, Second | RightHandSide], Symbol) :-
    symbol_after_dot([Second | RightHandSide], Symbol).

test_symbol_after_dot :-
    \+ symbol_after_dot([], _),
    symbol_after_dot([dot], null),
    symbol_after_dot([dot, t(a)], t(a)),
    symbol_after_dot([t(b), dot, nt(a)], nt(a)).

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
