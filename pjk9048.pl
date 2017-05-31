% Piotr Jander

% TODO
% one grammar at a time
% maybe let me model the parsing table first, and do the parsing
% Automaton has states transitions
% AutomatonRun has the stack additionally
% table is static, stack changes

action(State, Terminal, NewState).
goto(State, Nonterminal, NewState).
reduce(State, Number, Nonterminal).

% so the automaton (parsing table) is essentially a set of those relations
% make it simple: hardcode those relations for our simple grammar
% acceptance by empty stack? or by accept state

% but it the reduce, we need to remember how many symbols to pop off the stack

shift_in(0, a, 3).
shift_in(0, b, 4).
% shift_in()  % TODO accept
shift_in(2, a, 3).
shift_in(2, b, 4).
shift_in(3, a, 3).
shift_in(3, b, 4).

shift(state(From), t(Through), state(To)) :- shift_in(From, Through, To).

goto_in(0, a, 2).
goto_in(0, s, 1).
goto_in(2, a, 5).
goto_in(3, a, 6).

goto(state(From), nt(Through), state(To)) :- goto_in(From, Through, To).

reduce_in(4, 1, a).
reduce_in(5, 2, s).
reduce_in(6, 2, a).

succs(0, 0).
succs(N, succ(S)) :-
    M is N - 1,
    succs(M, S).

reduce(state(S), SuccNumber, nt(Nonterminal)) :-
    reduce_in(S, Number, Nonterminal),
    succs(Number, SuccNumber).

% accept(Automaton, Word).

accept(Word) :- automaton(Word, [state(0)]).

automaton([], Stack).  % TODO
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

% TODO we'd like debug messages when things dont go as expected
% TODO accept by empty stack or by accept state

% TODO how do we accept empty string??
