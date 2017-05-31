% s(Count)  -->  ablock(Count),bblock(Count),cblock(Count). 
    
% ablock(0)  -->  []. 
% ablock(succ(Count))  -->  [a],ablock(Count). 

% bblock(0)  -->  []. 
% bblock(succ(Count))  -->  [b],bblock(Count). 

% cblock(0)  -->  []. 
% cblock(succ(Count))  -->  [c],cblock(Count).

s  -->  ablock(Count),bblock(Count),cblock(Count). 

ablock(0)  -->  []. 
ablock(NewCount)  -->  [a],ablock(Count), 
                                            {NewCount  is  Count  +  1}. 

bblock(0)  -->  []. 
bblock(NewCount)  -->  [b],bblock(Count), 
                                            {NewCount  is  Count  +  1}. 

cblock(0)  -->  []. 
cblock(NewCount)  -->  [c],cblock(Count), 
                                            {NewCount  is  Count  +  1}.

block(_A, 0)  -->  []. 
block(A, succ(Count))  -->  [A], block(A, Count).

ss(Count) --> block(a, Count), block(b, Count), block(c, Count).

% verb --> [shoots].
% verb --> [kills].

succs(0, 0).
succs(N, succ(S)) :-
    M is N - 1,
    succs(M, S).

