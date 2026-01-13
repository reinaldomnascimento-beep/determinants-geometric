% -------------------------------------------------------------------------
% Determinants by Geometric Menorization
%
% This script implements a recursive construction of determinants based on
% geometric menorization, emphasizing the role of pivots in the formation
% of lower-order minors.
%
% The original algorithm and implementation date back to the period
% 2005-2006 and are preserved in full below.
%
% Additional documentation, conceptual framing, and didactic visualization
% were incorporated in 2026 in the context of the project:
%   "Determinants by Geometric Menorization"
%
% Author (original concept and code):
%   Reinaldo Mauricio do Nascimento
%
% Technical and conceptual collaboration (2026):
%   ChatGPT - OpenAI
%
% Related resources:
%   * Article (PDF - English)
%   * Artigo (PDF - Portugues)
%   * Interactive Web Visualization (HTML/JS)
%   * README
%
% Date of this header insertion: 2026-01
% -------------------------------------------------------------------------

function d=detp(A,op)
%                     DETP Linear algebra with pivoting
%========================================================================                     
%1 - detp(A) ==> Determinant calculus
%   [Aij]     ==> [Bij]         ==>  ...  ==>  [Cij]     ==> [P1]
%        n x n         n-1 x n-1                    2 x 2        1 x 1
%                                      n          
%                                      __         
%                             det(A) = || Pk ^ (2-k)
%             n                        k=1          
%             __         
%   detp(A) = || pk ==>  pn  =   Pn
%             k=1     p(n-1) = P(n-1)/Pn
%                        ... =   ...
%                        p2  = P2/(Pn*P(n-1)*...*P3)
%                        p1  = P1/(Pn*P(n-1)*...*P3*P2)
%========================================================================
%2 - detp(A,B) ==> Linear equations system B=AX
%========================================================================
%3 - detp(A,'p') ==> Equation of the "plane" using determinant(det(A)~=0)
%                           n  n
%    "Plane" [Aij]     ==> SumSumijXj=det(A); ij=cofactor.
%                 n x n    j=1i=1
%    detp(A,'p') ==> Linear equations system B=AX; B=[det(A)]
%                                                           n x 1
%==========================================================================
%4 - detp(A,-1) ==> Matrix inverse
%==========================================================================
%5 - detp(A,'ulp') ==> ULP factorization. Returns upper triangular matrix U,
%    lower triangular matrix L and permutation matrix P so that U*L = A*P. 
%==========================================================================    

%January,2008 by Reinaldo M. do Nascimento
%Revision: 0.1:- Code and help - 2005/07/08  18:00
%Revision: 0.2:- Code, help and new functions - 2008/02/05 19:30
%Revision: 0.3:- Code & ULP factorization - 2008/02/11 15:55

if diff(size(A)),error('Matrix A must be square.');else,n=length(A);end
if nargin==2 & op==-1
             op='inv';
             disp('                          inv (A)');
elseif nargin==2 & sym(op)=='p'
             B(n,1)=0;
             disp('            Equation of the "plane" using determinant');
elseif nargin==2 & size(op,1)==n & size(op,2)==1
             B=op;op='B=AX';
             disp('                 Linear equations system B=AX');
elseif nargin==2 & sym(op)=='ulp'
             B(n,1)=0;
             disp('                ULP factorization ==> U*L = A*P');
elseif nargin==2
       disp('Options==>detp(A), detp(A,-1), detp(A,B), detp(A,''p'') and detp(A,''ulp'').');
       error('Matrix dimensions must agree ==> B(length(A),1).');
else
             op='detA';
             disp('                          det (A)'); 
end

global LiFo
detA=1;LiFo=[];
%Check Pivot A(n,n), interchanges and det(A)====================
while n>1
       for k=n:-1:1
            if A(n,k)~=0
                 detA=detA*A(n,k);
                 if k~=n
                     A(:,[k n])=A(:,[n k]);
                     LiFo=[k n;LiFo];
                     detA=-detA;
                 end
            break
            elseif k==1
                    detA=0;
                    if sym(op)=='detA'
                              n=1;
                    elseif sym(op)=='ulp'
                                  B(n)=1;
                    else
                        A(n,n)=1;
                    end
            end
       end
%==========A======================a==========Pivoting (A to a)
%[A11 ... A1j ... A1n]  [a11 ... a1j ... A1n]
%[... ... ... ... ...]  [... ... ... ... ...] 
%[Ai1 ... Aij ... Ain]=>[ai1 ... aij ... Ain]=>aij=[Aij*Ann-Anj*Ain]/Ann
%[... ... ... ... ...]  [... ... ... ... ...] 
%[An1 ... Anj ... Ann]  [An1 ... Anj ... Ann]              
       for j=1:k-1                           
            if A(n,j)~=0
                q=A(n,j)/A(n,n);
                 for i=1:n-1
                      A(i,j)=A(i,j)-q*A(i,n);
                 end
            end
       end
       n=n-1;
end
detA=detA*A(1,1);
if detA==0 & sym(op)~='ulp' & sym(op)~='detA'
       disp('Warning: Matrix is singular to working precision.');
       if A(1,1)==0,A(1,1)=1;end
end
%======================================================Switch Case
switch op
    case 'detA'
        d=detA;
    case 'B=AX'
        d=B_to_b(A,B);
    case 'inv'
        d=Inverse(A);
    case 'p'
        B(:)=detA;
        d=struct('vector',B_to_b(A,B),'det',detA);
    case 'ulp'
        [U,L,P]=ulp(A,B);
        d=struct('U',U,'L',L,'P',P);
end
%Pivoting (B to b)========B=====a=========b==========================
%                       [ B1] [a1n]     [ b1] = [B1*ann-Bn*a1n]/ann
%                       [...] [...]     [...]
%                       [ Bi] [ain] ==> [ bi] = [Bi*ann-Bn*ain]/ann
%                       [...] [...]     [...]
%                       [ Bn] [ann]     [ bn] = [Bn]
function d=B_to_b(A,B)
%
for n=length(A):-1:2;
     if B(n)~=0
         q=B(n)/A(n,n);
          for i=1:n-1
               B(i)=B(i)-q*A(i,n);
          end
     end
end
d=b_to_X(A,B);
%================================================Back Substitution
%                              i-1
%                   X(i)=[b(i)-Sum(a(i,j)*X(j))]/p(i); p(i)=a(i,i)
%                              j=1
function d=b_to_X(A,B)
%
global LiFo
for i=1:length(A)
     X(i,1)=B(i);
           for j=1:i-1
                X(i)=X(i)-A(i,j)*X(j);
           end
     X(i)=X(i)/A(i,i);
end
%If necessary, reorder X=======================
if LiFo
    for i=1:size(LiFo,1)
         X([LiFo(i,1:2)])=X([LiFo(i,2:-1:1)]);
    end
end
d=X;
%=======================================================inv(A)
function d=Inverse(A)
%
d=[];
for n=length(A):-1:1
     B(n,1)=1;
     d=[B_to_b(A,B),d];
     B(n)=0;
end
%ULP factorization=============================
function [U,L,P]=ulp(A,B)
%
global LiFo
n=length(A);
%================================Upper & Lower triangular matrix
for j=n:-1:1
     L(j,j)=A(j,j);U(j,j)=1;pk=1;
     if A(j,j)~=0
          pk=A(j,j);
     elseif B(j)==1 | j==1%Optional(j==1)
              L(j,j)=1;U(j,j)=0;
     end
     for i=1:n
          if i>j
              L(i,j)=A(i,j);
          end
          if i<j
              U(i,j)=A(i,j)/pk;
          end
     end
end
%Permutation matrix======================================
P=eye(n,n);
if LiFo
    for i=size(LiFo,1):-1:1
         P(:,[LiFo(i,1:2)])=P(:,[LiFo(i,2:-1:1)]);
    end
end

%reinaldomn@terra.com.br
%Electrical Technician - CEFET-MG
%Electronic Technician - COTEMIG-MG
