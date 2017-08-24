%%?@Heisenberg model exact diagonalization
% periodic BC used
clc
clear

N = 4;
Nst = 2^N;
vec_st = 0:Nst-1;
vec_st_tool = zeros(Nst,1);
vec_st_toolxor = zeros(Nst,1);
H_1 = zeros(Nst,Nst);
H_2 = zeros(Nst,Nst);
H_2 = reshape(H_2,1,Nst^2);
temp = zeros(Nst,Nst);
b = 0;
%dec2bin(vec_st);
S = sparse(Nst,Nst);
for i = 1:N
   % nearest neibor with periodic BC
    j = mod(i,N) + 1;
    %if(j==0) j = N;end
    %Sz part
    vec_st_tool = (bitget(vec_st,i)==bitget(vec_st,j));
    H_1  = H_1+ diag((2*vec_st_tool-1)*0.25);
    % S+S- part
    %2^(i-1)+2^(j-1)
    vec_st_toolxor = bitxor(vec_st,2^(i-1)+2^(j-1));
    %index =  (vec_st+1)+ (vec_st_toolxor+1-1)*Nst;
    temp = (vec_st_toolxor<Nst).*vec_st_toolxor;
    S = S + sparse(vec_st+1, temp+1, (1-vec_st_tool).*ones(1,Nst)*0.5, Nst, Nst);
    %vec_st_toolxorH_2(index) = H_2(index) + 0.5.*(1-vec_st_tool);
end
%{
for a = 1:Nst
    for i = 1:N
    j = mod(i,N) + 1;
    
    if(bitget(a-1,i)~=bitget(a-1,j))
        b = bitxor(a-1,2^(i-1)+2^(j-1))
        S = S + sparse(a,b, 0.5,Nst,Nst);
    end
    end
    
end
%}

%H_2 = reshape(H_2,Nst,Nst);
H_1 = H_1 + full(S);
[V D] = eig(H_1);
eig = diag(D);
min(eig)
%sort(diag(D))
%2*bitget(vec_st,1)-1

%H_1 = hamiltonian_1d(N,);



%{
for a  = 1 :Nst
        b = vec_st_toolxor(a)+1
        H_2(a,b) = H_2(a,b) + 0.5*(1-vec_st_tool(a))
end
   %}
