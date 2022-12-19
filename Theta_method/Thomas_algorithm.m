function sol = Thomas_algorithm(a,b,c,d)
%b is the left diagonal array
%a is the center diagonal array
%c is the right diagonal array
% d is the answer array

N = length(d);

cn(1)=c(1)/b(1);
dn(1)=d(1)/b(1);

if all(abs(b)>=abs(a)+abs(c)) & any(abs(b)>abs(a)+abs(c))
     'Strict condition of main diagonal dominance satisfied.';
elseif all(abs(b)>=abs(a)+abs(c))
    'non-strict condition of main diagonal dominance satisfied.';
else
    disp ('Condition of main diagonal dominance NOT satisfied.')
    return
end

for i=2:N
    cn(i)=c(i)/(b(i)-a(i)*cn(i-1));
    dn(i)=(d(i)-a(i)*dn(i-1))/(b(i)-a(i)*cn(i-1));
end

sol(N) = dn(N);

for i=(N-1):-1:1;
    sol(i)=dn(i)-(cn(i)*sol(i+1));
end