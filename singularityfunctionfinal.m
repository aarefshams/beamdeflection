L=(input('L:'));               % Beam length
E =(input('E:'));              % Young's modulus
I =(input('I:'));              % Moment of inertia
bx1=(input('bx1:'));           %مکان تکیه گاه اول
bv1=(input('bv1:'));           %مقدار خیز در تکیه گاه اول
bx2=(input('bx2:'));           %مکان تکیه گاه دوم
bv2=(input('bv2:'));           %مقدار خیز در تکیه گاه دوم
bdv1=(input('bdv1:'));         %شیب خیز در تکیه گاه اول
bdv2=(input('bdv2:'));         %شیب خیز در تکیه گاه دوم
x=0:0.001:L;

%for P
P=(input('P:'));
ap=(input('aP:'));
dvp=@(p,a,x) ~isempty(p)*(p/2.*sing(x,a,2));
vp= @(p,a,x) ~isempty(p)*(p/6.*sing(x,a,3));
%(vp(P(i),ap(i),x));

%for M
M=(input('M:'));
aM=(input('aM:'));
dvM=@(M,a,x) ~isempty(M)*(M.*sing(x,a,1));
vM=@(M,a,x) ~isempty(M)*((M/2).*sing(x,a,2));
%(vM(M(i),aM(i),x));

%for w Distributed load
w=(input('w:'));
aw=(input('aw:'));
dvw=@(w,a,x) ~isempty(w)*((w/6).*sing(x,a,3));
vw=@(w,a,x) ~isempty(w)*((w/24).*sing(x,a,4));
%(vw(w(i),aw(i),x));

%for w Triangular Distributed load
w_=(input('w_:'));
aw_=(input('aw_:'));
b=L-aw_;
dvw_=@(w_,a,x) ~isempty(w_)*((w_/b*24).*sing(x,a,4));
vw_=@(w_,a,x) ~isempty(w_)*((w_/b*120).*sing(x,a,5));
%(vw_(w_(i),aw_(i),x));


if bdv1 == 0
    C1=bdv1-(  sum(dvp(P,ap,bx1))+  sum(dvM(M,aM,bx1))  +sum(dvw(w,aw,bx1))) + sum(dvw_(w_,aw_,bx1));
    C2=bv1-(  sum(vp(P,ap,bx1))+  sum(vM(M,aM,bx1))  +sum(vw(w,aw,bx1))) + sum(vw_(w_,aw_,bx1));
else 
    syms C1 C2 ;
    eqn1 = ((  sum(vp(P,ap,bx1))+  sum(vM(M,aM,bx1))  +sum(vw(w,aw,bx1))) + sum(vw_(w_,aw_,bx1))+C1*bx1+C2)==bv1;
    eqn2 = ((  sum(vp(P,ap,bx2))+  sum(vM(M,aM,bx2))  +sum(vw(w,aw,bx2))) + sum(vw_(w_,aw_,bx2))+C1*bx2+C2)==bv2;
    sol = solve([eqn1, eqn2], [C1,C2]);
    C1 = sol.C1;
    C2 = sol.C2;
    
end


V=0;
for i = 1:numel(P)
    V=V+(vp(P(i),ap(i),x));
end
for i = 1:numel(M)
    V=V+(vM(M(i),aM(i),x));
end
for i = 1:numel(w)
    V=V+(vw(w(i),aw(i),x));
end
for i = 1:numel(w_)
    V=V+(vw_(w_(i),aw_(i),x));
end
V= (V+C1*x+C2)/E*I;   

plot(x,V)
xlabel('Beam Length')
ylabel('Deflection')
title('Beam Deflection')

% Singularity function
function f=sing(x,a,n)
% This function evaluates the singularity function <x-a>^n.
if x<a
f= 0*x; %gives a list of zeros as long as x
else
f=(x-a).^n .* (x>=a); %(x>=a) is a list of ones and zeros.
end
end