** The following code is written in Form.
** The Form package can be downloaded from https://www.nikhef.nl/~form/
** It is very simple to install and run. The manual is available on-line at the same site.
**
** In this little template we want to verify Eqs.(29.22), (29.23), and (29.24) of Schwartz's book.
** If you install Form, you can run this code by simply typing (on a unix/linux machine):
** "form wz-wz.frm"
** The output will appear on the screen. If you want to save it in a file, just redirect the output, e.g.
** "form wz-wz.frm > wz-wz.out"
** I have also uploaded the output file for our inspection.
**
** We start by defining the building blocks of each amplitide, and then perform some simple algebra.
** Form knows about 4-vectors and their dot-product (denoted by a ".", e.g. p1.p2)
** If you use a different software, you will have to define 4-vectors and their dot-product.
** The print outs are in the corresponding output file "wz-wz.out".

nwrite statistics;

V p1,p2,p3,p4,ks,ku,p,k; VECTORS
I mu1,mu2,mu3,mu4,nu1,nu2,nu3,nu4,mu,nu,la,al,be; INDICES
S gwwz,gwwzz,s,t,u,mw,mz,[s-mw^2],[u-mw^2],[t-2*mw^2],[t-2*mz^2],[s+u],E,[t],[s],[u],n; 
** SYMBOLS (algebriac symbols)
V eps1,eps2,eps3,eps4;
CF Vwwz,Vwwzz,propw,As,Au,A4,A4test,g;
* COMMUTING FUNCTIONS: this means F*G=G*F if F and G are CF's.

** Define amplitude structure:
** [overall factor from couplings]*[polarization vectors]*[all the rest]

L Ms = -gwwz^2*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*As(p1,p2,p3,p4,mu1,mu2,mu3,mu4); 
** L for Local. This is like the input function. It's like y=f(x) <-> Local y = CFunction f(10)
** S CHANNEL
L Mu = gwwz^2*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*Au(p1,p2,p3,p4,mu1,mu2,mu3,mu4); 
** U CHANNEL
L M4 = gwwzz*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*A4(mu1,mu2,mu3,mu4); 4 POINT

** Express [all the rest] in terms of vertices and propagators:
** Vwwz: WWZ vertex
** Vwwzz: WWZZ vertex
** propw: W propagator

id As(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vwwz(-ks,p1,p2,nu,mu1,mu2)*Vwwz(-p3,ks,-p4,mu3,la,mu4)*propw(s,ks,nu,la); 
** id is IDENTIFY: make a substitution to that thing with the thing on the RHS
id Au(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vwwz(p1,-ku,-p4,mu1,nu,mu4)*Vwwz(-p3,ku,p2,mu3,la,mu2)*propw(u,ku,nu,la);
id A4(mu1?,mu2?,mu3?,mu4?)=Vwwzz(mu1,mu3,mu2,mu4);


** Feynman rules For vertices and propagators

id Vwwz(p1?,p2?,p3?,mu1?,mu2?,mu3?)=d_(mu1,mu2)*(p1(mu3)-p2(mu3))+d_(mu2,mu3)*(p2(mu1)-p3(mu1))+d_(mu3,mu1)*(p3(mu2)-p1(mu2));
id Vwwzz(mu?,nu?,al?,be?)=d_(al,mu)*d_(be,nu)+d_(al,nu)*d_(be,mu)-2*d_(al,be)*d_(mu,nu);
*** the ? is wildcard for pattern matching. if a term on the RHS has the adequate pattern that is in the wildcard in the LHS, then that corresponding part will be replaced in the RHS with the adequate substitution for n.
id propw(s,k?,nu?,la?)=1/[s-mw^2]*(-d_(nu,la)+k(nu)*k(la)/mw^2);
id propw(u,k?,nu?,la?)=1/[u-mw^2]*(-d_(nu,la)+k(nu)*k(la)/mw^2);


** Use that p*eps(p)=0, and substitute expression for polarization vectors

id p1.eps1=0;
id p2.eps2=0;
id p3.eps3=0;
id p4.eps4=0;

id eps1.p?=p1.p/mw+2*mw/[t-2*mw^2]*p3.p;
id eps2.p?=p2.p/mz+2*mz/[t-2*mz^2]*p4.p;
id eps3.p?=p3.p/mw+2*mw/[t-2*mw^2]*p1.p;
id eps4.p?=p4.p/mz+2*mz/[t-2*mz^2]*p2.p;


** Simplify result using simple kinematic relations and express it in terms of {s,t,u}

repeat; 
** START REPEAT
id p?.ks=p.p1+p.p2;
id p?.ku=p.p3-p.p2;
id p4.k?=-p3.k+p1.k+p2.k;

id p1.p1=mw^2;
id p2.p2=mz^2;
id p3.p3=mw^2;

id p1.p2=1/2*(s-mw^2-mz^2);
id p1.p3=-1/2*(t-2*mw^2);
id p2.p3=-1/2*(u-mw^2-mz^2);

id u=-s-t+2*mw^2+2*mz^2;

endrepeat;

bracket gwwz,gwwzz,[s-mw^2],[u-mw^2],[t-2*mw^2],[t-2*mz^2],mw,mz;
print;
.sort
*The .sort statement is a directive to FORM to execute a program block, sort the result (i.e., bring them in standard ordering), and prepare for further processing.
* .sort essentially means to continue

** Try further simplifications and expand denominators in the large-energy limit {s,t,u}>>{MW,MZ}

repeat;
id s*[s-mw^2]^-1=1+mw^2*[s-mw^2]^-1;
id u*[u-mw^2]^-1=1+mw^2*[u-mw^2]^-1;
id t*[t-2*mw^2]^-1=1+2*mw^2*[t-2*mw^2]^-1;
id t*[t-2*mz^2]^-1=1+2*mz^2*[t-2*mz^2]^-1;
endrepeat;
.sort

repeat;
id [s-mw^2]^-1=1/s*(1+mw^2/s+mw^4/s^2);
id [u-mw^2]^-1=1/u*(1+mw^2/u+mw^4/u^2);
id [t-2*mw^2]^-1=1/t*(1+2*mw^2/t+4*mw^4/t^2);
id [t-2*mz^2]^-1=1/t*(1+2*mz^2/t+4*mz^4/t^2);


** Keep only leading terms in the large-energy limit. To achieve that, rescale all invariants by a suitable energy factor (E^n) and then drop all terms proportional to negative powers of E. This requires to introduce different symbols for the rescaled invariants.

id t=E^2*[t];
id s=E^2*[s];
id u=E^2*[u];
id t^-1=E^-2*[t]^-1;
id s^-1=E^-2*[s]^-1;
id u^-1=E^-2*[u]^-1;
endrepeat;

id E^n?neg0_ = 0;

bracket gwwz,gwwzz,[s-mw^2],[u-mw^2],[t-2*mw^2],[t-2*mz^2],mw,mz,E;
print;
.sort


** Revert back to {s,t,u}. To compare with Schwartz's results, eliminate t in favor of {s,u}.

repeat;
id [s]=s;
id [s]^-1=s^-1;
id [t]=-s-u+(2*mw^2+2*mz^2)/E^2;
id [t]^-1=-[s+u]^-1*(1+2*(mw^2+mz^2)/E^2/[s+u]);
id [u]=u;
id [u]^-1=u^-1;
endrepeat;

id E^n?neg0_ = 0;

** Ms, Mu, And M4 in the final printout (see wz-wz.out) agree With Eq.(29.22), (29.23), And (29.24) of Schwratz's book (you just need to perform a little algebra to assemble terms with the same prefactor).

bracket gwwz,gwwzz,[s-mw^2],[u-mw^2],[t-2*mw^2],[t-2*mz^2],[s+u],mw,mz,E;
print;
.end