nwrite statistics;

V p1,p2,p3,p4,ks,kt,p,k;
I mu1,mu2,mu3,mu4,nu1,nu2,nu3,nu4,mu,nu,la,al,be;
S gwwz,gwwg,gwwww,g,s,t,u,mw,mz,[s-mz^2],[t-mz^2],[u-2*mw^2],[s+t],E,[t],[s],[u],n;
V eps1,eps2,eps3,eps4;
CF Vwwz,Vwwg,Vwwww,propz,propg,Asg,Asz,Atg,Atz,A4,A4test;

** Define amplitude structure:
** [overall factor from couplings]*[polarization vectors]*[all the rest]

L Msz = -gwwz^2*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*Asz(p1,p2,p3,p4,mu1,mu2,mu3,mu4);
L Msg = -gwwg^2*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*Asg(p1,p2,p3,p4,mu1,mu2,mu3,mu4);
L Ms = Msz+Msg;

L Mtz = -gwwz^2*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*Atz(p1,p2,p3,p4,mu1,mu2,mu3,mu4);
L Mtg = -gwwg^2*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*Atg(p1,p2,p3,p4,mu1,mu2,mu3,mu4);
L Mt = Mtz+Mtg;

L M4 = gwwww*eps1(mu1)*eps2(mu2)*eps3(mu3)*eps4(mu4)*A4(mu1,mu2,mu3,mu4);

L M = Ms+Mt+M4;

** Express [all the rest] in terms of vertices and propagators:
** Vwwz: WWZ vertex
** Vwwg: WWg vertex
** Vwwww: WWWW vertex
** propz: Z propagator
** propg: gamma propagator

id Asz(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vwwz(p1,p2,-ks,mu1,mu2,nu)*Vwwz(-p4,-p3,ks,mu4,mu3,la)*propz(s,ks,nu,la);
id Asg(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vwwg(p1,p2,-ks,mu1,mu2,nu)*Vwwg(-p4,-p3,ks,mu4,mu3,la)*propg(s,ks,nu,la);

id Atz(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vwwz(p1,-p3,-kt,mu1,mu3,nu)*Vwwz(-p4,p2,kt,mu4,mu2,la)*propz(t,kt,nu,la);
id Atg(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vwwg(p1,-p3,-kt,mu1,mu3,nu)*Vwwg(-p4,p2,kt,mu4,mu2,la)*propg(t,kt,nu,la);

id A4(mu1?,mu2?,mu3?,mu4?)=Vwwww(mu1,mu3,mu4,mu2);


** Feynman rules For vertices and propagators

id Vwwz(p1?,p2?,p3?,mu1?,mu2?,mu3?)=d_(mu1,mu2)*(p1(mu3)-p2(mu3))+d_(mu2,mu3)*(p2(mu1)-p3(mu1))+d_(mu3,mu1)*(p3(mu2)-p1(mu2));
id Vwwg(p1?,p2?,p3?,mu1?,mu2?,mu3?)=d_(mu1,mu2)*(p1(mu3)-p2(mu3))+d_(mu2,mu3)*(p2(mu1)-p3(mu1))+d_(mu3,mu1)*(p3(mu2)-p1(mu2));
id Vwwww(mu?,nu?,al?,be?)=-d_(be,al)*d_(mu,nu)-d_(nu,al)*d_(be,mu)+2*d_(mu,al)*d_(be,nu);

id propz(s,k?,nu?,la?)=[s-mz^2]^-1*(-d_(nu,la)+k(nu)*k(la)/mz^2);
id propz(t,k?,nu?,la?)=[t-mz^2]^-1*(-d_(nu,la)+k(nu)*k(la)/mz^2);
id propg(s,k?,nu?,la?)=s^-1*(-d_(nu,la));
id propg(t,k?,nu?,la?)=t^-1*(-d_(nu,la));


** Use that p*eps(p)=0, and substitute expression for polarization vectors

id p1.eps1=0;
id p2.eps2=0;
id p3.eps3=0;
id p4.eps4=0;

id eps1.p?=p1.p/mw+2*mw/[u-2*mw^2]*p4.p;
id eps2.p?=p2.p/mw+2*mw/[u-2*mw^2]*p3.p;
id eps3.p?=p3.p/mw+2*mw/[u-2*mw^2]*p2.p;
id eps4.p?=p4.p/mw+2*mw/[u-2*mw^2]*p1.p;


** Simplify result using simple kinematic relations and express it in terms of {s,t,u}

repeat;
id p?.ks=p.p1+p.p2;
id p?.kt=p.p1-p.p3;
id p4.k?=-p3.k+p1.k+p2.k;

id p1.p1=mw^2;
id p2.p2=mw^2;
id p3.p3=mw^2;

id p1.p2=1/2*(s-2*mw^2);
id p1.p3=-1/2*(t-2*mw^2);
id p2.p3=-1/2*(u-2*mw^2);

id u=-s-t+4*mw^2;

endrepeat;


** Try further simplifications and expand denominators in the large-energy limit {s,t,u}>>{MW,MZ}

repeat;
id s*[s-mz^2]^-1=1+mz^2*[s-mz^2]^-1;
id t*[t-mz^2]^-1=1+mz^2*[t-mz^2]^-1;
id u*[u-2*mw^2]^-1=1+2*mw^2*[u-2*mw^2]^-1;
endrepeat;
.sort

repeat;
id [s-mz^2]^-1=1/s*(1+mz^2/s+mz^4*s^-2);
id [t-mz^2]^-1=1/t*(1+mz^2/t+mz^4*t^-2);
id [u-2*mw^2]^-1=1/u*(1+2*mw^2/u+4*mw^4/u^2);


** Keep only leading terms in the large-energy limit. To achieve that, rescale all invariants by a suitable energy factor (E^n) and then drop all terms proportional to negative powers of E. This requires to introduce different symbols for the rescaled invariants.

id t=E^2*[t];
id s=E^2*[s];
id u=E^2*[u];
id t^-1=E^-2*[t]^-1;
id s^-1=E^-2*[s]^-1;
id u^-1=E^-2*[u]^-1;
endrepeat;

id E^n?neg0_ = 0;


** Revert back to {s,t,u}. To compare with Schwartz's results, eliminate u in favor of {s,t}.

repeat;
id [s]=s;
id [s]^-1=s^-1;
id [t]=t;
id [t]^-1=t^-1;
id [u]=-s-t+(4*mw^2)/E^2;
id [u]^-1=-[s+t]^-1*(1+4*mw^2/E^2/[s+t]);
endrepeat;

id E^n?neg0_ = 0;
id E = 1;
** id gwwg^2 = -gwwz^2 + g^2;
** id gwwz^2 = g^2*mw^2*mz^-2;
** id gwwww = g^2;

bracket gwwz,gwwg,gwwww,[s-mz^2],[t-mz^2],[u-2*mw^2],mw,mz,E;
print;
.end