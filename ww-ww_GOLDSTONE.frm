nwrite statistics;

V p1,p2,p3,p4,ks,kt,p,k;
I mu1,mu2,mu3,mu4,nu1,nu2,nu3,nu4,mu,nu,la,al,be;
S Gggh,Ggggg,s,t,u,mw,mh,[s-mh^2],[t-mh^2],[t-2*mw^2],[s+u],E,[t],[s],[u],n;
V eps1,eps2,eps3,eps4;
CF Vggh,Vgggg,proph,Ash,Ath,A4,g;
** g for Goldstone
** Define amplitude structure:
** [overall factor from couplings]*[all the rest]

L Msh = -Gggh^2*Ash(p1,p2,p3,p4,mu1,mu2,mu3,mu4);
L Mth = -Gggh^2*Ath(p1,p2,p3,p4,mu1,mu2,mu3,mu4);
L M4 = Ggggg^2*A4(mu1,mu2,mu3,mu4);
L M = Msh+Mth+M4;

* G for coupling, g for Goldstone
** Express [all the rest] in terms of vertices and propagators:
** Vggh: Goldstone-Goldstone-Higgs vertex
** proph: H propagator

id Ash(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vggh(-ks,p1,p2,mu1,mu2)*Vggh(-p3,ks,-p4,mu3,mu4)*proph(s,ks);
id Ath(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vggh(p1,-kt,-p3,mu1,mu3)*Vggh(-p4,kt,p2,mu4,mu2)*proph(t,kt);
id A4(p1?,p2?,p3?,p4?,mu1?,mu2?,mu3?,mu4?)=Vgggg(mu1,mu3,mu4,mu2);

** Feynman rules For vertices and propagators

id Vggh(p1?,p2?,k?,mu1?,mu2?)=mh^2/(2*mw);
id Vgggg(p1?,p2?,k?,mu1?,mu2?)=mh^2*g^2/mw^2;

id proph(s,k?)=1/[s-mh^2];
id proph(t,k?)=1/[t-mh^2];





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
id s*[s-mh^2]^-1=1+mh^2*[s-mh^2]^-1;
id t*[t-mh^2]^-1=1+mh^2*[t-mh^2]^-1;
id t*[t-2*mw^2]^-1=1+2*mw^2*[t-2*mw^2]^-1;
endrepeat;
.sort

repeat;
id [s-mh^2]^-1=1/s*(1+mh^2/s);
id [t-mh^2]^-1=1/t*(1+mh^2/t);
id [t-2*mw^2]^-1=1/t*(1+2*mw^2/t);


** Keep only leading terms in the large-energy limit. To achieve that, rescale all invariants by a suitable energy factor (E^n) and then drop all terms proportional to negative powers of E. This requires to introduce different symbols for the rescaled invariants.

id t=E^2*[t];
id s=E^2*[s];
id u=E^2*[u];
id t^-1=E^-2*[t]^-1;
id s^-1=E^-2*[s]^-1;
id u^-1=E^-2*[u]^-1;
endrepeat;

id E^n?neg0_ = 0;


** Revert back to {s,t,u}. To compare with Schwartz's results, eliminate t in favor of {s,u}.

repeat;
id [s]=s;
id [s]^-1=s^-1;
id [t]=-s-u+(4*mw^2)/E^2;
id [t]^-1=-[s+u]^-1*(1+4*mw^2/E^2/[s+u]);
id [u]=u;
id [u]^-1=u^-1;
endrepeat;

id E^n?neg0_ = 0;
id E = 1;
** id gwwh^2 = g^2*mw^2;

bracket Gggh,Ggggg,[s-mh^2],[t-mh^2],[t-2*mw^2],mw,mh,E;
print;
.end