#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["E:\\MMA_data"];


max\[Mu] = 
 319/2; \[Omega]d = 10.; R = 10.(*can be 10, decided by \[Phi]nvi \
generated before*); kf = 8; xi0 = 2; ef = 1/2 kf xi0;
g = Which[kf xi0 == 2, 0.221*10/1.35, kf xi0 == 4, 0.221*4.35/1.35, 
   kf xi0 == 6, 0.32*1.365, kf xi0 == 8, 0.31, kf xi0 == 10, 0.238, 
   kf xi0 == 12, 0.21*0.91, kf xi0 == 14, 0.159, kf xi0 == 16, 
   0.275/2];
Nstates = 100; rstep = 0.01;(*calculation parameters*)
ris = Range[0., R, rstep];(*r coordinates*)
iprint = 0(*1输出每步流程，0不输出，对于多核并行不建议输出*);


\[Alpha][i_,m_]:=N@BesselJZero[m,i];(*Bessel函数第i个零点*)
\[ScriptCapitalN][i_,m_]:=(R BesselJ[m+1,\[Alpha][i,m]])/Sqrt[2];(*归一化系数*)
\[Phi][i_,m_,r_]:=1/\[ScriptCapitalN][i,m] BesselJ[m,\[Alpha][i,m]/R r];
\[Phi]d[i_,m_,r_]:=1/\[ScriptCapitalN][i,m] D[BesselJ[m,\[Alpha][i,m]/R r],r];
\[Phi]ni[\[Mu]_,switch_]:=If[switch==-1,ParallelTable[Quiet[\[Phi][i,Abs[\[Mu]-1/2],ris]],{i,Nstates}],ParallelTable[Quiet[\[Phi][i,Abs[\[Mu]+1/2],ris]],{i,Nstates}]];
(*do not subtract Ef*)
ffermi[\[Iota]E_,Ttempe_]:=If[Ttempe==0,If[\[Iota]E>0,0,1],1/(If[(1.75(\[Iota]E))/Ttempe>=-200.,If[(1.75(\[Iota]E))/Ttempe<=200.,Exp[(1.75(\[Iota]E))/Ttempe],Infinity],0.]+1)](*temp =t/tc;*)
fprecswitch[prec_]:=N[1.,prec];
processeigen[HBdG_]:=Module[{data,dataeigbef},If[iprint==1,Print["Calculating eigenvalues of HBdG",Now],Null];data=Eigensystem[HBdG];dataeigbef=Transpose[data];dataeigbef={#[[1]],Normalize[#[[2]]]}&/@Sort[dataeigbef,Function[{x,y},If[x[[1]]>0\[And]y[[1]]<0,True,Abs[x[[1]]]<Abs[y[[1]]]]]];dataeigbef];
HBdG[local\[Mu]_]:=Module[{\[CapitalDelta]mat,\[Mu],hbdginitmat},\[Mu]=local\[Mu];If[iprint==1,Print["Calculating \[CapitalDelta]ij",Now],Null];
\[CapitalDelta]mat=(*Table*)\[Phi]niv[\[Mu],-1] . If[cycle==1,SparseArray[Table[{i,i}->rstep ris[[i]] Abs[\[CapitalDelta]init[[i]]],{i,Length@ris}]],SparseArray[Table[{i,i}->rstep ris[[i]]Abs[\[CapitalDelta]now[[i]]],{i,Length@ris}]]] . Transpose[\[Phi]niv[\[Mu],1]];If[iprint==1,Print["Generating Full HBdG",Now],Null];hbdginitmat=hbdginitmatall[[Abs[\[Mu]+If[\[Mu]>0,1/2,-1/2]]]];
hbdginitmat[[1;;Nstates,Nstates+1;;2Nstates]]=\[CapitalDelta]mat;hbdginitmat[[Nstates+1;;2Nstates,1;;Nstates]]=Transpose@\[CapitalDelta]mat;
hbdginitmat];
HBdGinit[local\[Mu]_]:=Module[{hbdginit,Tminusmat,Tplusmat,\[Mu]},\[Mu]=local\[Mu];
hbdginit=Table[0.,{i,2Nstates},{j,2Nstates}];hbdginit=SparseArray[Table[{i,i}->If[i<Nstates+1,(1/(2kf xi0) (\[Alpha][i,Abs[\[Mu]-1/2]]/R)^2-ef),(-1/(2kf xi0) (\[Alpha][i-Nstates,Abs[\[Mu]+1/2]]/R)^2+ef)],{i,1,2Nstates}]]]
(*str=OpenRead["Nstates=100,R=10,rstep=0.01,max\[Mu]=319,\[Phi]dnivm"];\[Phi]temp=Read[str];
Table[\[Phi]dniv[i,-1]=\[Phi]temp[[(2i+1)/2]],{i,1/2,319/2,1}];
Close[str];
str=OpenRead["Nstates=100,R=10,rstep=0.01,max\[Mu]=319,\[Phi]dnivp"];\[Phi]temp=Read[str];
Table[\[Phi]dniv[i,1]=\[Phi]temp[[(2i+1)/2]],{i,1/2,319/2,1}];
Close[str];*)
E\[Eta]bef[dataeigbef_,\[Eta]_]:=dataeigbef[[\[Eta],1]];
f\[CapitalDelta]r[\[Mu]_,dataeigbef_]:=Module[{c1,fE,fmE,tab1,c2,t,tab2,tab3,cutnumberp,cutnumbern},If[iprint==1,Print["Gernerating new \[CapitalDelta](r)",Now],Null];(*\[CapitalDelta]x*)N\[Eta]=2Nstates;cutnumberp=Length@Select[dataeigbef[[All,1]],#>0\[And]#<\[Omega]d&];

cutnumbern=Length@Select[dataeigbef[[All,1]],#<0\[And]#>-\[Omega]d&];If[cutnumberp==0\[Or]cutnumbern==0,0.,c1=Join[dataeigbef[[1;;cutnumberp,2,1;;Nstates]],dataeigbef[[Nstates+1;;Nstates+cutnumbern,2,1;;Nstates]]];
c2=Join[dataeigbef[[1;;cutnumberp,2,Nstates+1;;2Nstates]],dataeigbef[[Nstates+1;;Nstates+cutnumbern,2,Nstates+1;;2Nstates]]];
fE=SparseArray[Table[{\[Eta],\[Eta]}->If[\[Eta]<=cutnumberp,ffermi[E\[Eta]bef[dataeigbef,\[Eta]],temp],ffermi[E\[Eta]bef[dataeigbef,\[Eta]-cutnumberp+Nstates],temp]],{\[Eta],1,cutnumberp+cutnumbern,1}]];fmE=SparseArray[Table[{\[Eta],\[Eta]}->If[\[Eta]<=cutnumberp,ffermi[-E\[Eta]bef[dataeigbef,\[Eta]],temp],ffermi[-E\[Eta]bef[dataeigbef,\[Eta]-cutnumberp+Nstates],temp]],{\[Eta],1,cutnumberp+cutnumbern,1}]];
tab1=Transpose[c1] . fE . Conjugate[c2];
tab2=Transpose[c1] . fmE . Conjugate[c2];
tab3=Transpose[c2] . Conjugate[c2];(*definition changed*)
(g Diagonal[Transpose[\[Phi]niv[\[Mu],-1]] . (tab2-tab1) . Conjugate[\[Phi]niv[\[Mu],1]]])/(2Pi)]]

fcycle[]:=Module[{},While[If[cycle==1,True,(Norm[\[CapitalDelta]data[[-2]]-\[CapitalDelta]data[[-1]]]/Norm[\[CapitalDelta]data[[-1]]]>0.01)\[And]cycle<50],t1=Now(*Max[\[CapitalDelta]data[[-1]]]*);dataeigbefall=ParallelTable[processeigen[HBdG[\[Mu]]],{\[Mu],-1/2,-max\[Mu],-1}];\[CapitalDelta]now=ParallelSum[f\[CapitalDelta]r[\[Mu],dataeigbefall[[Abs[\[Mu]+If[\[Mu]>0,1/2,-1/2]]]]],{\[Mu],-1/2,-max\[Mu],-1}];
AppendTo[\[CapitalDelta]data,\[CapitalDelta]now];
(*If[cycle>1,Print[{cycle,dataeigbefall[[1]][[All,1]]//Max(*the max eigenvalue*),Max@\[CapitalDelta]now(*max energy gap*),Norm[\[CapitalDelta]data[[-2]]-\[CapitalDelta]data[[-1]]]/Norm[\[CapitalDelta]data[[-1]]],ListLinePlot[Transpose[{Select[ris,#\[LessEqual]1&],#[[1;;Length@Select[ris,#\[LessEqual]1&]]]}]&/@\[CapitalDelta]data,PlotRange\[Rule]All]}],Null]*);(*Print[Now];*)cycle++;(*str=OpenWrite["g=0.3338(2Pi),maxmu=319,Nstates=100,R=10,rstep=0.2,ef=8,temp=0.01,\[CapitalDelta]now,dataeigenbefall"];
Write[str,{\[CapitalDelta]now,dataeigbefall}];
Close[str];*)tdiff=Now-t1;];dataeigbefall=ParallelTable[processeigen[HBdG[\[Mu]]],{\[Mu],-1/2,-max\[Mu],-1}];(*when \[CapitalDelta] converges, calculate current dataeigbefall*)];


str = OpenRead[
  "Nstates=100,R=10,rstep=0.01,max\[Mu]=319,\[Phi]nivm"]; \[Phi]temp =
  Read[str];
Table[\[Phi]niv[i, -1] = (\[Phi]temp[[(2 i + 1)/2]](*[[All,1;;
    501]]*)), {i, 1/2, 319/2, 1}];
Close[str];
Table[\[Phi]niv[\[Mu], 1] = \[Phi]niv[\[Mu] + 1, -1], {\[Mu], 1/2, 
   max\[Mu] - 1, 1}];
\[Phi]niv[max\[Mu], 1] = Chop[\[Phi]ni[max\[Mu], 1], 1.*^-50];
Table[\[Phi]niv[i,-1]=\[Phi]niv[-i,1];\[Phi]niv[i,1]=\[Phi]niv[-i,-1],{i,-1/2,-max\[Mu],-1}];


(* nsc calculation *)
cycle = 1;
\[CapitalDelta]init = Tanh[3 ris];(*initial gap*)
\[CapitalDelta]data = {\[CapitalDelta]init};
(*改变kf xi0 以及ef要重新初始化矩阵*)
hbdginitmatall = 
  ParallelTable[
   HBdGinit[i], {i, -1/2, - max\[Mu], -1}](*产生多个\[Mu]的初始化矩阵*);
negeigbefall=Table[processeigen[HBdG[i]], {i, -1/2, -max\[Mu], -1}]


(* sc calculation *)
Table[temp = 0.01; cycle = 1;
 \[CapitalDelta]init = Tanh[3 ris];(*initial gap*)
 \[CapitalDelta]data = {\[CapitalDelta]init}; kf = 1;
 xi0 = xi; ef = kf xi0 /2; 
 g = Which[kf xi0 == 2, 0.221*10/1.35, kf xi0 == 4, 0.221*4.35/1.35, 
   kf xi0 == 6, 0.32*1.365, kf xi0 == 8, 0.31, kf xi0 == 10, 0.238, 
   kf xi0 == 12, 0.21*0.91, kf xi0 == 14, 0.159, kf xi0 == 16, 
   0.275/2]; 
 hbdginitmatall = 
  ParallelTable[HBdGinit[i], {i, -1/2, -max\[Mu], -1}]; fcycle[]; 
 Table[dataeigbefall[[i]][[1, 1]], {i, 7}], {xi, {16}}]
