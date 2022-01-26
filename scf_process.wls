#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["/home/physicsmaster/Documents/FXW/STM_votex"];


maxmu=159/2;R=30.;temp=0.01;rstep=0.01;kf=1.;xi0=7.6;ef=12.;omgd=ecut=10.;
Nstates=60;
ris=Range[0.,R,rstep];(*r coordinates*)iprint=0(*1输出每步流程，0不输出，对于多核并行不建议输出*);
cycle=1;
rangecut=8.;


alpha[i_,m_]:=N@BesselJZero[m,i];(*Bessel函数第i个零点*)
scN[i_,m_]:=(R BesselJ[m+1,alpha[i,m]])/Sqrt[2];(*归一化系数*)
basisrange[kfxi_,ef_,maxmu_:319/2,Nstates_:100]:=Flatten[#,1]&/@Transpose[{{#-(Nstates/2-1),#+Nstates/2}&/@(DeleteDuplicatesBy[Select[Flatten[ParallelTable[{mu,i,ef/(kfxi)^2 (alpha[i,Abs[mu-1/2]]/R)^2-ef},{mu,1/2,maxmu},{i,1,1000,1}],1],Abs[#[[3]]]<0.1&],First][[All,2]]),Range[1/2,maxmu,1]}];
jrange[ecut_,ef_,kfxi_,mu_,sigma_]:=Block[{low,high,exceed,count=0},{low=Floor@#[[1]],high=Ceiling@#[[2]]}&@{1+kfxi R Sqrt[1-ecut/ef]/Pi+1/4-(mu+sigma 1/2)/2+(mu+sigma 1/2)^2/(2 Pi kfxi R Sqrt[1-ecut/ef]),1+kfxi R Sqrt[1+ecut/ef]/Pi+1/4-(mu+sigma 1/2)/2+(mu+sigma 1/2)^2/(2 Pi kfxi R Sqrt[1+ecut/ef])};
If[high-low<Nstates-1,Print["ecut is too small!"],If[high-low>Nstates,While[high-low>=(Nstates-1)\[And]count<10,exceed=high-low-Nstates;low=(low+exceed/2);
high=(high-exceed/2);count++;];{IntegerPart@low,IntegerPart@high,mu},{IntegerPart@low,IntegerPart@high,mu}]]];
phi[i_,m_,r_]:=1/scN[i,m] BesselJ[m,alpha[i,m]/R r];
phid[i_,m_,r_]:=1/scN[i,m] D[BesselJ[m,alpha[i,m]/R r],r];
(*phini[mu_,switch_]:=If[switch\[Equal]-1,ParallelTable[phi[i,Abs[mu-1/2],ris],{i,Nstates}],ParallelTable[phi[i,Abs[mu+1/2],ris],{i,Nstates}]];*)
phini[mu_,switch_,basisrange_]:=If[switch==-1,ParallelTable[Quiet[phi[i,Abs[mu-1/2],ris]],{i,basisrange[[1]],basisrange[[2]]}],ParallelTable[Quiet[phi[i,Abs[mu+1/2],ris]],{i,basisrange[[1]],basisrange[[2]],brstep}]]
(*do not subtract Ef*)
ffermi[iotaE_,Ttempe_]:=If[Ttempe==0,If[iotaE>0,0,1],1/(If[(1.75(iotaE))/Ttempe>=-200.,If[(1.75(iotaE))/Ttempe<=200.,Exp[(1.75(iotaE))/Ttempe],Infinity],0.]+1)](*temp=t/tc;*);
fprecswitch[prec_]:=N[1.,prec];
processeigen[HBdG_]:=Module[{data,dataeigbef},If[iprint==1,Print["Calculating eigenvalues of HBdG",Now],Null];data=Eigensystem[HBdG];dataeigbef=Transpose[data];dataeigbef={#[[1]],Normalize[#[[2]]]}&/@Sort[dataeigbef,Function[{x,y},If[x[[1]]>0\[And]y[[1]]<0,True,Abs[x[[1]]]<Abs[y[[1]]]]]];dataeigbef];

HBdG[localmu_]:=Module[{deltamat,mu,hbdginitmat},mu=localmu;If[iprint==1,Print["Calculating deltaij",Now],Null];
deltamat=(*Table*)phiniv[mu,-1] . If[cycle==1,SparseArray[Table[{i,i}->rstep ris[[i]] Abs[deltainit[[i]]],{i,Length@ris}]](*DiagonalMatrix[rstep*ris*Abs[deltainit]]*),SparseArray[Table[{i,i}->rstep ris[[i]]Abs[deltanow[[i]]],{i,Length@ris}]](*DiagonalMatrix[rstep*ris*Abs[deltanow]]*)] . Transpose[phiniv[mu,1]];If[iprint==1,Print["Generating Full HBdG",Now],Null];hbdginitmat=hbdginitmatall[[mu+1/2]];
hbdginitmat[[1;;Nstates,Nstates+1;;2Nstates]]=deltamat;hbdginitmat[[Nstates+1;;2Nstates,1;;Nstates]]=Transpose@deltamat;
hbdginitmat];

HBdGinit[br_]:=Module[{hbdginit,Tminusmat,Tplusmat,mu},mu=br[[3]];
(*with experiment parameters*)hbdginit=DiagonalMatrix[Flatten[{Table[(ef/(kf xi0)^2 (alpha[i,Abs[mu-1/2]]/R)^2-ef),{i,br[[1]],br[[2]],brstep}],Table[(-ef/(kf xi0)^2 (alpha[i,Abs[mu+1/2]]/R)^2+ef),{i,br[[1]],br[[2]],brstep}]},1]]];
Eetabef[dataeigbef_,eta_]:=dataeigbef[[eta,1]];
fdeltar[mu_,dataeigbef_]:=Module[{c1,fE,fmE,tab1,c2,t,tab2,tab3,cutnumberp,cutnumbern},If[iprint==1,Print["Gernerating new delta(r)",Now],Null];(*deltax*)Neta=2Nstates;cutnumberp=Length@Select[dataeigbef[[All,1]],#>0\[And]#<omgd&];
cutnumbern=Length@Select[dataeigbef[[All,1]],#<0\[And]#>-omgd&];If[cutnumberp==0\[Or]cutnumbern==0,0.,c1=Join[dataeigbef[[1;;cutnumberp,2,1;;Nstates]],dataeigbef[[Nstates+1;;Nstates+cutnumbern,2,1;;Nstates]]];
c2=Join[dataeigbef[[1;;cutnumberp,2,Nstates+1;;2Nstates]],dataeigbef[[Nstates+1;;Nstates+cutnumbern,2,Nstates+1;;2Nstates]]];
fE=SparseArray[Table[{eta,eta}->If[eta<=cutnumberp,ffermi[Eetabef[dataeigbef,eta],temp],ffermi[Eetabef[dataeigbef,eta-cutnumberp+Nstates],temp]],{eta,1,cutnumberp+cutnumbern,1}]];fmE=SparseArray[Table[{eta,eta}->If[eta<=cutnumberp,ffermi[-Eetabef[dataeigbef,eta],temp],ffermi[-Eetabef[dataeigbef,eta-cutnumberp+Nstates],temp]],{eta,1,cutnumberp+cutnumbern,1}]];
tab1=Transpose[c1] . fE . Conjugate[c2];
tab2=Transpose[c1] . fmE . Conjugate[c2];
tab3=Transpose[c2] . Conjugate[c2];(*definition changed*)(g*Diagonal[Transpose[phiniv[mu,-1][[All,1;;(rangecut/rstep+1)]]] . (tab2-tab1) . Conjugate[phiniv[mu,1][[All,1;;(rangecut/rstep+1)]]]])]];

fcycle[]:=Module[{t1},If[AllTrue[Flatten[{g,deltainit,deltadata,cycle}],NumberQ],While[If[cycle==1,True,(Norm[deltadata[[-2,1;;(R/rstep)]]-deltadata[[-1,1;;(R/rstep)]]]/Norm[deltadata[[-2,1;;(R/rstep)]]]>0.0001)\[And]cycle<50],t1=Now(*Max[deltadata[[-1]]]*);dataeigbefall=(*ParallelTable*)ParallelTable[processeigen[HBdG[mu]],{mu,1/2,maxmu,1}];deltanow=Abs@Flatten@Append[(*ParallelSum*)ParallelSum[fdeltar[mu,dataeigbefall[[(2 mu+1)/2]]],{mu,1/2,maxmu,1}],Table[1.,Length[Range[0,R,rstep]]-Length[Range[0,rangecut,rstep]]]];
AppendTo[deltadata,deltanow];
(*If[cycle>1,Print[{cycle,dataeigbefall[[1]][[All,1]]//Max(*the max eigenvalue*),Max@deltanow(*max energy gap*),Norm[deltadata[[-2]]-deltadata[[-1]]]/Norm[deltadata[[-1]]],ListLinePlot[Transpose[{Select[ris,#\[LessEqual]1&],#[[1;;Length@Select[ris,#\[LessEqual]1&]]]}]&/@deltadata,PlotRange\[Rule]All]}],Null]*);(*Print[Now];*)cycle++;(*str=OpenWrite["g=0.3338(2Pi),maxmu=319,Nstates=100,R=10,rstep=0.2,ef=8,temp=0.01,deltanow,dataeigenbefall"];
Write[str,{deltanow,dataeigbefall}];
Close[str];*)tdiff=Now-t1;Print[{Length[deltadata],ef,g,temp,If[cycle==1,0.,Norm[deltadata[[-2]]-deltadata[[-1]]]/Norm[deltadata[[-1]]]],tdiff,ListLinePlot[Transpose[{ris,deltadata[[-1]]}],PlotRange->Full],ListLinePlot[Transpose[{ris,#}]&/@deltadata,PlotRange->{{0,1.2},{0,1.5}}]}]];dataeigbefall=(*ParallelTable*)ParallelTable[processeigen[HBdG[mu]],{mu,1/2,maxmu,1}](*when delta converges,calculate current dataeigbefall*);,Print["Check your constants!"]]];


(*calculate Nstates rangecut*)
Nstatesrange=Nstates;
sortvector[vec_]:=Flatten[{Sort[Select[vec,#1>0&]],Sort[Select[vec,#1<0&],Function[{x,y},Abs[x]<Abs[y]]]}];
deltaana[xi0_:10,xi1_:2,r_]:=If[r==0,0.`,1/(1+(xi0 Exp[-(r/xi1)])/r)];
deltainit=Table[deltaana[10,2,r],{r,0,R,rstep}];
br=Table[jrange[ecut,ef,kf*xi0,mu,-1],{mu,1/2,maxmu}];
brstep=1;


(*read basis functions*)
(* str=OpenRead["expcutbasis,kfxi=7.6,ef=12.0,ecut=10,Nstatesrange=60,Nstates=60,R=30.0,rstep=0.01,maxmu=79.5,phinivm"(*"Nstates=100,R=10,rstep=0.01,maxmu=319,phinivm"*)]; *)
str = OpenRead["expcutbasis,kfxi=7.6,ef=12.0,ecut=10,Nstatesrange=60,Nstates=60,R=30.0,rstep=0.01,maxμ=79.5,ϕnivm"];
phitemp=Read[str];
Table[phiniv[i,-1]=(phitemp[[(2 i+1)/2]][[All,1;;(R/rstep+1)]]),{i,1/2,maxmu,1}];
Table[phiniv[i,1]=(phitemp[[(2 i+1)/2+1]][[All,1;;(R/rstep+1)]]),{i,1/2,maxmu-1,1}];
phiniv[maxmu,1]=Chop[phini[maxmu,1,br[[-1]](*{1,100,319/2}*)],10.^-50];
(*add the maxmu term to +1 beseel function*)
(*Close[str];*)
(*ClearAll[phitemp,str];*)


hbdginitmatall = ParallelTable[HBdGinit[br[[i]]], {i, 1, Length[br]}];(*产生多个\[Mu]的初始化矩阵*)


dataeigbefall = Table[processeigen[HBdG[mu]], {mu, 1/2, maxmu, 1}];



