(* ::Package:: *)

(* ::Input:: *)
(*(*Greece Pops*)t0=7;(*Life cyrcle of COVID19*)N0=10816286;*)
(*P=((86440/tmax)+79.554)/N0;*)
(*m=120297/(N0*tmax);*)
(*tsic=30;*)
(*atr[t_]:=(0.215) ((1+(t/1825))/10);*)
(*sig[t_]:=1/t0*(1+t/120);*)
(*b1=2.415;*)
(*b2=2.415/5;*)
(*k2=1-atr[t];*)
(*k1=1/3;*)
(*tquar=60;*)
(*rec=1/t0;*)
(*qu=9.5/10;*)
(*ttq[t_]:=(1/tquar)+((Exp[(t-210)/60]))/(10^9);*)
(*dis[t_]:=1000 (((0.93)*(1+Cos[(t/900)*2*Pi])/2)*(1/(Exp[(t-275)/30]+1))+0.07);*)
(*Plot[dis[t],{t,1,100}]*)
(*T=0.30;*)
(*(*Tourist influx parameter*)*)
(*tt[t_]:=Abs[Sin[Pi*(t/365)]]*(1/(Exp[(t-450)/15]+1));*)
(**)
(*(*Seasonality*)*)
(*tmax=365;*)
(*mt[t_]:=Exp[7t/(8tmax)];(*Mutation*)*)
(*tp[t_]:=mt[t]*(1+(0.85)Sin[6*Pi*((t-106)/(365))]);*)
(*Plot[tp[t],{t,1,800}]*)
(*i0=1;*)
(*e0=0;*)
(*s0=0;*)
(*s0=N0-i0-e0;*)
(**)


(* ::Input:: *)
(*(*vaccination*)*)
(*v[t_]:=1-(3(1/(Exp[(-t+395)/20]+1)+1/(Exp[(-t+460)/100]+1))/8);*)
(*sig2[t_]:=sig[t]*v[t];*)
(*atr2[t_]:=atr[t]*v[t];*)


(* ::Input:: *)
(*sol=NDSolve[{s'[t]==P+1000*tt[t]-m*s[t]-T*(b1*tp[t]*e[t-2t0]+b2*tp[t]*i[t-2t0])*s[t]/(s[t]+i[t]+e[t]+r[t]+d[t])-dis[t]*s[t]*i[t]/(s[t]+i[t]+e[t]+r[t]+d[t])+rec*r[t]+ttq[t]*d[t],e'[t]==50*tt[t]-(m+k1+sig2[t])*e[t]+T*(b1*tp[t]*e[t-2t0]+b2*tp[t]*i[t-2t0])*s[t]/(s[t]+i[t]+e[t]+r[t]+d[t])-dis[t]*e[t]*i[t]/(s[t]+i[t]+e[t]+r[t]+d[t]),*)
(*i'[t]==sig2[t]*e[t]-m*i[t]-atr2[t]*i[t]-k2*i[t]-qu*i[t],*)
(*r'[t]==k1*e[t]+k2*i[t]-m*r[t]-rec*r[t],d'[t]==dis[t]*(e[t]+s[t])*i[t]/(s[t]+i[t]+e[t]+r[t]+d[t])+qu*i[t]-ttq[t]*d[t],s[1]==s0,e[1]==e0,i[1]==i0,r[1]==0,d[1]==0},{s[t],e[t],i[t],r[t],d[t]},{t,1,1000},MaxSteps->Infinity];*)
(*ssol[t_]:=s[t]/.sol[[1,1]];*)
(*esol[t_]:=e[t]/.sol[[1,2]];*)
(*isol[t_]:=i[t]/.sol[[1,3]];*)
(*rsol[t_]:=r[t]/.sol[[1,4]];*)
(*dsol[t_]:=d[t]/.sol[[1,5]];*)


(* ::Input:: *)
(*(*spreading rate without prevention*)*)
(*R0=((b1+b2)/(2 m+sig[t]+1+k1));*)
(*Print["R0 ="]*)
(*R0//N*)
(*(*spreading rate with prevention*)*)
(**)
(*RT=(T (b1+b2)/(2 m+sig[t]+1+k1));*)
(*Print["RT ="]*)
(*RT//N*)
(**)
(*Print["Number of Deaths (Yearly):"]*)
(*Integrate[atr2[t]*isol[t],{t,1,800}]//N*)
(*sol;*)
(*Plot[ssol[t],{t,1,500},PlotStyle->RGBColor[0.7,0.73,0.19]]*)
(*Ist=Plot[isol[t],{t,1,800},PlotStyle->RGBColor[0.2,0.7,0.1],PlotRange->All,PlotLabel->"Infected from Day 10",AxesLabel->{"Days after Day 10","Number of infected people per day"}]*)
(*Plot[dsol[t],{t,1,500},PlotStyle->RGBColor[0.24,0.73,0.8]]*)
(*Plot[esol[t],{t,1,500},PlotStyle->RGBColor[0.64,0.33,0.59],PlotRange->All]*)
(*Plot[rsol[t],{t,1,500},PlotStyle->RGBColor[0.4,0.3,0.9]]*)


(* ::Input:: *)
(*vacset=Import["GreeceStat.xlsx","Data"];*)
(*vacset=Flatten[vacset];*)
(*vacset=Drop[vacset,1];*)
(*Est=ListPlot[vacset]*)
(*Show[Ist,Est]*)


Nn[t_]=(ssol[t]+isol[t]+rsol[t]+esol[t]+dsol[t])^(-1);
sn[t_]=ssol[t]*Nn[t];
rn[t_]=rsol[t]*Nn[t];
in[t_]=isol[t]*Nn[t];
en[t_]=esol[t]*Nn[t];
dn[t_]=dsol[t]*Nn[t];



Sen[t_]:=-(sn[t]*Log[sn[t]]+en[t]*Log[en[t]]+rn[t]*Log[rn[t]]+in[t]*Log[in[t]]+dn[t]*Log[dn[t]])
DSen[t_]:=Sen[t]-Sen[t-1]
Plot[{DSen[t],Sen[t]},{t,2,1000},PlotLegends->{"Change of Entropy","Entropy"},AxesLabel->Automatic]
Plot[Sen[t],{t,2,1000}];


ISS=Table[isol[t],{t,1,1000}];
ISNS=Table[(ISS[[t+1]]-ISS[[t]]),{t,1,999}];


Ost=ListPlot[ISNS]


Show[Ist,Est,Ost]


(*Export["Entropy.dat",ISS]*)
