(* ::Package:: *)

(* ::Input:: *)
(*Define symmetric quantization rho*)


CleanSlate;
$MinPrecision=MachinePrecision;
  \[Rho]z[z_] := z/(1 + Sqrt[1 - z])^2; 
 \[Rho]zb[z_] := \[Rho]z[Conjugate[z]]; 
 (*error estimates*)
    \[Rho]intErrorEstimateG[d_, DDs_, z_, \[Gamma]_] := (2^(\[Gamma] + 4*d)*\[Beta]^(\[Gamma] - 4*d)*Gamma[1 - \[Gamma] + 4*d, \[Beta]*DDs])/(Gamma[1 - \[Gamma] + 4*d]*(1 - r^2)^\[Gamma])/.  {r -> Abs[\[Rho]z[z]], \[Beta] -> -Log[Abs[\[Rho]z[z]]]}; 
      \[Rho]intErrorEstimateFt[d_, DDs_, z_, \[Gamma]_] := 
    extPrefac[d, 1-z] *\[Rho]intErrorEstimateG[d, DDs, z, \[Gamma]] + extPrefac[d, z] *\[Rho]intErrorEstimateG[d, DDs, 1 - z, \[Gamma]] ; 
     (*coefficients and prefactors*)
 extPrefac[\[CapitalDelta]\[Phi]_, x_] :=  ((x)*( Conjugate[x]))^\[CapitalDelta]\[Phi]; 
	kfunct[\[Beta]_, x_] :=  x^(\[Beta]/2)*Hypergeometric2F1[\[Beta]/2, \[Beta]/2, \[Beta], x]; 
    ConformalBlock[DD_, l_, z_] := ConformalBlock[DD, l, z] =((-1)^l/2^l)*((z*Conjugate[z])/(z - Conjugate[z]))*(kfunct[DD + l, z]*kfunct[DD - l - 2, Conjugate[z]] - 
       kfunct[DD + l, Conjugate[z]]*kfunct[DD - l - 2, z]); 

      
      (*random sample of z around (1/2+I0)*)
           Sample[nz_,var_,seed_] := Module[{imax}, SeedRandom[seed];Table[Abs[RandomVariate[NormalDistribution[0, var]]]+
           1/2+I Abs[RandomVariate[NormalDistribution[0, var]]],{imax,1,nz}]]; 
(*Exponential factor to renormalize*)
dimExpFactor[a_,b_]:=Log[4(1-Sqrt[1/2 -a -b])^2/(1/2 + a +b)];

renomFactor[dim_]:=Exp[-(dim-1)dimExpFactor[0,0]];
(*renomFactor[dim_]:=1;*)
(*Conformal Blocks*)
(*Generates block for one operator*)
qQGen[\[CapitalDelta]\[Phi]_,\[CapitalDelta]_,L_,zsample_]:=renomFactor[\[CapitalDelta]] (extPrefac[\[CapitalDelta]\[Phi], 1-zsample]    ConformalBlock[\[CapitalDelta], L , zsample]- extPrefac[\[CapitalDelta]\[Phi], zsample] ConformalBlock[\[CapitalDelta], L,1- zsample])2^(L);
(*Generates block for a list of operators*)
qQGenDims[\[CapitalDelta]\[Phi]_,\[CapitalDelta]L_,z_]:=qQGen[\[CapitalDelta]\[Phi],#1[[1]],#1[[2]], z]&/@\[CapitalDelta]L;
(*Generates block for the identity*)
qQId[\[CapitalDelta]\[Phi]_,zsample_]:=extPrefac[\[CapitalDelta]\[Phi], zsample] -extPrefac[\[CapitalDelta]\[Phi], 1-zsample]  ;

(*functionals*)
chi2Functional[qq0_,id_,w_,rhovec_]:=Block[{nu,s,r},
nu = Dimensions[w][[1]]-Length[rhovec];
r=(qq0.rhovec-id);
s=r.w.r;
Return[s/nu]];

chi2Functional[qq0_,id_,w_]:=Block[{nu,s,r,rhovec=(cweightedLeastSquares[qq0,id,w])[[1]];},
nu = Dimensions[w][[1]]-Length[rhovec];
r=(qq0.rhovec-id);
s=r.w.r;
Return[s/nu]];

logDetFunctional[qq0_,id_]:=Block[{pp},
pp=Join[qq0,id];
Return[Log[Det[pp]^2]]];

(*assorted functions*)
deltaFree[n_]:={2#,2#-2}&/@Range[1,n,1];
opeFreeRen[n_]:=(renomFactor[2#])^(-1) 2((2#-2)!)^2/(2(2#-2))!&/@Range[1,n,1];
selectiveMinors[mat_,elems_]:=
Det[mat[[elems]]];
spinAppender[\[CapitalDelta]_]:=Transpose[{\[CapitalDelta],Range[0,2Length[\[CapitalDelta]] -2,2]}];

(*MC routine*)

(*\[CapitalDelta]\[Phi]0_ Initial external dimension, \[CapitalDelta]LOriginal_Initial spectrum, Ndit_Number of Iterations, prec_precision, betad_1/Temperature, seed_ ,sigmaMC_sigma for the MC step, dcross_regularization radius, 
lmax_order of the highest operator to vary (redundant at this stage), idTag_string for identifying runs, initialOps_ This tells the routine whether to vary all operators from the begining or just let the newly added one vary on its own for a certain number of steps,
opsToVary_array with the *)

MetroGoFixedSelectiveDir[\[CapitalDelta]\[Phi]0_,deltaExtMax_,\[CapitalDelta]LOriginal_,Ndit_,prec_,betad_,seed_,sigmaMC_,dcross_,lmax_,idTag_,initialOps_,opsToVary_,sigmaz_,Nz_,elems_]:=
Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, nAccept=0,nReject=0,Ident,OPEcoeff,ActionTot,  TotD ,DDldataold,QQold,\[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L,dw,smearedaction,\[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]0,\[CapitalDelta]\[Phi]old,nops=Length[\[CapitalDelta]LOriginal]}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Join[Sample[Nz[[1]],sigmaz[[1]],seed],Sample[Nz[[2]],sigmaz[[2]],seed+1]]; 
  Print[Dimensions[zsample]];
Idsample = qQId[\[CapitalDelta]\[Phi], zsample];
Print[Dimensions[Idsample]];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
    Print[Dimensions[QQ0]];
     
    (*Monte Carlo Iteration*)
TotD =   Reap[ Do[
$MinPrecision=prec;
(*Save previous values*)
          \[CapitalDelta]LOld=\[CapitalDelta]L;
          QQold=QQ0;  
          \[CapitalDelta]\[Phi]old=\[CapitalDelta]\[Phi];
(*putative Unitarity bound
If[\[CapitalDelta]L[[1,1]]<1,\[CapitalDelta]L[[1,1]]=\[CapitalDelta]L[[1,1]]+1/2];*)
(*let every successive run start by varying only the new operator*)
        If[it<Ndit/10&& nops!=initialOps, 
        dimToVary=opsToVary[[-1]],  
        dimToVary = opsToVary[[RandomInteger[{1,Length[opsToVary]}]]]];
       (*Shift one dimension by a random amount*)       
         If[dimToVary==0,
         (*Vary external*)
         \[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]+ RandomVariate[NormalDistribution[0,sigmaMC/2]];
	 If[Abs[\[CapitalDelta]\[Phi]- \[CapitalDelta]\[Phi]0] > deltaExtMax, \[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]old;];
         Idsample = qQId[\[CapitalDelta]\[Phi], zsample];    
         QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];,
         (*Vary exchanged*)
          \[CapitalDelta]L[[dimToVary,1]] = \[CapitalDelta]L[[dimToVary,1]]+ RandomVariate[NormalDistribution[0,sigmaMC]];
          If[\[CapitalDelta]L[[1,1]]<1,\[CapitalDelta]L[[1,1]]=\[CapitalDelta]L[[1,1]]+1/2];
          QQ0[[dimToVary]] = qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]],\[CapitalDelta]L[[dimToVary]][[2]],zsample];
          ];
(*Reevaluate coefficients*)
           
          
    (*Coefficients for LES and action thence*)
          PP= Join[QQ0,{Idsample}]; 
          Actionnew = Log[(selectiveMinors[PP//Transpose,#]&/@elems)^2]//Total; 
         
QQsave=QQ0;
(*Brot noch schmieren? *)
If[dcross!=0,
smearedaction=Reap[Table[
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]+dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[(selectiveMinors[PP//Transpose,#]&/@elems)^2]//Total ];
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]-dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[(selectiveMinors[PP//Transpose,#]&/@elems)^2]//Total ];
          QQ0=QQsave,
          {dimToVary,1,lmax}]];

 Actionnew =Actionnew +Total[smearedaction[[2]]//Flatten] 
 ];
         
          metcheck = Exp[(-betad)*(Actionnew - Action)];
          rr = RandomReal[{0, 1}];
          If[metcheck>rr, 
          Action = Actionnew; nAccept= nAccept+1,
          \[CapitalDelta]L=\[CapitalDelta]LOld;QQ0=QQold;\[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]old; nReject= nReject+1,
          Print["Error"]; Break[]];
          
$MinPrecision=10;
   dw=Join[{\[CapitalDelta]\[Phi]},\[CapitalDelta]L[[All,1]]];
          Sow[ {it, dw, N[Action,10]}],
     {it, 1, Ndit}]]; 
     Print[{nAccept,nReject}];
$MinPrecision=3;
      Export["Res-fixed_Param_Nit="<>ToString[Ndit]<>"Nz="<>ToString[Nz]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[betad,3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[dcross,3]]<>"seed="<>ToString[seed]<>"id="<>idTag<>".txt", TotD[[2]]];]

(*MC routine-Chi2*)
MetroGoFixedSelectiveDirChi2[\[CapitalDelta]\[Phi]0_,deltaExtMax_,\[CapitalDelta]LOriginal_,Nz_,Ndit_,prec_,betad_,seed0_,sigmaMC_,dcross_,lmax_,idTag_,sigmaz_,tol_,opsToVary_]:=Block[{itd, DDldata, sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample,  PP0, PP1, lr, nr, Errvect, Factor, Factor0, seed=seed0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,\[CapitalDelta]\[Phi]old,  
    IdsampleEx,zOPE,QQOPE,converge=False,Calc,coeffTemp,Ident,OPEcoeff,ActionTot,  TotD ,DDldataold,QQold,resultsOld,\[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L,dw,errSample,\[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]0,results,nzeros=Length[\[CapitalDelta]LOriginal],fracvio=100,nzerosnew,fracvionew,res,sigmamods=ConstantArray[1,Length[\[CapitalDelta]LOriginal]]}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = qQId[\[CapitalDelta]\[Phi],zsample];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
     
          errSample=\[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample,0];
    (*Monte Carlo Iteration*)
TotD =   Reap[ Do[
$MinPrecision=prec;
If[fracvio<=tol, 
If[converge,Print["Convergence!"];Break[],
Print["resetting seed"]; 
seed=seed+1;  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = qQId[\[CapitalDelta]\[Phi],zsample];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
     
          errSample=\[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample,0];
                    res=(results[[1]].QQ0-Idsample)/errSample;
nzeros=Count[results[[1]],0];
fracvio=Count[Abs[res]<1//Thread,False]/Nz;
          converge=True],converge=False];
          \[CapitalDelta]LOld=\[CapitalDelta]L;
          QQold=QQ0; 
          \[CapitalDelta]\[Phi]old=\[CapitalDelta]\[Phi];
          resultsOld=results; 
(*This If avoids varying operators which still don't enter the OPE*)
       (*Shift one dimension by a random amount*)       
            If[Length[opsToVary]-nzeros<=0,      dimToVary = opsToVary[[RandomInteger[{1,Length[opsToVary]}]]],  dimToVary = opsToVary[[RandomInteger[{1,Length[opsToVary]-nzeros}]]]];
       (*Shift one dimension by a random amount*)       
         If[dimToVary==0,
         (*Vary external*)
         \[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi](1+ RandomVariate[NormalDistribution[0,sigmaMC]]);
	 If[Abs[\[CapitalDelta]\[Phi]- \[CapitalDelta]\[Phi]0] > deltaExtMax, \[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]old;];
         Idsample = qQId[\[CapitalDelta]\[Phi], zsample];
         errSample=\[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample,0];
         QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];,
         (*Vary exchanged*)
          \[CapitalDelta]L[[dimToVary,1]] = \[CapitalDelta]L[[dimToVary,1]](1+ RandomVariate[NormalDistribution[0,sigmaMC]]);
          If[\[CapitalDelta]L[[1,1]]<1,\[CapitalDelta]L[[1,1]]=\[CapitalDelta]L[[1,1]]+1/2];
          QQ0[[dimToVary]] = qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]],\[CapitalDelta]L[[dimToVary]][[2]],zsample];];
           results=cweightedLeastSquares[(QQ0//Transpose)/errSample,Idsample/errSample,IdentityMatrix[Nz]];
           
(*Need to create two different sets of points for this to work. Leaving out for now
  zsample = Sample[Nz,sigmaz,seed+1]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];

errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
*)

res=(results[[1]].QQ0-Idsample)/errSample;
nzerosnew=Count[results[[1]],0];
fracvionew=Count[Abs[res]<1//Thread,False]/Nz;

(*Debugging*)
Print[nzerosnew];
Print[fracvionew];

          If[fracvionew<=fracvio && nzeros>=nzerosnew, Print["Accepted"]; fracvio = fracvionew;nzeros=nzerosnew;sigmamods[[dimToVary]] =sigmamods[[dimToVary]] (101/100),
          \[CapitalDelta]L=\[CapitalDelta]LOld;QQ0=QQold;\[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]old ;  results=resultsOld; sigmamods[[dimToVary]] =sigmamods[[dimToVary]] (99/100)];
$MinPrecision=10;
dw=Join[{\[CapitalDelta]\[Phi]},\[CapitalDelta]L[[All,1]]];
          Sow[ {it, dw,results[[1]], {nzeros,fracvio}}],
     {it, 1, Ndit}]]; 
$MinPrecision=3;
      Export["Res-chi_Param_Nit="<>ToString[Ndit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[betad,3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[dcross,3]]<>"seed="<>ToString[seed0]<>"Nz="<>ToString[Nz]<>"id="<>idTag<>".txt", TotD[[2]]];]

cweightedLeastSquares[qq0_,id_,w_]:=Block[{rhovec,nu,s,r,n=1,qq0bis},
rhovec=Inverse[Transpose[qq0].w.qq0].Transpose[qq0] . w.id;
nu = Dimensions[w][[1]]-Length[rhovec];
r=(qq0.rhovec-id);
s=r.w.r;
While[Or@@(rhovec<0//Thread),
qq0bis=qq0[[;;,1;;-1-n]];
rhovec=Inverse[Transpose[qq0bis].w.qq0bis].Transpose[qq0bis] . w.id;
nu = Dimensions[w][[1]]-Length[rhovec];
r=(qq0bis.rhovec-id);
s=r.w.r;n=n+1];
If[n>1,
Return[{Join[rhovec,ConstantArray[0,n-1]], Sqrt[s/nu]}],
Return[{rhovec, Sqrt[s/nu]}]
];
]
weightedLeastSquares[qq0_,id_,w_]:=Block[{rhovec,nu,s,r},
rhovec=Inverse[Transpose[qq0].w.qq0].Transpose[qq0] . w.id;
nu = Dimensions[w][[1]]-Length[rhovec];
r=(qq0.rhovec-id);
s=r.w.r;
Return[{rhovec,(Diagonal[Inverse[Transpose[qq0].w.qq0]])^(-1/2), Sqrt[s/nu]}];
]
metroReturnAvg[\[CapitalDelta]\[Phi]_,deltaExtMax_,prec_,nit_,\[Beta]_,\[CapitalDelta]L_,seed_,initialOps_,idtag_,sigmaMC_,opsToVary_,sigmaz_,nz_,elems_]:=Block[{data,exact=Join[{1},deltaFree[Length[\[CapitalDelta]L]][[;;,1]]]},
MetroGoFixedSelectiveDir[\[CapitalDelta]\[Phi],deltaExtMax,\[CapitalDelta]L,nit,prec,\[Beta],seed,sigmaMC,1/3,Length[\[CapitalDelta]L],ToString[Length[\[CapitalDelta]L]]<>idtag,initialOps,opsToVary,sigmaz,nz,elems];
data= Get["Res-fixed_Param_Nit="<>ToString[nit]<>"Nz="<>ToString[nz]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".txt"];
(*Export["Plot-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".pdf",ListPlot[Table[data[[All,2]][[All,i]],{i,1,Length[\[CapitalDelta]L]+1}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},\[CapitalDelta]L[[;;,2]]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
Export["rel-error-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,Length[\[CapitalDelta]L]+1}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},\[CapitalDelta]L[[;;,2]]],PlotRange->All,PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
Export["zoomed-rel-error-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,Length[\[CapitalDelta]L]+1}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotRange->{-1/10,1/10},PlotLegends->Join[{"ext"},\[CapitalDelta]L[[;;,2]]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
*)
{Mean[data[[nit-100;;nit,2]]],(*StandardDeviation[data[[nit-100;;nit,2]]],*)data[[-1]]}];

metroReturnAvgChi2[\[CapitalDelta]\[Phi]_,deltaExtMax_,prec_,nit_,Nz_,\[Beta]_,\[CapitalDelta]L_,seed_,initialOps_,idtag_,sigmaz_,sigmaMC_,tol_,opsToVary_]:=Block[{data,exact=Join[{1},deltaFree[Length[\[CapitalDelta]L]][[;;,1]]]},
MetroGoFixedSelectiveDirChi2[\[CapitalDelta]\[Phi],deltaExtMax,\[CapitalDelta]L,Nz,nit,prec,\[Beta],seed,sigmaMC,1/3,Length[\[CapitalDelta]L],idtag,sigmaz,tol,opsToVary];
data= Get["Res-chi_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"Nz="<>ToString[Nz]<>"id="<>idtag<>".txt"];
(*Export["rel-error-chi2_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,Length[\[CapitalDelta]L]}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->\[CapitalDelta]L[[;;,2]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
Export["zoomed-rel-error-chi2_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,Length[\[CapitalDelta]L]}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotRange->{-1/10,1/10},PlotLegends->\[CapitalDelta]L[[;;,2]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
{Mean[data[[nit-100;;nit,2]]],StandardDeviation[data[[nit-100;;nit,2]]],Mean[data[[nit-100;;nit,3]]],StandardDeviation[data[[nit-100;;nit,3]]]}*)
Return[{data[[-1,2]],data[[-1,3]]}];
];

(*Plotters*)
logdetPlotnAv[filename_]:=Block[{data,exact,numDims},
data= Get[filename];
numDims=Length[data[[1,2]]];
exact=Join[{1},deltaFree[numDims-1][[;;,1]]];
Export[filename<>"nat-plot"<>".pdf",ListPlot[Table[data[[All,2]][[All,i]],{i,1,numDims}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},deltaFree[numDims-1][[;;,2]]],PlotLabel->filename]];
Export[filename<>"rel-error"<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,numDims}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},deltaFree[numDims-1][[;;,2]]],PlotLabel->filename]];
Export[filename<>"zoomed-error"<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,numDims}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotRange->{-1/10,1/10},PlotLegends->Join[{"ext"},deltaFree[numDims-1][[;;,2]]],PlotLabel->filename]];
{Mean[data[[-100;;-1,2]]],(*StandardDeviation[data[[nit-100;;nit,2]]],*)data[[-1]]}];

chi2PlotnAv[filename_]:=Block[{data,exact,numDims},
data= Get[filename];
numDims=Length[data[[1,2]]];
exact=Join[{1},deltaFree[numDims-1][[;;,1]]];
Export[filename<>"rel-error"<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,numDims}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},deltaFree[numDims-1][[;;,2]]],PlotLabel->filename]];
Export[filename<>"zoomed-rel-error"<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,numDims}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotRange->{-1/10,1/10},PlotLegends->Join[{"ext"},deltaFree[numDims-1][[;;,2]]],PlotLabel->filename]];
(*{Mean[data[[-100;;-1,2]]],StandardDeviation[data[[-100;;-1,2]]],Mean[data[[-100;;-1,3]]],StandardDeviation[data[[-100;;-1,3]]]}*)
Return[{data[[-1,2]],data[[-1,3]]}];
];

ccheckMetroWeightedBis[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_]:=Block[{itd, DDldata, sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,res,nzeros}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample =qQId[\[CapitalDelta]\[Phi],zsample];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=\[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample,0];
results=cweightedLeastSquares[(QQ0//Transpose)/errSample,Idsample/errSample,IdentityMatrix[Nz]];

  zsample = Sample[Nz,sigmaz,seed+1]; 
Idsample = qQId[\[CapitalDelta]\[Phi],zsample];

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];

nzeros=Count[results[[1]],0];

errSample=\[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1-nzeros,1]],zsample,0];
res=(results[[1]].QQ0-Idsample)/errSample;
Export["histogram-res-dist.pdf",Histogram[res,Round[Nz/50]]];
finalcheck=Abs[res]<1//Thread;
Return[{results,Count[finalcheck,True]/Nz, nzeros}];
];



mcIterator[\[CapitalDelta]\[Phi]_,initialOps_,finalOps_,\[CapitalDelta]Linitial_,\[Beta]_,nz_,prec_,seedO_,nits_,runid_,sigmaz_,sigmaMC_,maxReps_]:=Block[{\[CapitalDelta]L=\[CapitalDelta]Linitial,results,repCount=0,checks,it,seed=seedO,nzeros=finalOps},
it=initialOps;

SetOptions[RandomReal,WorkingPrecision->100];
results=Reap[
While[it<=finalOps,
\[CapitalDelta]L[[1;;it,1]]=Sow[metroReturnAvg[\[CapitalDelta]\[Phi],prec,nits[[it-initialOps+1]],\[Beta][[it-initialOps+1]],\[CapitalDelta]L[[1;;it]],seed+it,initialOps,runid,sigmaMC]][[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
Print[checks[[2;;3]]];
If[(checks[[2]]==1) &&(checks[[3]]<=nzeros+1) ,
nzeros=checks[[3]];it=it+1;repCount=0,
If[repCount<maxReps,\[CapitalDelta]L[[1;;it,1]]=\[CapitalDelta]L[[1;;it,1]](1+Table[RandomReal[{-1/100,1/100}],{i,1,it}]);
Print["Rejected"];seed=seed+finalOps;repCount=repCount+1,Print["Andato_a"];Break[]]];
];
];
Export["averages_n_checks"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt", results];
]

mcIteratorFullThing[\[CapitalDelta]\[Phi]0_,deltaExtMax_,initialOps_,finalOps_,\[CapitalDelta]Linitial_,\[Beta]_,nz_,prec_,seedO_,nits_,runid_,sigmaz_,sigmaMC_,maxReps_,sigmaChi_,opsToVary_,sigmazLogDet_,nzLogDet_,elems_]:=
Block[{\[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]0,\[CapitalDelta]L=\[CapitalDelta]Linitial,results,repCount=0,checks,it,seed=seedO,nzeros=finalOps,holdMyBeer},
it=initialOps;

SetOptions[RandomReal,WorkingPrecision->100];
results=Reap[
While[it<=finalOps,
holdMyBeer=Sow[metroReturnAvg[\[CapitalDelta]\[Phi],deltaExtMax,prec,nits[[it-initialOps+1]],\[Beta][[it-initialOps+1]],\[CapitalDelta]L[[1;;it]],seed+it,initialOps,runid,sigmaMC,opsToVary[[it-initialOps+1]],sigmazLogDet[[it-initialOps+1]],nzLogDet[[it-initialOps+1]],elems[[it-initialOps+1]]]][[1]];
\[CapitalDelta]L[[1;;it,1]]=holdMyBeer[[2;;-1]];
\[CapitalDelta]\[Phi]=holdMyBeer[[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
Print[checks[[2;;3]]];
holdMyBeer=Sow[metroReturnAvgChi2[\[CapitalDelta]\[Phi],deltaExtMax,prec,nits[[it-initialOps+1]],nz,1,\[CapitalDelta]L[[1;;it]],seed+2it,initialOps,runid,sigmaz,sigmaChi[[it-initialOps+1]],0,opsToVary[[it-initialOps+1]]]][[1]];
\[CapitalDelta]L[[1;;it,1]]=holdMyBeer[[2;;-1]];
\[CapitalDelta]\[Phi]=holdMyBeer[[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
Print[checks[[2;;3]]];

If[(checks[[2]]==1) &&(checks[[3]]<=nzeros+1) ,
nzeros=checks[[3]];it=it+1;repCount=0,
If[repCount<maxReps,\[CapitalDelta]L[[1;;it,1]]=\[CapitalDelta]L[[1;;it,1]](1+Table[RandomReal[{-1/100,1/100}],{i,1,it}]);
Print["Rejected"];seed=seed+finalOps;repCount=repCount+1,Print["Andato_a"];Break[]]];
];
];
Export["averages_n_checks"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt", results];
]

(*no check baby*)
mcIteratorNoCheck[\[CapitalDelta]\[Phi]_,initialOps_,finalOps_,\[CapitalDelta]Linitial_,\[Beta]_,nz_,prec_,seedO_,nits_,runid_,sigmaz_,sigmaMC_,maxReps_]:=Block[{\[CapitalDelta]L=\[CapitalDelta]Linitial,results,repCount=0,checks,it,seed=seedO,nzeros=finalOps},
it=initialOps;
SetOptions[RandomReal,WorkingPrecision->100];
results=Reap[
While[it<=finalOps,
\[CapitalDelta]L[[1;;it,1]]=Sow[metroReturnAvg[\[CapitalDelta]\[Phi],prec,nits[[it-initialOps+1]],\[Beta][[it-initialOps+1]],\[CapitalDelta]L[[1;;it]],seed+it,initialOps,runid,sigmaMC]][[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
If[True,
nzeros=checks[[3]];it=it+1;repCount=0,
If[repCount<maxReps,\[CapitalDelta]L[[1;;it,1]]=\[CapitalDelta]L[[1;;it,1]](1+Table[RandomReal[{-1/100,1/100}],{i,1,it}]);
Print["Rejected"];Print[checks[[2;;3]]];seed=seed+finalOps;repCount=repCount+1,Print["Andato_a"];Break[]]];
];
];
Export["averages_n_checks"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt", results];
]

mcRefiner[\[CapitalDelta]\[Phi]_,nreps_,\[CapitalDelta]Linitial_,\[Beta]_,nz_,prec_,seedO_,nits_,runid_,sigmaz_,sigmaMC_]:=Block[{\[CapitalDelta]L=\[CapitalDelta]Linitial,results,repCount=0,checks,it,seed=seedO},
it=1;
SetOptions[RandomReal,WorkingPrecision->100];
results=Reap[
While[it<=nreps,
\[CapitalDelta]L[[;;,1]]=Sow[metroReturnAvg[\[CapitalDelta]\[Phi],prec,nits,\[Beta],\[CapitalDelta]L,seed,Length[\[CapitalDelta]L],runid,sigmaMC[[it]]]][[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L,prec,seed+1,nz,sigmaz]];
it=it+1;
];
];
Export["averages_n_checks-refiner"<>ToString[nreps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt", results];
]
mcPlotRef[nreps_,nz_,prec_,seed_,runid_]:=Block[{data,dims,opes,mcDimserr,mcDims,mcOpes,refDims,refOpes},
data= Get["averages_n_checks-refiner"<>ToString[nreps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt"];
dims=data[[1,Range[1,2nreps -1,2]]];
opes=data[[1,Range[2,2nreps ,2]]];
mcDims=ListPlot[dims[[;;,1,1;;4]]//Transpose,PlotLegends->{"l=0","l=2","l=4","l=6"},PlotRange->All];
mcDimserr=Table[ListPlot[((dims[[;;,1,i]]) - 2i)/2i,PlotRange->{-0.1,0.1}],{i,1,4}];
refDims=Plot[{2,4,6,8},{x,0,nreps},PlotStyle->Dashed];
Export["dims-plot-refiner"<>ToString[nreps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".pdf",Show[mcDims,refDims],PlotLabel->"Dims_from-refiner"<>ToString[nreps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]];
Table[Export["dims-err-plot-refiner"<>ToString[nreps]<>runid<>ToString[i]<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".pdf",Show[mcDimserr[[i]]],PlotLabel->"dim_error_from-refiner"<>ToString[nreps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]],{i,1,4}];
]


fullMC[debugging_,\[CapitalDelta]\[Phi]_,initialOps_,finalOps_,\[CapitalDelta]Linitial_,\[Beta]_,nz_,prec_,seedO_,nits_,runid_,sigmaz_,sigmaMC_,maxReps_,tol_]:=Block[{\[CapitalDelta]L=\[CapitalDelta]Linitial,results,repCount=0,checks,it,seed=seedO,nzeros=finalOps,logdetConv=True},
it=initialOps;
SetOptions[RandomReal,WorkingPrecision->100];
results=Reap[
While[it<=finalOps,
\[CapitalDelta]L[[1;;it,1]]=Sow[metroReturnAvg[\[CapitalDelta]\[Phi],prec,nits[[it-initialOps+1]],\[Beta][[it-initialOps+1]],\[CapitalDelta]L[[1;;it]],seed+it,initialOps,runid,sigmaMC[[1]]]][[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
If[debugging,Print[checks]];
If[(checks[[2]]==1) &&(checks[[3]]<=nzeros+1) ,
nzeros=checks[[3]];it=it+1;repCount=0,
If[repCount<maxReps,\[CapitalDelta]L[[1;;it,1]]=\[CapitalDelta]L[[1;;it,1]](1+Table[RandomReal[{-1/100,1/100}],{i,1,it}]);
Print["Rejected"];Print[checks[[2;;3]]];seed=seed+finalOps;repCount=repCount+1,
Print["Failed logdet mc"];logdetConv=False;Break[]]];
];
];
Export["averages_n_checks"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt", results];
If[nzeros<finalOps-initialOps&&logdetConv,
metroReturnAvgChi2[prec,nits[[1]],nz,1,\[CapitalDelta]L,seed,initialOps,runid,sigmaz,sigmaMC[[2]],tol],
Print["Bad logdet. Skipping nov min."]];
]

mcPlotDimsAndOPEs[initialOps_,finalOps_,nz_,prec_,seed_,runid_]:=Block[{data,dims,opes,mcDimserr,mcDims,mcOpes,refDims,refOpes,nrecs=(finalOps-initialOps) +1},
data= Get["averages_n_checks"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".txt"];
dims=data[[1,Range[1,2nrecs -1,2]]];
opes=data[[1,Range[2,2nrecs ,2]]];
mcDims=ListPlot[dims[[;;,1,1;;4]]//Transpose,PlotLegends->{"l=0","l=2","l=4","l=6"},PlotRange->All];
mcDimserr=Table[ListPlot[((dims[[;;,1,i]]) - 2i)/2i,PlotRange->{-0.1,0.1}],{i,1,4}];
refDims=Plot[{2,4,6,8},{x,0,nrecs},PlotStyle->Dashed];
mcOpes=ListPlot[opes[[;;,1,1,1;;4]]//Transpose,PlotLegends->{"l=0","l=2","l=4","l=6"},PlotRange->All];
refOpes=Plot[{2,1/3},{x,0,nrecs},PlotStyle->Dashed];
Export["dims-plot"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".pdf",Show[mcDims,refDims],PlotLabel->"Dims_from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]];
Export["opes-plot"<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".pdf",Show[mcOpes,refOpes],PlotLabel->"OPEs_from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]];
Table[Export["dims-err-plot"<>ToString[i]<>"from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]<>"nz="<>ToString[nz]<>".pdf",Show[mcDimserr[[i]]],PlotLabel->"dim_error_from"<>ToString[initialOps]<>"to"<>ToString[finalOps]<>runid<>"prec="<>ToString[prec]<>"seed="<>ToString[seed]],{i,1,4}];
Return[opes[[;;,2]]]
]

gradientComparisonLog[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,epsilon_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,gradientLog,func0}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;Length[\[CapitalDelta]L]
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
func0=logDetFunctional[QQ0,{Idsample}];
gradientLog=(Table[
    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]L[[i,1]]+epsilon; QQ0[[i]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[i]][[1]],\[CapitalDelta]L[[i]][[2]],zsample];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]LOriginal[[i,1]];
logDetFunctional[QQ0,{Idsample}],{i,1,Length[\[CapitalDelta]L]}]-func0)/(epsilon);
Return[{func0,gradientLog}];

];

gradientComparisonLogSmeared[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,epsilon_,dcross_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,gradientLog,func0}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
    QQsave=QQ0;
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
func0=logDetFunctional[QQ0,{Idsample}];
smearedaction=Reap[Table[
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]+dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]-dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];QQ0=QQsave;,{dimToVary,1,Length[\[CapitalDelta]L]}]];
          
 func0 =func0 +Total[smearedaction[[2]]//Flatten]; 
gradientLog=Table[
    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]L[[i,1]]+epsilon; QQ0[[i]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[i]][[1]],\[CapitalDelta]L[[i]][[2]],zsample];
   Actionnew= logDetFunctional[QQ0,{Idsample}];
smearedaction=Reap[Table[
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]+dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]-dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];QQ0=QQsave;,{dimToVary,1,Length[\[CapitalDelta]L]}]];
\[CapitalDelta]L[[i,1]]=\[CapitalDelta]LOriginal[[i,1]];
 Actionnew =Actionnew +Total[smearedaction[[2]]//Flatten]; 
(Actionnew-func0)/(epsilon),{i,1,Length[\[CapitalDelta]L]}];


Return[{func0,gradientLog}];

];

gradientComparisonChi[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,epsilon_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,gradient,func0}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
results=weightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];
func0=chi2Functional[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)],results[[1]]];
gradient=(Table[ QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]L[[i,1]]+epsilon; QQ0[[i]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[i]][[1]],\[CapitalDelta]L[[i]][[2]],zsample];results=weightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]LOriginal[[i,1]];
chi2Functional[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)],results[[1]]],{i,1,Length[\[CapitalDelta]L]}]-func0)/(epsilon);
Return[{func0,gradient}];
];

gradientComparisonChiSmeared[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,epsilon_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,gradient,func0}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
results=weightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];
func0=chi2Functional[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)],results[[1]]];
gradient=(Table[ QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]L[[i,1]]+epsilon; QQ0[[i]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[i]][[1]],\[CapitalDelta]L[[i]][[2]],zsample];results=weightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];\[CapitalDelta]L[[i,1]]=\[CapitalDelta]LOriginal[[i,1]];
chi2Functional[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)],results[[1]]],{i,1,Length[\[CapitalDelta]L]}]-func0)/(epsilon);
Return[{func0,gradient}];
];


(*PLotters*)
genLog[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,dcross_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,gradientLog,func0}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = qQId[\[CapitalDelta]\[Phi],zsample];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  
    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
    QQsave=QQ0;
       Actionnew= logDetFunctional[QQ0,{Idsample}];
       If[dcross!=0,
smearedaction=Reap[Table[
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]+dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]-dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];QQ0=QQsave;,{dimToVary,1,Length[\[CapitalDelta]L]}]];
 Actionnew =Actionnew +Total[smearedaction[[2]]//Flatten]; 
 ];
Return[Actionnew];
];


biDimChi2DependencePlotter[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,range_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,deltaArray}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
results=cweightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];
ContourPlot[\[CapitalDelta]L[[1,1]]=2+x;\[CapitalDelta]L[[2,1]]=4+y;QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];chi2Functional[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)],results[[1]]],{x,-range,range},{y,-range,range},PlotLegends->Automatic]
]
biDimChi2DependenceDataGen[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,range_,npoints_,opstovary_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,deltaArray}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
results=cweightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];
ParallelTable[{\[CapitalDelta]L[[opstovary[[1]],1]]=\[CapitalDelta]LOriginal[[opstovary[[1]],1]]+x range/npoints, \[CapitalDelta]L[[opstovary[[2]],1]]=\[CapitalDelta]LOriginal[[opstovary[[2]],1]]+ y  range/npoints, QQ0[[opstovary]] = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[opstovary]],zsample];
results=cweightedLeastSquares[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)]];chi2Functional[(QQ0//Transpose),Idsample,DiagonalMatrix[errSample^(-2)],results[[1]]]},{x,-npoints,npoints},{y,-npoints,npoints}]
];

(*This function is analogous to ccheckMetroWeightedBis, but it takes an arbitrary set of OPE coefficients instead of solving for them.*)
manCheck[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,rhovec_]:=Block[{itd, DDldata, sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini,
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,res},
    (*precision*)
    SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed];
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  
Export["histogram-res-dist.pdf",Histogram[res,Round[Nz/50]]];
    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],0],{i,1,Nz}];
res=(rhovec.QQ0-Idsample)/errSample;
finalcheck=Abs[res]<1//Thread;
Return[{Sqrt[(res.res)/(Nz-Length[rhovec])],If[And@@finalcheck,And@@finalcheck,finalcheck]}];
];

(*REsidual plotter*)
resFunc[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,rhovec_,prec_,x_,y_]:=Block[{itd, DDldata, sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini,
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,res},
    (*precision*)
    SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

  zsample = x + I y;
Idsample = SetPrecision[(zsample*Conjugate[zsample])^\[CapitalDelta]\[Phi] -
        ((1 - zsample)*(1 - Conjugate[zsample]))^\[CapitalDelta]\[Phi],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample= \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample,0];
res=(rhovec.QQ0-Idsample)/errSample;
Return[res];
];
biDimChi2DependenceDataGen[\[CapitalDelta]\[Phi]_,\[CapitalDelta]LOriginal_,prec_,seed_,Nz_,sigmaz_,opstovary_,x_,y_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
    zsample, Idsample, PP0, PP1, lr, nr, Errvect, Factor, Factor0, ppm, DDldataEx, PPEx, QQEx, Idsampleold, ip, nvmax, QQFold,  
    \[CapitalDelta]LOld,dimToVary,PP,QQsave,\[CapitalDelta]L=\[CapitalDelta]LOriginal,dw,smearedaction,\[Rho],rhovec,eqs,rhosol,last,check,results,indices,rhopos,meanrho,sigmarho,finalcheck,errSample,deltaArray}, 
    (*precision*)
SetOptions[{RandomReal,RandomVariate,NSolve},WorkingPrecision->prec];
$MaxPrecision=prec;
$MinPrecision=prec;

    SeedRandom[seed];
  zsample = Sample[Nz,sigmaz,seed]; 
Idsample = SetPrecision[Table[(zsample[[zv]]*Conjugate[zsample[[zv]]])^\[CapitalDelta]\[Phi] -
        ((1 - zsample[[zv]])*(1 - Conjugate[zsample[[zv]]]))^\[CapitalDelta]\[Phi], {zv, 1, Nz}],prec];
    \[CapitalDelta]L = \[CapitalDelta]LOriginal;
  \[CapitalDelta]L[[All,1]] = SetPrecision[\[CapitalDelta]L[[All,1]],prec];
  

    QQ0 = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample];
errSample=Table[ \[Rho]intErrorEstimateFt[\[CapitalDelta]\[Phi],\[CapitalDelta]LOriginal[[-1,1]],zsample[[i]],1],{i,1,Nz}];
\[CapitalDelta]L[[opstovary[[1]],1]]=\[CapitalDelta]LOriginal[[opstovary[[1]],1]]+x; \[CapitalDelta]L[[opstovary[[2]],1]]=\[CapitalDelta]LOriginal[[opstovary[[2]],1]]+ y ;
 QQ0[[opstovary]] = qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[opstovary]],zsample];
results=cweightedLeastSquares[(QQ0//Transpose)/errSample,Idsample/errSample,IdentityMatrix[Nz]];
chi2Functional[(QQ0//Transpose)/errSample,Idsample/errSample,IdentityMatrix[Nz],results[[1]]]
];


(* Fullthing Tester
*)

minops=4;
maxops=7;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
\[CapitalDelta]L[[2,1]]=4;
nits=2 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 300;
sigmaChiList=Table[1/1000,{i,minops,maxops}];
opsToVary=Table[Drop[Range[0,opa],{3}],{opa,minops,maxops}];
sigmaz=Table[{opa/20,opa/20},{opa,minops,maxops}];
Nz=Table[{5,opa},{opa,minops,maxops}];
elems=Table[Table[RandomSample[Range[5+opa],opa+1],{i,1,opa}],{opa,minops,maxops}];
ParallelTable[
\[Beta]list=Table[1/((2i-6+temps)(i)),{i,minops,maxops}];
mcIteratorFullThing[11/10,1/4,minops,maxops,\[CapitalDelta]L,\[Beta]list,135,100,29+10temps,nits,"testingtesting-external-temps="<>ToString[temps],1/10,1/10,0,sigmaChiList,opsToVary,sigmaz,Nz,elems],{temps,1,4}]


(* Plotting 

files = {
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00333sigmaMC=0.100dcross=0seed=6357id=print-test-600-5-t={8}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00417sigmaMC=0.100dcross=0seed=6353id=print-test-600-5-t={6}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=3000Nz={700, 100}prec=100beta=0.00417sigmaMC=0.100dcross=0seed=16id=print-test-800-5-t={5}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00556sigmaMC=0.100dcross=0seed=6349id=print-test-600-5-t={4}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00833sigmaMC=0.100dcross=0seed=6345id=print-test-600-5-t={2}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00476sigmaMC=0.100dcross=0seed=6351id=print-test-600-5-t={5}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00370sigmaMC=0.100dcross=0seed=6355id=print-test-600-5-t={7}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.00667sigmaMC=0.100dcross=0seed=6347id=print-test-600-5-t={3}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=4000Nz={590, 10}prec=100beta=0.0111sigmaMC=0.100dcross=0seed=6343id=print-test-600-5-t={1}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=3000Nz={700, 100}prec=100beta=0.00500sigmaMC=0.100dcross=0seed=14id=print-test-800-5-t={4}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=3000Nz={700, 100}prec=100beta=0.00625sigmaMC=0.100dcross=0seed=12id=print-test-800-5-t={3}.txt",
"datsync/uluviano/Res-fixed_Param_Nit=3000Nz={700, 100}prec=100beta=0.00357sigmaMC=0.100dcross=0seed=18id=print-test-800-5-t={6}.txt"
};
logdetPlotnAv[#]&/@files

*)



(* multipoint tester

ops=4;
a=
ParallelTable[
nz=20;
\[CapitalDelta]L=deltaFree[ops];
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1+ 2/10);
elemss=Table[Range[1+5j-5,5j],{j,1,nz/5}];
MetroGoFixedSelectiveDir[1,\[CapitalDelta]L,2000,100,1/(2(2+t)),61+2t,1/10,0,4,"print-test-800-5-t="<>ToString[{t}],4,{1,2,3,4},{1/5,1/5},{nz-10,10},elemss]//Timing ,{t,1,6}];
Print[a];
*)

(* Landscape fun
\[CapitalDelta]L=deltaFree[4];
aaa10=ParallelTable[\[CapitalDelta]L[[1,1]]=y/100;{x/100,y/100,Total@Table[genLog[x/100,\[CapitalDelta]L,100,i,5,1/5,0], {i, 1, 100}]},{x,1,200},{y,100,300}];


Export["landscape-nosmear.csv",Flatten[aaa10,1]]


Export["landscape-nosmear-500-5.pdf",ListPlot3D[Flatten[aaa10,1]]]



ListPlot3D[Flatten[aaa10,1]] *)
