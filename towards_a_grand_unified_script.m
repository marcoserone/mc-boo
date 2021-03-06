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
 extPrefac[\[CapitalDelta]\[Phi]_, x_] := extPrefac[\[CapitalDelta]\[Phi], x] = ((x)*( Conjugate[x]))^\[CapitalDelta]\[Phi]; 
	kfunct[\[Beta]_, x_] := kfunct[\[Beta], x] = x^(\[Beta]/2)*Hypergeometric2F1[\[Beta]/2, \[Beta]/2, \[Beta], x]; 
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

(*MC routine*)
(*\[CapitalDelta]\[Phi]0_ Initial external dimension, \[CapitalDelta]LOriginal_Initial spectrum, Ndit_Number of Iterations, prec_precision, betad_1/Temperature, seed_ ,sigmaMC_sigma for the MC step, dcross_regularization radius, 
lmax_order of the highest operator to vary (redundant at this stage), idTag_string for identifying runs, initialOps_ This tells the routine whether to vary all operators from the begining or just let the newly added one vary on its own for a certain number of steps,
opsToVary_array with the *)
MetroGoFixedSelectiveDir[\[CapitalDelta]\[Phi]0_,\[CapitalDelta]LOriginal_,Ndit_,prec_,betad_,seed_,sigmaMC_,dcross_,lmax_,idTag_,initialOps_,opsToVary_,sigmaz_,Nz_,elems_]:=Block[{itd, DDldata,  sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
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
        If[it<Ndit/10&& nops!=initialOps,dimToVary=Nz-1,  dimToVary = opsToVary[[RandomInteger[{1,Length[opsToVary]}]]]];
       (*Shift one dimension by a random amount*)       
         If[dimToVary==0,
         (*Vary external*)
         \[CapitalDelta]\[Phi]=\[CapitalDelta]\[Phi]+ RandomVariate[NormalDistribution[0,sigmaMC/2]];
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
smearedaction=Reap[Table[
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]+dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[(selectiveMinors[PP//Transpose,#]&/@elems)^2]//Total ];
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]-dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[(selectiveMinors[PP//Transpose,#]&/@elems)^2]//Total ];
          QQ0=QQsave,
          {dimToVary,1,lmax}]];

 Actionnew =Actionnew +Total[smearedaction[[2]]//Flatten] ;
         
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
      Export["Res-fixed_Param_Nit="<>ToString[Ndit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[betad,3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[dcross,3]]<>"seed="<>ToString[seed]<>"id="<>idTag<>".txt", TotD[[2]]];]

(*MC routine-Chi2*)
MetroGoFixedSelectiveDirChi2[\[CapitalDelta]\[Phi]0_,\[CapitalDelta]LOriginal_,Nz_,Ndit_,prec_,betad_,seed0_,sigmaMC_,dcross_,lmax_,idTag_,sigmaz_,tol_,opsToVary_]:=Block[{itd, DDldata, sigmaD, Action=100000000, Actionnew=0, Action0, DDldatafixed, QQ0, QQ1, str, Lmax, Nvmax, rr, metcheck, sigmaDini, 
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
      Export["Res-chi_Param_Nit="<>ToString[Ndit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[betad,3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[dcross,3]]<>"seed="<>ToString[seed0]<>"Nz="<>ToString[Nz]<>"id="<>idTag<>".txt", TotD[[2]]];
      Return[converge];
      ]

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
metroReturnAvg[\[CapitalDelta]\[Phi]_,prec_,nit_,\[Beta]_,\[CapitalDelta]L_,seed_,initialOps_,idtag_,sigmaMC_,opsToVary_,sigmaz_,nz_,elems_]:=Block[{data,exact=Join[{1},deltaFree[Length[\[CapitalDelta]L]][[;;,1]]]},
MetroGoFixedSelectiveDir[\[CapitalDelta]\[Phi],\[CapitalDelta]L,nit,prec,\[Beta],seed,sigmaMC,1/3,Length[\[CapitalDelta]L],ToString[Length[\[CapitalDelta]L]]<>idtag,initialOps,opsToVary,sigmaz,nz,elems];
data= Get["Res-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".txt"];
Export["Plot-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".pdf",ListPlot[Table[data[[All,2]][[All,i]],{i,1,Length[\[CapitalDelta]L]+1}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},\[CapitalDelta]L[[;;,2]]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
Export["rel-error-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,Length[\[CapitalDelta]L]+1}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->Join[{"ext"},\[CapitalDelta]L[[;;,2]]],PlotRange->All,PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
Export["zoomed-rel-error-fixed_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>idtag<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-exact[[i]])/exact[[i]],{i,1,Length[\[CapitalDelta]L]+1}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotRange->{-1/10,1/10},PlotLegends->Join[{"ext"},\[CapitalDelta]L[[;;,2]]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
{Mean[data[[nit-100;;nit,2]]],(*StandardDeviation[data[[nit-100;;nit,2]]],*)data[[-1]]}];


metroReturnAvgChi2[\[CapitalDelta]\[Phi]_,prec_,nit_,Nz_,\[Beta]_,\[CapitalDelta]L_,seed_,initialOps_,idtag_,sigmaz_,sigmaMC_,tol_,opsToVary_]:=Block[{data},
MetroGoFixedSelectiveDirChi2[\[CapitalDelta]\[Phi],\[CapitalDelta]L,Nz,nit,prec,\[Beta],seed,sigmaMC,1/3,Length[\[CapitalDelta]L],idtag,sigmaz,tol,opsToVary];

data= Get["Res-chi_Param_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"Nz="<>ToString[Nz]<>"id="<>idtag<>".txt"];
Export["rel-error-chi2_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-2i)/(2i),{i,1,Length[\[CapitalDelta]L]}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotLegends->\[CapitalDelta]L[[;;,2]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
Export["zoomed-rel-error-chi2_Nit="<>ToString[nit]<>"prec="<>ToString[prec]<>"beta="<>ToString[N[\[Beta],3]]<>"sigmaMC="<>ToString[N[sigmaMC,3]]<>"dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]<>"id="<>ToString[Length[\[CapitalDelta]L]]<>".pdf",ListPlot[Table[(data[[All,2]][[All,i]]-2i)/(2i),{i,1,Length[\[CapitalDelta]L]}],Joined->True,GridLines->Automatic,PlotStyle->Thin,PlotRange->{-1/10,1/10},PlotLegends->\[CapitalDelta]L[[;;,2]],PlotLabel->ToString[Length[\[CapitalDelta]L]]<>"Nit="<>ToString[nit]<>" prec="<>ToString[prec]<>" beta="<>ToString[N[\[Beta],3]]<>" sigmaMC="<>ToString[N[1/10,3]]<>" dcross="<>ToString[N[1/3,3]]<>"seed="<>ToString[seed]]];
(*{Mean[data[[nit-100;;nit,2]]],StandardDeviation[data[[nit-100;;nit,2]]],Mean[data[[nit-100;;nit,3]]],StandardDeviation[data[[nit-100;;nit,3]]]}*)
Return[{data[[-1,2]],data[[-1,3]],converge}];
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

mcIteratorFancySchmancy[\[CapitalDelta]\[Phi]_,initialOps_,finalOps_,\[CapitalDelta]Linitial_,\[Beta]_,nz_,prec_,seedO_,nits_,runid_,sigmaz_,sigmaMC_,maxReps_,sigmaChi_]:=
Block[{\[CapitalDelta]L=\[CapitalDelta]Linitial,results,repCount=0,checks,it,seed=seedO,nzeros=finalOps,chi2Results},
it=initialOps;

SetOptions[RandomReal,WorkingPrecision->100];
results=Reap[
While[it<=finalOps,
\[CapitalDelta]L[[1;;it,1]]=Sow[metroReturnAvg[\[CapitalDelta]\[Phi],prec,nits[[it-initialOps+1]],\[Beta][[it-initialOps+1]],\[CapitalDelta]L[[1;;it]],seed+it,initialOps,runid,sigmaMC]][[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
Print[checks[[2;;3]]];
chi2Results=Sow[metroReturnAvgChi2[prec,2 nits[[it-initialOps+1]],nz,\[Beta][[1]],\[CapitalDelta]L[[1;;it]],seed+it,initialOps,runid,sigmaz,sigmaChi[[it-initialOps+1]],0]];
\[CapitalDelta]L[[1;;it,1]]=chi2Results[[1]];
checks=Sow[ccheckMetroWeightedBis[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[1;;it]],prec,seed+1,nz,sigmaz]];
Print[checks[[2;;3]]];

If[((checks[[2]]==1) &&(checks[[3]]<=nzeros+1)) ||chi2Results[[3]],
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
smearedaction=Reap[Table[
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]+dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];
           QQ0[[dimToVary]] =qQGen[\[CapitalDelta]\[Phi],\[CapitalDelta]L[[dimToVary]][[1]]-dcross,\[CapitalDelta]L[[dimToVary]][[2]],zsample];  PP = Join[QQ0, {Idsample}]; 
          Sow[ Log[Det[PP]^2]];QQ0=QQsave;,{dimToVary,1,Length[\[CapitalDelta]L]}]];
 Actionnew =Actionnew +Total[smearedaction[[2]]//Flatten]; 
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



spinAppender[\[CapitalDelta]_]:=Transpose[{\[CapitalDelta],Range[0,2Length[\[CapitalDelta]] -2,2]}];



opeFreeRen[9]//N


\[CapitalDelta]L=deltaFree[20];
\[CapitalDelta]L[[1;;4,1]]=\[CapitalDelta]L[[1;;4,1]] (1+ 1/2);
\[CapitalDelta]L[[4;;20,1]]=\[CapitalDelta]L[[4;;20,1]] (1+ 1/10);
nits=15{300,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100};
\[Beta]list={1/7,1/9,1/11,1/13,1/14,1/15,1/20,1/20,1/20,1/20,1/20,1/20,1/20,1/20,1/20,1/20,1/20};
mcIteratorNoCheck[1,4,20,\[CapitalDelta]L,\[Beta]list,20,100,123,nits,"let's_see_how_far_we_go",1/100,1/10,1]


\[CapitalDelta]L=deltaFree[20];
\[CapitalDelta]L[[;;,1]]=\[CapitalDelta]L[[;;,1]] (1+ 1/1000);
metroReturnAvg[1,100,20,1/20,\[CapitalDelta]L,123,20,"manyops",1/10]


\[CapitalDelta]L=deltaFree[20];
\[CapitalDelta]L[[;;,1]]=\[CapitalDelta]L[[;;,1]] (1- 1/1000);
metroReturnAvg[1,100,20,1/20,\[CapitalDelta]L,123,20,"manyops-neg",1/10]


\[CapitalDelta]L=deltaFree[5];
\[CapitalDelta]L[[;;,1]]=\[CapitalDelta]L[[;;,1]] (1+ 1/1000);
metroReturnAvg[1,100,20,1/20,\[CapitalDelta]L,123,20,"fewops",1/10]


\[CapitalDelta]L=deltaFree[5];
\[CapitalDelta]L[[;;,1]]=\[CapitalDelta]L[[;;,1]] (1- 1/1000);
metroReturnAvg[1,100,20,1/20,\[CapitalDelta]L,123,20,"fewops-neg",1/10]


\[CapitalDelta]L=deltaFree[5];
\[CapitalDelta]L[[;;,1]]=\[CapitalDelta]L[[;;,1]] (1- 1/1000);
metroReturnAvg[1,100,20,1/20,\[CapitalDelta]L,123,20,"fewops-neg",1/10]


minops=4;
maxops=10;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
mcIteratorFancySchmancy[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,100,100,564,nits,"testing-fancy",1/10,1/10,3]


minops=4;
maxops=10;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
sigmaChiList=Table[1/1000,{i,minops,maxops}];
sigmaChiList[[-3]]=10^(-7/2);
sigmaChiList[[-2]]=10^(-4);
sigmaChiList[[-1]]=10^(-4);
mcIteratorFancySchmancy[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,200,100,54,nits,"testing-fancy-autonomous-sigma",1/10,1/10,3,sigmaChiList]


minops=4;
maxops=10;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
sigmaChiList=Table[1/1000,{i,minops,maxops}];
sigmaChiList[[-3]]=10^(-7/2);
sigmaChiList[[-2]]=10^(-4);
sigmaChiList[[-1]]=10^(-9/2);
ParallelTable[
mcIteratorFancySchmancy[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,200,100,10+i,nits,"testing-seed"<>ToString[i],1/10,1/10,0,sigmaChiList]//Timing,{i,1,4}]


minops=4;
maxops=12;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
sigmaChiList=Table[1/1000,{i,minops,maxops}];
sigmaChiList[[-6]]=10^(-7/2);
sigmaChiList[[-5]]=10^(-4);
sigmaChiList[[-4]]=10^(-9/2);
sigmaChiList[[-3]]=10^(-5);
sigmaChiList[[-2]]=10^(-11/2);
sigmaChiList[[-1]]=10^(-6);
mcIteratorFancySchmancy[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,200,100,736,nits,"testing-conv"<>ToString[2],1/10,1/10,0,sigmaChiList]//Timing


minops=4;
maxops=12;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=16 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
sigmaChiList=Table[1/1000,{i,minops,maxops}];
sigmaChiList[[-6]]=10^(-7/2);
sigmaChiList[[-5]]=10^(-9/2);
sigmaChiList[[-4]]=10^(-9/2);
sigmaChiList[[-3]]=10^(-5);
sigmaChiList[[-2]]=10^(-11/2);
sigmaChiList[[-1]]=10^(-6);
ParallelTable[
mcIteratorFancySchmancy[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,200,100,999+i,nits,"testing-more-nits"<>ToString[i],1/10,1/10,0,sigmaChiList]//Timing,{i,1,4}]


minops=4;
maxops=12;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
sigmaChiList=Table[1/1000,{i,minops,maxops}];
sigmaChiList[[-6]]=10^(-7/2);
sigmaChiList[[-5]]=10^(-4);
sigmaChiList[[-4]]=10^(-9/2);
sigmaChiList[[-3]]=10^(-5);
sigmaChiList[[-2]]=10^(-11/2);
sigmaChiList[[-1]]=10^(-6);
ParallelTable[
mcIteratorFancySchmancy[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,200,100,734+i,nits,"testing-seed"<>ToString[i],1/10,1/10,0,sigmaChiList]//Timing,{i,1,4}]


minops=4;
maxops=10;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
mcIterator[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,100,100,9,nits,"limit-50-10-check-bis",1/50,1/10,3]


minops=4;
maxops=25;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
mcIterator[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,100,100,96,nits,"limit-50-10-check-bis",1/50,1/10,3]


minops=4;
maxops=50;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=16 (ConstantArray[100,maxops-minops+1]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-3),{i,minops,maxops}];
mcIteratorNoCheck[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,100,120,49,nits,"limit-50-10-check",1/50,1/10,3]


mcPlotDimsAndOPEs[4,30,100,120,49,"limit-50-10-check"]


minops=4;
maxops=25;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;minops,1]]=\[CapitalDelta]L[[1;;minops,1]] (1+ 1/2);
\[CapitalDelta]L[[minops+1;;maxops,1]]=\[CapitalDelta]L[[minops+1;;maxops,1]] (1+ 1/10);
nits=15 (ConstantArray[100,maxops-minops]);
nits[[1]] = 3000;
\[Beta]list=Table[1/(2i-2),{i,minops,maxops}];
mcIterator[1,minops,maxops,\[CapitalDelta]L,\[Beta]list,100,200,96,nits,"limit-50-10-check-highprec",1/50,1/10,3]


mcPlotDimsAndOPEs[4,40,20,100,499,"limit-50-30"]


maxops=4;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;maxops,1]]=\[CapitalDelta]L[[1;;maxops,1]] (1+ 1/2);
metroReturnAvg[1,100,3000,1/5,\[CapitalDelta]L,123,4,"testing-refine",1/10]



\[CapitalDelta]L[[1;;maxops,1]]={1.88590382856212932468624550126776117182679919757466050162807461805214006888373493768357294078838465`100.,4.063693655015287303110863867082401606483915566359442620710942313969832545774009536291364649118885776`100.,5.864603015065234000680756400806591219496727148075402434101686987019290634364885266817155054317919521`100.,7.642572459836518252507454327186402226028077417843381469262517437296138083036980874757984863806232652`100.}


\[CapitalDelta]L



metroReturnAvg[1,100,1000,1/5,\[CapitalDelta]L,123,4,"testing-refine",1/20]


\[CapitalDelta]L[[1;;maxops,1]]={1.914381497099365473944104879771293645925141151509430257383874397652935398936797598847095411652859041`100.,3.994432400887540321889568266594671506154399012833844287257440820635946280422457598961326793941662677`100.,5.986731131910766633638016450802184656216994464283898119051427379750289894673427396581942730426455941`100.,7.973643479599193229645221207893662195313760229508709201904318053971311848721894358446905692428530763`100.};


metroReturnAvg[1,100,1000,1/5,\[CapitalDelta]L,123,4,"more-refined",1/30]


\[CapitalDelta]L[[1;;maxops,1]]=metroReturnAvg[1,100,1000,1/5,\[CapitalDelta]L,123,4,"more-refined",1/40][[1]]


\[CapitalDelta]L[[1;;maxops,1]]=metroReturnAvg[1,100,1000,1/5,\[CapitalDelta]L,123,4,"more-refined",1/50][[1]]


\[CapitalDelta]L[[1;;maxops,1]]=metroReturnAvg[1,100,1000,1/5,\[CapitalDelta]L,123,4,"more-refined",1/60][[1]]


metroReturnAvg[1,100,100,1/5,\[CapitalDelta]L,123,4,"more-refined",1/70]


maxops=4;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;maxops,1]]=\[CapitalDelta]L[[1;;maxops,1]] (1+ 1/2);
sigmaMC=Table[1/(10i),{i,1,5}];
mcRefiner[1,5,\[CapitalDelta]L,1/5,40,100,5,1500,"repet-test",1/10,sigmaMC]


Table[mcPlotRef[10,40,100,10+i,"biggertest"<>ToString[i]],{i,1,6}]


maxops=4;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;maxops,1]]=\[CapitalDelta]L[[1;;maxops,1]] (1+ 6/10);
sigmaMC=Table[1/(5i),{i,2,10}];
ParallelTable[
mcRefiner[1,5,\[CapitalDelta]L,1/5,40,100,5+35i,3000,"init6-10"<>ToString[i],1/10,sigmaMC]//Timing,{i,10,15}]


maxops=4;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;maxops,1]]=\[CapitalDelta]L[[1;;maxops,1]] (1+ 1/2);


\[CapitalDelta]L


maxops=4;
\[CapitalDelta]L=deltaFree[maxops];
\[CapitalDelta]L[[1;;maxops,1]]=\[CapitalDelta]L[[1;;maxops,1]] (1+ 7/10);
sigmaMC=Table[1/(5i),{i,2,10}];
ParallelTable[
mcRefiner[1,5,\[CapitalDelta]L,1/5,40,100,5+35i,3000,"init6-10"<>ToString[i],1/10,sigmaMC]//Timing,{i,10,15}]


Table[mcPlotRef[5,40,100,5+35i,"init6-10"<>ToString[i]],{i,10,15}]





data=Get["averages_n_checksfrom4to25limit-50-10-check-highprecprec=200seed=196nz=100.txt"];


ListPlot[data[[1,Range[1,21,2],1]][[;;,1;;4]]//Transpose,GridLines->Automatic]


dims=data[[1,Range[1,21,2],1]];


\[CapitalDelta]L=spinAppender[{1.9311855622833880771515155329449895730209136881291167362071058293029744741007286656519072515812943559661543276377386720709374031217605451566884720901916649346700406581626452217663922983007089452549988`200., 3.9987804311486566923804811145989301809914029802042064191754102318776709439126885696111095975694581854912633288646631512870671641362137246501723898213930806790041371560867995189614646894802664258478625`200., 6.0022885200158557318398591589982473739596513159760237279234873420524360905387554128029104939152572143336684873754353890621961540076143266153713464813497552257167157668657056293542795970604329205769705`200., 8.0030850514644125011003440782124448872287234386220607018025562612733113151466469321768820412689341685281736731556401313757100422728964075331876934689484624563856925456910668242706835653039519709646858`200., 10.004509602684044539776799946007910448858527572649204661552421103675618238134551277043622383341518322310176826145356538603408703696282048638621466208883234189211153782905356313780046945772048497315235`200., 12.007335484505358053952027104992695991384822954285227590038781424670471536102357295987316443994196078613584646737707538985736728127557300618148210275306884413764439828810723878140469005121909534490246`200., 15.360101171519475928732428127840039145422958939086343310680668244092448190912677967407121736664535821930072152301178832267251536850882530482078909090412664329590620071977885278491309491171479139285041`200.}];


\[CapitalDelta]L//N


metroReturnAvgChi2[100,2000,1000,1,\[CapitalDelta]L,5,4,"chi2-test-8ops",1/10,1/1000,0]


metroReturnAvgChi2[100,2000,200,1,\[CapitalDelta]L,5,4,"dd",1/10,1/1000,0]


\[CapitalDelta]L={{2,0},{4,2},{6,0}};
metroReturnAvg[3/2,100,3000,1/5,\[CapitalDelta]L,13,3,"wtf",1/10]


f[x_,y_]:=genLog[x,{{y,0},{2x+2,2},{2x+2,0},{2x+4,2},{2x+4,0}},100,123,6,1/100];


ContourPlot[f[x,y],{x,3/2,7/4},{y,3-1/10,7/2+1/10}]


\[CapitalDelta]L=deltaFree[4];
Plot[genLog[x,\[CapitalDelta]L,100,123,5,1/100,1/3]+genLog[x-1/3,\[CapitalDelta]L,100,123,5,1/100,1/3]+genLog[x+1/3,\[CapitalDelta]L,100,123,5,1/100,1/3],{x,0,2},PlotRange->All]


ops=10;
\[CapitalDelta]L=deltaFree[ops];
Plot[genLog[x,\[CapitalDelta]L,100,123,ops+1,1/100,1/3],{x,0,2},PlotRange->All]


ops=10;
\[CapitalDelta]L=deltaFree[ops];
Plot[genLog[x,\[CapitalDelta]L,100,123,ops+1,1/10,1/3],{x,0,2},PlotRange->All]


Plot[genLog[x,\[CapitalDelta]L,100,123,11,1/100,1/3]+genLog[x-1/3,\[CapitalDelta]L,100,123,11,1/100,1/3]+genLog[x+1/3,\[CapitalDelta]L,100,123,11,1/100,1/3],{x,1-1/5,1+1/5},PlotRange->All]


ops=8;
\[CapitalDelta]L=deltaFree[ops];
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1);
\[CapitalDelta]L[[2,1]]=4;
metroReturnAvg[(130)/100,100,500,1/3,\[CapitalDelta]L,433,ops,"var_smaller-d"<>ToString[115],1/10,{0}][[1,1]]


ops=8;
\[CapitalDelta]L=deltaFree[ops];
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1);
\[CapitalDelta]L[[2,1]]=4;
Table[
metroReturnAvg[(132+2j)/100,100,500,1/3,\[CapitalDelta]L,433i,ops,"var_smaller-d"<>ToString[j],1/10,{0}][[1,1]],{i,1,4},{j,0,5}]


\[CapitalDelta]0={{1.789774094691055416085529247810831323634803587320569253152758301727872153808068519554052882773873738`100.,1.78906299321541646433362000219275977989315697492224272518790706626029979003093571762174208369391017`100.,1.788842351699641197037721584831857505577383663988898270867234073333014294938735682079599022853748759`100.,1.793281294902012941701857872897648723248627903342977864207166108512837436936395629662461406443717471`100.,1.789873036222850863582178281188195223895123385603781499674320445281756023345182506272434756236556298`100.,1.789195919974158784937073008772110484880971662220071918288978961073554384232245562380744960438935039`100.},{1.793804212979880956690018788112996930448806134632435209043743890173328937574973956820928021129020984`100.,1.791311001323192335176181722494074778026141669139086280120402117444725829055244564832259461276735362`100.,1.79271379304105709084386035662175351260556546109175488823005177596017100870556570577542167693713696`100.,1.791442524569980996364179379811378413569783687963511068727648929889027438545005114231601451769498199`100.,1.638228260771591931665884903318126722998038356954315866340134314472711753442499875246493437840217608`100.,1.791520338174535915281087771139244626310259778947695554893435807574017483415215594335926844880739714`100.},{1.789851256435161827103170835335120359259039011765238575164457640245909367163546052851743223539935377`100.,1.789597497896943759864323965725982271696659273794395215371341880608890695668218599919084318892566362`100.,1.789136651427685885305716764368566892329838323544553534825196676504329235383527792671175885193274359`100.,1.789841315620759076372732147494827614618932305194031928340890623610680168323979085254656723839919101`100.,1.790295752076915085340523011610165859403971447749729859182864388153901556956489099152489114637437449`100.,1.789811479930626157427540419456578400159047537702608761208733475295320485409543259373294278813999269`100.},{1.788574198528809137602901607511492756958976159712745458913045928796111491752970047348911190304304988`100.,1.793134015341897396982923418618256270141500887178714558068203288282874892524353765915852456591590014`100.,1.637933576774804640434707061358077725516557193715150458823344346862987932580643230718540618384615747`100.,1.637977792602548164448111010179055217512667144840695508235722800937052896464472040506860237774111342`100.,1.791715616875078407065550069009368494228165045152980321789016009508230536923104039964919033433020447`100.,1.791316294051550615112614411820260455186484825079029319699611518920884756512010534371266771686637668`100.}}//Mean//Mean;


ops=4;
\[CapitalDelta]L=deltaFree[ops];
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1+ 3/10);
metroReturnAvg[10/10,100,3000,1/5,\[CapitalDelta]L,439,ops,"testing_double_sigma-ref",1/10,{1,2,3,4},{1/100,1/100},{4,1}]


 ops=4;
\[CapitalDelta]L=deltaFree[ops];[[1]]
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1+ 1/10);
\[CapitalDelta]L[[2,1]]=4;
 metroReturnAvgChi2[11/10,100,1000,500,1,\[CapitalDelta]L,1562,ops,"amoaver2",1/10,10^(-3),0,{0,1,3,4},{1/100,2}]





\[CapitalDelta]L=deltaFree[4];


\[CapitalDelta]L=deltaFree[4];
aaa10=ParallelTable[\[CapitalDelta]L[[1,1]]=y/100;{x/100,y/100,genLog[x/100,\[CapitalDelta]L,100,123,5,1/100,1/3]+genLog[x/100,\[CapitalDelta]L,100,123,5,2,1/3]},{x,1,200},{y,100,300}];


aaa//Dimensions[[1]]


bbb//Dimensions


Export["ext_and_scalar_landscp-smeared.csv",Flatten[aaa10,1]]


ListPlot3D[Flatten[aaa10,1]]





[[1]]


qQGenDims[\[CapitalDelta]\[Phi],\[CapitalDelta]L,zsample[[1]]]


faaark=RandomReal[{0,1},{10,10}];
faaarkbig=RandomReal[{0,1},{11,10}];
Det[faaark]//RepeatedTiming
Minors[faaarkbig,10]//Total


Minors[faaarkbig,10]^2//Log//Total


{{2.827355419131852082964114183416335089310789687512488402307342289013031588115624378931341791119212712`100.*^-53,3.831670875874290370299877927232665741066506192914136493301558490891361784841220851418953028164786144`100.*^-30,7.98506799079893431447538238417605904920585671581430110175104303464835936117703413561060982134533011`100.*^-29,6.182490880300291110931720233967866559666495658516723204073460978136001817428486630964715980251487478`100.*^-27,1.319573107314472843736889381910167425402807404812831074350082265092456621195347584094569472167610807`100.*^-28,1.311028157370531245002419647381856104966058148278138303095969036711783954145598489890925056470736092`100.*^-32}}//Flatten//Total


ops=6;
\[CapitalDelta]L=deltaFree[ops];
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1+ 1/10);
\[CapitalDelta]L[[2,1]]=4;
elemss=Table[RandomSample[Range[13],7],{i,1,7}];
metroReturnAvg[115/100,100,4000,1/(49),\[CapitalDelta]L,43,ops,"in_search_of_ext",1/10,{0,1,3,4,5,6},{1/100,3/2},{10,3},elemss]











ops=4;
\[CapitalDelta]L=deltaFree[ops];(*
\[CapitalDelta]L[[1;;ops,1]]=\[CapitalDelta]L[[1;;ops,1]] (1+ 5/100);
\[CapitalDelta]L[[2,1]]=4;*)
elemss=Table[RandomSample[Range[20],5],{i,1,7}];
metroReturnAvg[120/100,100,1000,1/(49),\[CapitalDelta]L,69,ops,"in_search_of_ext",1/100,{0},{1/100,3/2},{12,8},elemss]



