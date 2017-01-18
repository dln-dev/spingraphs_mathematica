(* ::Package:: *)

BeginPackage["spingraphs`"]
Unprotect@@Names["spingraphs`*"];
ClearAll@@Names["spingraphs`*"];

NullMatrix::usage="NullMatrix[n], returns square matrix of dimension n"
RandSquareMat::usage="RandSquareMat[n], returns square matrix of dimension n containing random reals from (0,1)"
Pauli::usage="Pauli[i, N, j], returns jth generalized Pauli matrix acting on site i in a space of size N"
GraphProduct::usage="GraphProduct[G1, G2], returns cartesian graph product of two graphs G1 and G2"
WeightedForm::usage="WeightedForm[M], changes zeros in matrix M to infinity"
WeightedFormDot::usage="WeightedFormDot[M], changes zeros followed by dots in matrix M to infinity"
UnweightedForm::usage="UnweightedForm[M], changes Inifinity in matrix M to zeros"
MakeUnified::usage="MakeUnified[M], changes all values in matrix M to 1, if they are not 0. DEPRECATED, USE unitize"

ApproxMatrixExp::usage="ApproxMatrixExp[A, N], calculates MatrixExp[A] as a sum up to Nth power"
ApproxNumMatrixExp::usage="ApproxNumMatrixExp[A, NN, p], calculates MatrixExp[A] numerically with precision p up to NNth power"

HamiltonGeneral::usage="HamiltonGeneral[N, P], creates general Hamiltonian for a chain of size N with periodic boundary if P=1, or open boundary if P=0. x-, y-, and z-interaction strength is governed by Jx, Jy, Jz, local potential can be added via Hz. Lambda is 1 if not set otherwise"
HamiltonXX::usage"HamiltonXX[N, P], creates XX Hamiltonian for a chain of size N, periodic boundary if P=1, open boundary for P=0"

g2::usage="basic 2 node graph";
g3::usage="basic 3 node graph";
SpinStar::usage="SpinStar[N], creates spin star network with one central spin and N arms";
Star::usage="Star[h1, h2, h3, h4], creates spin star network with 3 arms and given local potentials. Central spin is the first";
WeightedSpinStar::usage="WeightedSpinStar[N], creates spin star with N arms and local potentials, using Hz for the potentials";
StarSquare::usage="StarSquare, returns the squared spin star switch with local potentials";
StarSquareWeighted::usage="StarSquareWeighted, returns the squared spin star switch with local potentials in weighted form";
FindDisconnectedGraphParts::usage="FindDisconnectedGraphParts[M], returns M^vertices. If matrix is a block-matrix, the graph is disconnected";

NikoAdjMat::usage="NikoAdjMat[mat], create matrix from a permutation matrix mat such that its matrix exponential recreates the permutation matrix mat";
NikoGraph::usage="NikoGraph[mat], creates graph such that the matrix exponential of its adjacency matrix yields the permutation matrix mat";

Begin["`Private`"]

lambda=1;
Jx[i_,N_]:=lambda/2 Sqrt[i (N-i)];
Jy[i_,N_]:=lambda/2 Sqrt[i (N-i)];
Jz[i_,N_]:=0*i+0*N;
Hz[i_,N_]:=0(*i+2*);

NullMatrix[n_]:=IdentityMatrix[n]-IdentityMatrix[n];
RandSquareMat[n_]:=RandomReal[{0,1},{n,n}];
Pauli[i_,N_,j_]:=(Ret={1};Do[Ret=KroneckerProduct[Ret,IdentityMatrix[2]],{k,1,i}];Ret=KroneckerProduct[Ret,PauliMatrix[j]];Do[Ret=KroneckerProduct[Ret,IdentityMatrix[2]],{l,i+1,N-1}];Ret);
GraphProduct[G1_,G2_]:=AdjacencyGraph[KroneckerProduct[AdjacencyMatrix[G1],IdentityMatrix[VertexCount[G2]]]+KroneckerProduct[IdentityMatrix[VertexCount[G1]],AdjacencyMatrix[G2]]];
WeightedForm[M_]:=M /. 0->Infinity;
WeightedFormDot[M_]:=M /. 0.->Infinity;
UnweightedForm[M_]:=M/. Infinity->0;
MakeUnified[M_]:=M/. !0->1;

ApproxMatrixExp[A_,N_]:=(Result=IdentityMatrix[Dimensions[A]];Do[Result+=MatrixPower[-1 A,n]/Factorial[n],{n,1,N}];Return[Result]);
ApproxNumMatrixExp[A_,NN_,p_]:=(Result=IdentityMatrix[Dimensions[A]];Do[Result+=MatrixPower[N[-A,p],n]/Factorial[n],{n,1,NN}];Return[Result]);

HamiltonGeneral[N_,P_]:=1/2 Sum[Jx[i+1,N]Pauli[i,N,1].Pauli[i+1,N,1]+Jy[i+1,N]Pauli[i,N,2].Pauli[i+1,N,2]+Jz[i+1,N]Pauli[i,N,3].Pauli[i+1,N,3]+2Hz[i,N] Pauli[i,N,3],{i,0,N-2}]+1/2 P (Jx[N-1,N]Pauli[N-1,N,1].Pauli[0,N,1]+Jy[N-1,N]Pauli[N-1,N,2].Pauli[0,N,2]+Jz[N-1,N]Pauli[N-1,N,3].Pauli[0,N,3]);
HamiltonXX[N_,P_]:=1/2Sum[Pauli[i,N,1].Pauli[i+1,N,1]+Pauli[i,N,2].Pauli[i+1,N,2],{i,0,N-2}]+1/2 P (Pauli[N-1,N,1].Pauli[0,N,1]+Pauli[N-1,N,2].Pauli[0,N,2]);

g2:=Graph[{1<->2}];
g3:=Graph[{5<->6,6<->7}];
SpinStar[N_]:=Graph[Table[0<->i,{i,N}]];
Star[h1_,h2_,h3_,h4_]:={{h1,1.41,1,1},{1.41,h2,0,0},{1,0,h3,0},{1,0,0,h4}};
WeightedSpinStar[N_]:=SpinStar[4]+Sum[Hz[i,N] Pauli[i,N,3],{i,0,N-1}];
StarSquare=WeightedAdjacencyMatrix[WeightedAdjacencyGraph[WeightedFormDot[KroneckerProduct[IdentityMatrix[4],Star[1.84,-2.61,0.766,0.766]]+KroneckerProduct[Star[1.84,-2.61,0.766,0.766],IdentityMatrix[4]]]]];
StarSquareWeighted=WeightedFormDot[KroneckerProduct[IdentityMatrix[4],Star[1.84,-2.61,0.766,0.766]]+KroneckerProduct[Star[1.84,-2.61,0.766,0.766],IdentityMatrix[4]]];
FindDisconnectedGraphParts[M_]:=MatrixPower[M,Dimensions[M]];
(*FindMatrixBlocks[M_]:=*)

skewStarMat={{0,0,1,1,0,0},{0,0,1,1,0,0},{1,1,0,0,-1,-1},{1,1,0,0,1,1},{0,0,-1,1,0,0},{0,0,-1,1,0,0}};
skewStar=WeightedAdjacencyGraph[WeightedForm[skewStarMat],{VertexLabels->"Name",EdgeLabels->"EdgeWeight"}];
skewStarMat2={{0,0,1,1,0,0},{0,0,1,1,0,0},{1,1,0,0,1,-1},{1,1,0,0,-1,1},{0,0,1,-1,0,0},{0,0,-1,1,0,0}};
skewStar2=WeightedAdjacencyGraph[WeightedForm[skewStarMat2],{VertexLabels->"Name",EdgeLabels->"EdgeWeight"}];
skewStarMat3={{0,-1,1,0,0,0},{-1,0,0,1,1,1},{1,0,0,-1,1,1},{0,1,-1,0,0,0},{0,1,1,0,0,0},{0,1,1,0,0,0}};
skewStar3=WeightedAdjacencyGraph[WeightedForm[skewStarMat3],{VertexLabels->"Name",EdgeLabels->"EdgeWeight"}];
skewStarMat4={{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}};
skewStar4=WeightedAdjacencyGraph[WeightedForm[skewStarMat4],{VertexLabels->"Name",EdgeLabels->"EdgeWeight"}];
skewStarMatUnity={{0,0,1,1,0,0},{0,0,1,1,0,0},{1,1,0,0,1,1},{1,1,0,0,1,1},{0,0,1,1,0,0},{0,0,1,1,0,0}};
skewStarUnity=WeightedAdjacencyGraph[WeightedForm[skewStarMatUnity],{VertexLabels->"Name",EdgeLabels->"EdgeWeight"}];

NikoAdjMat[mat_]:=(MAT=NullMatrix[Dimensions[mat][[1]]];Do[MAT+=Arg[Eigensystem[mat][[1]][[i]]]*Outer[Times,Eigensystem[mat][[2]][[i]],Eigensystem[mat][[2]][[i]]],{i,1,Dimensions[mat][[1]]}];MAT/Pi(*+IdentityMatrix[Dimensions[mat][[1]]]*));
NikoGraph[mat_]:=WeightedAdjacencyGraph[WeightedForm[NikoAdjMat[mat]],VertexLabels->"Name",EdgeLabels->"EdgeWeight"];


BroadenEdges[pts_List, e_]:=Block[{s=0.015,weight=PropertyValue[{HEC4,e},EdgeWeight]},{Arrowheads[{{s,0.1},{s,0.9}}],AbsoluteThickness[EdgeWeight*1.5],Arrow[pts]}];



End[]
Protect@@Names["spingraphs`*"];
EndPackage[]
