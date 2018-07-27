#include "itensor/all.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
using namespace itensor;
/* Returns MPO for Hamiltonian.
 * N refers to the number of qubits in the system.
 * 
 * positions is a vector containing the relevant position values. Namely:
 *   <First Special Qubit Site, Second ..., First Special Coupling, Second...>
 * 
 * weights is a vector containing the relevant weights. Namely:
 *   <First Special Qubit Weight, Second..., Weight of Normal Couplings, Weight of Special...>
 */

//Diagonalizes an MPO to find the gap.
Real gapEDMPO(MPO Ham) {
  ITensor prod(1.0);

  for(int i = 1; i <= Ham.N(); i++) {
    prod *= Ham.A(i);
  }

  ITensor U, D;

  diagHermitian(prod, U, D);

  Index ci = commonIndex(U, D);

  return (D.real(ci(ci.m() - 1), prime(ci)(ci.m() - 1)) - D.real(ci(ci.m()), prime(ci)(ci.m())));
}

MPO getHam(int N, std::vector<int> positions, std::vector<Real> weights, SpinHalf sites, Real s){
    int spec1 = positions[0];
    int spec2 = positions[1];
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N; j++){
        Real couplingWeight = weights[2];
        if(j == positions[2] or j == positions[3]){couplingWeight = weights[3];}
        int k = j+1;
        if(j==N){k = 1;}
        ampo += s*couplingWeight*4, "Sz", j, "Sz", k;
        ampo += (s-1)*2, "Sx", j;
    }
    ampo += s*weights[0]*2, "Sz", spec1;
    ampo += s*weights[1]*2, "Sz", spec2;
    Print(ampo);
    MPO Ham = MPO(ampo);

    Print(gapEDMPO(Ham));

    return Ham;
}

//Creates a Hamiltonian as a Single ITensor (method used for troubleshooting)
ITensor getTensorHam(int N, std::vector<int> positions, std::vector<Real> weights, SpinHalf sites, Real s){
    int spec1 = positions[0];
    int spec2 = positions[1];
    ITensor Ham;
    for(auto b: range1(N-1)){
        Real couplingWeight = weights[2];
        //if(b == positions[2] or b == positions[3]){couplingWeight = weights[3];}
        //int k = b+1;
        //if(b==N){k = 1;}
        auto term = s*couplingWeight * sites.op("Sx", b) * sites.op("Sx",b+1);
        //term += (s-1) * sites.op("S+",b) * sites.op("S-", b+1);
        //term += (s-1) * sites.op("S-",b) * sites.op("S+", b+1);
        //if(b == spec1){term += s * weights[0] * sites.op("Sz",b) * sites.op("Id",k);}
        //else if(b == spec2){term += s * weights[1] * sites.op("Sz", b) * sites.op("Id",k);}
        printfln("Term Rank 1: ", term.r());
        for(auto j: range1(b-1)){term *= sites.op("Id",j);}
        printfln("Term Rank 2: ", term.r());
        if(b <= N-1){
            for(auto j: range(b+2, N+1)){term *= sites.op("Id",j);}
        }
        printfln("Term Rank 3: ", term.r());
        printfln("Rank Ham: ", Ham.r());
        Ham += term;
        printfln("Just finished: ", b);
    }
    return Ham;
}

// Converts MPS to an ITensor
ITensor mpsToTensor(MPS matrixps){
    int N = matrixps.N();
    ITensor result = ITensor(1);
    for(int i = 1; i <= N; i++){
        result = result * matrixps.A(i);
    }
    return result;
}
ITensor mpoToTensor(MPO operator1){
    int N = operator1.N();
    ITensor result = ITensor(1);
    for(int i = 1; i<= N; i++){
        result = result * operator1.A(i);
    }
    return result;
}
//Returns the entropy of a given quantum state
Real stateEntropy(ITensor T){
    IndexSet indices = T.inds();
    Index ind1 = indices.index(1);
    Index ind2 = indices.index(2);
    ITensor U(ind1, ind2), S, V;
    svd(T,U,S,V);
    Index i1 = commonIndex(S,V);
    Index i2 = commonIndex(U,S);
    Real entropy = 0;
    for(int i = 1; i<= i1.m(); i++){
        Real l = S.real(i1(i),i2(i));
        entropy += -l*log(l);
    }
    return entropy;
}
/*Finds entanglement entropy of MPS state across "bond".
 * */
Real getEntropy(MPS state, int bond){
    state.position(bond);
    ITensor psi = state.A(bond)*state.A(bond+1);
    auto U = state.A(bond);
    ITensor S,V;
    auto decomposition = svd(psi,U,S,V);
    Real entropy = 0;
    for(auto lambda: decomposition.eigs()){
        entropy += -lambda*log(lambda);
    }
    return entropy;
}
/*Returns the size of the largest entanglement entropy
 * across any of the bonds of a linear MPS.
 */
Real maxEntropy(MPS state){
    Real maxEntropy = 0;
    for(int i = 1; i< state.N(); i++){
        Real entropy = getEntropy(state,i);
        if(entropy > maxEntropy){maxEntropy=entropy;}
    }
    return maxEntropy;
}
/* Returns the maximum entropy across all possible partitions of 
 * a ring MPS. */
Real maxEntropy2(int N, std::vector<int> p, std::vector<Real> w, SpinHalf sites, Real s){
    Real entropy = 0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j < 4; j++){
            p.at(j) += i;
            if(p.at(j) > N){p.at(j) = p.at(j) - N;}
        }
        MPO Ham = getHam(N,p,w,sites,s);
        MPS psi = MPS(sites);
        auto sweeps = Sweeps(13);
        sweeps.maxm() = 50,50,100,100,300;
        sweeps.noise() = 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 0;
        sweeps.cutoff() = 1E-10;
        dmrg(psi,Ham,sweeps,"Quiet");
        Real ent = maxEntropy(psi);
        if(ent > entropy){entropy = ent;}
    }
    return entropy;
}

/*Given the Sites and the Hamiltonian, returns the
 * gap between the ground and first excited state and the
 * entropy of the ground state in that order.
 */
std::vector<Real> gapAndEntropy(int N, std::vector<int> p, std::vector<Real> w, SpinHalf sites, Real s){
    auto psi0 = MPS(sites);
    auto EnStates = std::vector<MPS>(1);
    auto En = std::vector<Real>(1);
    auto sweeps = Sweeps(13);
    sweeps.maxm() = 50,50,100,100,200,300;
    sweeps.cutoff() = 1E-10;
    sweeps.noise() = 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 0;
    MPO Ham = getHam(N,p,w,sites,s);
    auto E0 = dmrg(psi0,Ham,sweeps,{"Quiet=",true});
    EnStates.at(0) = psi0;
    En.at(0) = E0;
    for(int i=1; i<=12; i++){
        MPS psiI = MPS(sites);
        Real Ei = dmrg(psiI,Ham,EnStates,sweeps,{"Quiet=",true,"Weight=",20});
        EnStates.push_back(psiI);
        En.push_back(Ei);
    }
    Real entropy = maxEntropy2(N,p,w,sites,s);
    std::sort(En.begin(),En.end());
    Real gap = En[1]-En[0];
    auto results = std::vector<Real>(2);
    results.at(0) = gap;
    results.at(1) = entropy;
    return results;
}
bool helper1(std::vector<Real> v1, std::vector<Real> v2){
    return (v1[0] < v2[0]);
}
bool helper2(std::vector<Real> v1, std::vector<Real> v2){
    return (v1[1] > v2[1]);
}
std::vector<Real> minGapAndEntropy(int N, std::vector<int> positions, std::vector<Real> weights, Real step){
    auto results = std::vector<std::vector<Real>>();
    auto sites = SpinHalf(N);
    for(Real s=0; s<=1.0; s+= step){
        results.push_back(gapAndEntropy(N, positions, weights, sites, s));
    }
    std::sort(results.begin(),results.end(), helper1);
    auto output = std::vector<Real>();
    output.push_back(results[0][0]);
    std::sort(results.begin(),results.end(), helper2);
    output.push_back(results[0][1]);
    return output;
}
//Creates title.txt file listing the entropy and energy gap over time
void timeToText(string title,int N, std::vector<int> positions, std::vector<Real> weights, Real step){
    ofstream myfile;
    myfile.open(title);
    auto sites = SpinHalf(N);
    for(Real s=0; s<=1; s+= step){
        auto results = gapAndEntropy(N, positions, weights, sites, s);
        myfile << s << " " << results.at(0) << " " << results.at(1);
        myfile << "\n";
    }
}

//creates title.txt file listing the maximum entropy and minimum gap as qubit number increases
//
void qubitCountToText(string title, int UpperQubitNumber, std::vector<Real> weights, Real step){
    ofstream myfile;
    myfile.open(title);
    if(UpperQubitNumber >= 6){
        for(int n = 6; n<= UpperQubitNumber; n++){
            int mypositions[] = {1,n/2+1,n/2,n/2+1};
            std::vector<int> positions(mypositions, mypositions+4);
            auto results = minGapAndEntropy(n,positions,weights,step);
            myfile << n << " " << results.at(0) << " "<< results.at(1);
            myfile << "\n";
            printfln("???????????????????????????????????????????????????????????????");
            printfln("This many qubits: ",n);
            printfln("???????????????????????????????????????????????????????????????");
        }
    }
}


int main(int argc, char* argv[]) {
    int N = 6;
    int mypositions[] = {1,N/2+1,N/2,N/2+1};
    Real myweights[] = {0.75,-1,-1,-0.5};
    std::vector<int> positions(mypositions,mypositions+4);
    std::vector<Real> weights(myweights,myweights+4);
    SpinHalf spins = SpinHalf(N);
    timeToText("SixQubitEvolution.txt",N,positions,weights,0.01);
    //qubitCountToText("NQubitEvolution.txt",16,weights,0.02);
    /*for(Real s = 0.75; s<= 1; s+=0.01){
        MPO Ham = getHam(N, positions, weights, spins, s);
        ITensor Hamiltonian = mpoToTensor(Ham);
        ITensor U, d;
        diagHermitian(Hamiltonian, U, d);
        IndexSet sett = d.inds();
        Index ind1 = sett.index(1);
        Index ind2 = sett.index(2);
        auto energies = std::vector<Real>();
        for(int i = 1; i<= ind1.m(); i++){
            energies.push_back(d.real(ind1(i),ind2(i)));
        }
        std::sort(energies.begin(),energies.end());
        printfln("Gap: ", energies.at(1)-energies.at(0));
    }*/

}
