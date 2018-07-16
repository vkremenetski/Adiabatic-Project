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
MPO getHam(int N, std::vector<int> positions, std::vector<Real> weights, SpinHalf sites, Real s){
    int spec1 = positions[0];
    int spec2 = positions[1];
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N; j++){
        Real couplingWeight = weights[2];
        if(j == positions[2] or j == positions[3]){couplingWeight = weights[3];}
        int k = j+1;
        if(j==N){k = 1;}
        ampo += s*couplingWeight, "Sz", j, "Sz", k;
        ampo += (s-1), "Sx", j;
    }
    ampo += s*weights[0], "Sz", spec1;
    ampo += s*weights[1], "Sz", spec2;
    MPO Ham = MPO(ampo);
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
        auto sweeps = Sweeps(5);
        sweeps.maxm() = 50,50,100,100,300;
        sweeps.cutoff() = 1E-9;
        sweeps.noise() = 0.2;
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
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-9;
    sweeps.noise() = 0.15;
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
    for(Real s=0.3+step; s<=0.5; s+= step){
        results.push_back(gapAndEntropy(N, positions, weights, sites, s));
    }
    std::sort(results.begin(),results.end(), helper1);
    auto output = std::vector<Real>();
    output.push_back(results[0][0]);
    std::sort(results.begin(),results.end(), helper2);
    output.push_back(results[0][1]);
    return output;
}
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
    int N = 8;
    int mypositions[] = {1,N/2+1,N/2,N/2+1};
    Real myweights[] = {3,-4,-4,-2};
    std::vector<int> positions(mypositions,mypositions+4);
    std::vector<Real> weights(myweights,myweights+4);
    //timeToText("EightQubEvo.txt",N,positions,weights,0.01);
    qubitCountToText("NQubitEvolution.txt",18,weights,0.01);
}
