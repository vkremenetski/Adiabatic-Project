#include "itensor/all.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
using namespace std;
using namespace itensor;

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

template <typename MPS, typename Real>
std::vector<std::pair<MPS,Real>> zipp(std::vector<MPS> A, std::vector<Real> B){
    auto zipped = std::vector<std::pair<MPS,Real>>(0);
    for(int i = 0; i<A.size(); i++){
        zipped.push_back(std::make_pair(A.at(i),B.at(i)));
    }

    return zipped;
}

bool helper3(std::pair<MPS, Real> p1, std::pair<MPS, Real> p2){
    return (p1.second < p2.second);
}

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
        ampo += s*couplingWeight*4, "Sz", j, "Sz", k;
        ampo += (s-1)*2, "Sx", j;
    }
    ampo += s*weights[0]*2, "Sz", spec1;
    ampo += s*weights[1]*2, "Sz", spec2;
    //Print(ampo);
    MPO Ham = MPO(ampo);
    //Print(gapEDMPO(Ham));
    return Ham;
}


/*In the event that the states given are not the true ground/1st excited states, but the given states are believed to span these eigenstates,
 * input the given states and the Hamiltonian, and the method will return an ITensor with the eigenstates expressed as a linear combination of
 * the given states.
 */
ITensor getTrueEigenstates(MPO Ham, std::vector<MPS> givenStates){
    int dim = givenStates.size();
    Index i1 = Index("i1",dim);
    // Index i2 = Index("i2",dim);
    Index i2 = prime(i1);
    ITensor Ham2 = ITensor(i1,i2);
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            MPS state1 = givenStates.at(i);
            MPS state2 = givenStates.at(j);
            Ham2.set(i1(i+1), i2(j+1), overlap(state1, Ham, state2));
        }
    }
    ITensor U,d;
    diagHermitian(Ham2, U, d);
    return U;
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
        auto sweeps = Sweeps(22);
        sweeps.maxm() = 50,50,100,100,400;
        sweeps.noise() = 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10, 3e-11, 0;
        sweeps.cutoff() = 1E-11;
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
    auto sweeps = Sweeps(22);
    sweeps.maxm() = 50,50,100,100,200,400;
    sweeps.cutoff() = 1E-11;
    sweeps.noise() = 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10, 3e-11, 0;
    MPO Ham = getHam(N,p,w,sites,s);
    auto E0 = dmrg(psi0,Ham,sweeps,{"Quiet=",true});
    EnStates.at(0) = psi0;
    En.at(0) = E0;
    for(int i=1; i<2; i++){
        MPS psiI = MPS(sites);
        Real Ei = dmrg(psiI,Ham,EnStates,sweeps,{"Quiet=",true,"Weight=",20});
        EnStates.push_back(psiI);
        En.push_back(Ei);
    }
    auto EnS = zipp(EnStates, En);
    ITensor U = getTrueEigenstates(Ham, EnStates);
    IndexSet is = U.inds();
    Index i1 = is.index(1);
    Index i2 = is.index(2);
    auto a = U.cplx(i1(1),i2(1));
    auto b = U.cplx(i1(1),i2(2));
    auto GS = (a*EnStates.at(0)).plusEq(b*EnStates.at(1)); 
    Real entropy = maxEntropy2(N,p,w,sites,s);
    std::sort(EnS.begin(),EnS.end(), helper3);
    Real gap = EnS[1].second - EnS[0].second;
    printfln("DMRG GS: ", EnS[0].second);
    printfln("Our GS: ", overlap(GS, Ham, GS));
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
    for(Real s=0.5; s<=0.8; s+= step){
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
    for(Real s=0.68; s<0.69+step; s+= step){
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
    if(UpperQubitNumber >= 4){
        for(int n = 4; n<= UpperQubitNumber; n++){
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
//Finds the ground states of the Hamiltonian at time1 and time2, then prints out the overlap intermediate ground states have with these two.
//
void overlapToText(string title, Real time1, Real time2, int N, std::vector<int> positions, std::vector<Real> weights, Real step){
    ofstream myfile;
    myfile.open(title);
    auto sites = SpinHalf(N);
    auto psi1 = MPS(sites);
    auto psi2 = MPS(sites);
    auto sweeps = Sweeps(26);
    sweeps.maxm() = 50,50,100,100,200,300,500,1000,2000;
    sweeps.cutoff() = 1E-12;
    sweeps.noise() = 1, 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10, 3e-11, 1e-11, 3e-12, 1e-12, 0;
    MPO Ham1 = getHam(N,positions, weights, sites, time1);
    MPO Ham2 = getHam(N,positions, weights, sites, time2);
    Real t1 = 0.6835;
    Real t2 = 0.6855;
    dmrg(psi1, Ham1, sweeps, "Quiet");
    dmrg(psi2, Ham2, sweeps, "Quiet");
    for(Real s = t1; s<= t2; s += step){
        auto middlePsi = MPS(sites);
        MPO middleHam = getHam(N,positions,weights,sites,s);
        //Real middleEntropy = maxEntropy2(N,positions,weights,sites,s);
        dmrg(middlePsi, middleHam, sweeps, "Quiet");
        auto innProd1 = overlap(middlePsi, psi1);
        auto innProd2 = overlap(middlePsi, psi2);
        myfile << s << " " << innProd1 << " " << innProd2;
        /*auto complexInProd1 = overlapC(middlePsi, psi1);
        auto complexInProd2 = overlapC(middlePsi, psi2);
        auto approximatePsi = (complexInProd1 * psi1).plusEq(complexInProd2*psi2);
        normalize(approximatePsi);
        Real middleEntropy = maxEntropy(middlePsi);
        Real approximateEntropy = maxEntropy(approximatePsi);
        myfile << " " << middleEntropy << " " << approximateEntropy;*/
        myfile << "\n";
    }
    /*Real leftEntropy = maxEntropy2(N,positions,weights,sites,time1);
    Real rightEntropy = maxEntropy2(N,positions,weights,sites,time2);
    printfln("Left-Hand Entropy: ", leftEntropy);
    printfln("Right-Hand Entropy: ", rightEntropy);*/
}

int main(int argc, char* argv[]) {
    int N = 7;
    int mypositions[] = {1,N/2+1,N/2,N/2+1};
    Real myweights[] = {0.75,-1,-1,-0.5};
    std::vector<int> positions(mypositions,mypositions+4);
    std::vector<Real> weights(myweights,myweights+4);
    SpinHalf spins = SpinHalf(N);
    auto x = gapAndEntropy(N,positions, weights, spins, 0.5);
    //overlapToText("Overlap14q.txt", 0.65, 0.75, N, positions, weights, 0.00002);
    //timeToText("ThirteenQubitEvolution2.txt",N,positions,weights,0.0002);
    //qubitCountToText("NQubitEvolution.txt",16,weights,0.01);
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
