#include "itensor/all.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
using namespace itensor;
AutoMPO getHamAMPO(int N, int spec1, int spec2, int breakPoint, Real weightFactor){
    auto sites = SpinHalf(N);
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; j++){
        ampo += weightFactor,"Sz",j,"Sz",j+1;
    }
    return ampo;
}
int main(int argc, char* argv[]) {
    printfln("Stuff");
}
