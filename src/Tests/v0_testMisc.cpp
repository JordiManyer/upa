
#include <iostream>
#include "myMath.h"
#include "debugIO.h"

using namespace std;
using namespace upa;

int main() {

    for (int i = 1; i < 5; ++i) {
        double w[i], z[i];
        GLquad(i,-1.0,1.0,w,z);

        cout << " > Order = " << i << endl;
        cout << "   > Weights = " << endl;
        printArray(i,w);
        cout << "   > Coords = " << endl;
        printArray(i,z);
        cout << endl;
    }


}

