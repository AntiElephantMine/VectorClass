/* A few notes about dealing with doubles

    1 i. The double compare function needs to be changed. It is not a transitive equality operation. A method to
         fix this would be to snap the doubles on to a grid and return true if two doubles snap on to the same section
         of the grid. I don't know how to implement this just yet - wait until I've read more about float comparisons.
    2. Doubles can get stored as negative zero, so adding +.0 when outputting the vector prevents displaying "-0".

*/

#include "Vector.h"

using namespace std;

const double epsilon = 1e-6; //double tolerance

bool double_equals(double a, double b){
    return abs(a-b) < epsilon;
}

double approximate(double a, double b){
    return double_equals(a,b) ? b : a;
}

Vector<3> cross_product(const Vector<3>& lhs, const Vector<3>& rhs){
    double newX = (lhs[1]*rhs[2]) - (lhs[2]*rhs[1]);
    double newY = (lhs[2]*rhs[0]) - (lhs[0]*rhs[2]);
    double newZ = (lhs[0]*rhs[1]) - (lhs[1]*rhs[0]);

    return Vector<3>({newX, newY, newZ});
}
