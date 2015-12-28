#ifndef MATHSVectors
#define MATHSVectors

#include <string>
#include <iostream>
#include <array>
#include <initializer_list>
#include <stdexcept> //for the exceptions
#include <cmath> //sqrt, sin, cos, abs
#include <algorithm>
#include <numeric> //accumulate algorithm


//FORWARD DECLARATIONS:

//main templates
template <unsigned N> class Vector;
template <unsigned N> std::ostream& operator<< (std::ostream&, const Vector<N>&);
template <unsigned N> std::istream& operator>> (std::istream&, Vector<N>&);
template <unsigned N> bool operator==(const Vector<N>&, const Vector<N>&);
template <unsigned N> bool operator!=(const Vector<N>&, const Vector<N>&);
template <unsigned N> Vector<N> operator+(const Vector<N>&, const Vector<N>&);
template <unsigned N> Vector<N> operator-(const Vector<N>&, const Vector<N>&);
template <unsigned N> Vector<N> operator*(const Vector<N>&, double);
template <unsigned N> Vector<N> operator*(double, const Vector<N>&);
template <unsigned N> Vector<N> operator/(const Vector<N>&, double);
template <unsigned N> double dot_product(const Vector<N>&, const Vector<N>&);
template <unsigned N> Vector<N> normalised(const Vector<N>&); //radians
template <unsigned N> Vector<N> rotated_vector(const Vector<N>&, size_t, size_t, double); //radians

//in Vector.cpp
bool double_equals(double, double);
double approximate(double, double);
Vector<3> cross_product(const Vector<3>&, const Vector<3>&);




//CLASS DEFINITION:

template <unsigned N> class Vector{
    friend bool operator== <N>(const Vector&, const Vector&); /*Friendship granted here so base-array can be used for algorithms*/
    friend double dot_product<N>(const Vector&, const Vector&);
public:
    Vector() = default;
    template <typename... Args> Vector(Args...); /*if we use Args&& the this is the better match for copy-constructor. So use either const Args&&, const Args&, Args&, or Args.*/

    explicit operator bool() const;

    double& operator[](size_t p);
    const double& operator[](size_t p) const;

    double& at(size_t p);
    const double& at(size_t p) const;

    Vector& operator+=(const Vector&);
    Vector& operator-=(const Vector&);
    Vector& operator*=(double);
    Vector& operator/=(double);

    double magnitude() const;
private:
    std::array<double, N> x = {};
};





//NON-MEMBER TEMPLATE FUNCTION DEFINITIONS:

template <unsigned N> std::ostream& operator<<(std::ostream& os, const Vector<N>& rhs){
    os << "[";

    for(unsigned it = 0; it != N; ++it){
        os << approximate(rhs[it], 0);
        if(it != N-1) os << ", ";
    }

    os << "]";
    return os;
}

template <unsigned N> std::istream& operator>>(std::istream& is, Vector<N>& rhs){
    Vector<N> errorRet = rhs;

    for(unsigned it = 0; it != N; ++it)
        is >> rhs[it];

    if(!is)
        rhs = errorRet;

    return is;
}

template <unsigned N> bool operator==(const Vector<N>& lhs, const Vector<N>& rhs){
    return std::equal(lhs.x.begin(), lhs.x.end(), rhs.x.begin(), double_equals);
}

template <unsigned N> bool operator!=(const Vector<N>& lhs, const Vector<N>& rhs){
    return !(lhs == rhs);
}

template <unsigned N> Vector<N> operator+(const Vector<N>& lhs, const Vector<N>& rhs){
    Vector<N> sum = lhs;
    sum += rhs;
    return sum;
}

template <unsigned N> Vector<N> operator-(const Vector<N>& lhs, const Vector<N>& rhs){
    Vector<N> sum = lhs;
    sum -= rhs;
    return sum;
}

template <unsigned N> Vector<N> operator*(const Vector<N>& rhs, double d){
    Vector<N> product = rhs;
    product *= d;
    return product;
}

template <unsigned N> Vector<N> operator*(double d, const Vector<N>& rhs){
    return rhs*d;
}

template <unsigned N> Vector<N> operator/(const Vector<N>& rhs, double d){
    Vector<N> remain = rhs;
    remain /= d;
    return remain;
}

template <unsigned N> double dot_product(const Vector<N>& lhs, const Vector<N>& rhs) {
    return std::inner_product(lhs.x.begin(), lhs.x.end(), rhs.x.begin(), 0.0);
}

template <unsigned N> Vector<N> normalised(const Vector<N>& lhs) {
    Vector<N> ret = lhs;

    double length = ret.magnitude();
    if(!double_equals(length, 0.0)) ret /= length;
    return ret;
}

template <unsigned N> Vector<N> rotated_vector(const Vector<N>& lhs, size_t i, size_t j, double angle){
    Vector<N> ret = lhs;

    ret[i] = lhs[i]*cos(angle) - lhs[j]*sin(angle);
    ret[j] = lhs[j]*cos(angle) + lhs[i]*sin(angle);

    return ret;
}






//MEMBER FUNCTIONS:

template <unsigned N> template <typename... Args> inline Vector<N>::Vector(Args... li) {
    /*okay in the sense that it allows braced initialisation from only the right sized list.
    Not okay in the sense that is_constructible is true for any combination of arguments. See
    stackoverflow */
    static_assert(sizeof...(li) == N, "Wrong number of constructor arguments.");
    x = {static_cast<double>(li)...};
}

template <unsigned N> inline Vector<N>::operator bool() const {
    /*We expect more non-zero vectors than zero vectors. This algorithm should exit early for non-zero vectors.*/
    return !std::all_of(x.begin(), x.end(), [](double d){ return double_equals(d,0.0);});
}

template <unsigned N> inline double& Vector<N>::operator[](size_t p){
    return x[p];
}

template <unsigned N> inline const double& Vector<N>::operator[](size_t p) const {
    return x[p];
}

template <unsigned N> inline double& Vector<N>::at(size_t p){
    if(p >= N) throw std::out_of_range(std::string("Invalid coordinate specified for function ") + __func__ + ".");
    return x[p];
}

template <unsigned N> inline const double& Vector<N>::at(size_t p) const {
    if(p >= N) throw std::out_of_range(std::string("Invalid coordinate specified for function ") + __func__ + ".");
    return x[p];
}

template <unsigned N> inline Vector<N>& Vector<N>::operator+=(const Vector<N>& rhs){
    for(unsigned it = 0; it != N; ++it)
        x[it] += rhs.x[it];
    return *this;
}

template <unsigned N> inline Vector<N>& Vector<N>::operator-=(const Vector<N>& rhs){
    for(unsigned it = 0; it != N; ++it)
        x[it] -= rhs.x[it];
    return *this;
}

template <unsigned N> inline Vector<N>& Vector<N>::operator*=(double d){
    for(double& it : x)
        it *= d;
    return *this;
}

template <unsigned N> inline Vector<N>& Vector<N>::operator/=(double d){
    if(d == 0) throw std::domain_error(std::string("Division by zero in function ") + __func__ + ".");
    for(double& it : x)
        it /= d;
    return *this;
}

template <unsigned N> inline double Vector<N>::magnitude() const {
    return sqrt(std::accumulate(x.begin(), x.end(), 0.0, [](double sum, double d) {return sum + d*d;}));
}

#endif
