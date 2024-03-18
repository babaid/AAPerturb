#include "geometry.h"
#include "molecules.h"


#include<cmath>
#include<string>
#include<iostream>
#include<array>


using Vector3 = std::array<double, 3>;

// Function to rotate a vector around an arbitrary axis
void rotateCoordinatesAroundAxis(Vector3& vector, const Vector3& p, const Vector3& axis, double angle) {
    double cs = std::cos(angle);
    double si = std::sin(angle);
    double t = 1 - cs;

    double x = vector[0];
    double y = vector[1];
    double z = vector[2];

    double a = p[0];
    double b = p[1];
    double c = p[2];

    double u = axis[0];
    double v = axis[1];
    double w = axis[2];

    double x_new = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z ))*t + x*cs + (- c*v + b*w - w*y + v*z)*si;
    double y_new = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z ))*t + y*cs + (c*u - a*w + w*x - u*z)*si;
    double z_new = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z ))*t + z*cs + (- b*u + a*v - v*x + u*y)*si;

    vector[0] = x_new;
    vector[1] = y_new;
    vector[2] = z_new;

    //bool memoryleak = false; // Anil's contribution.
}


Vector3 findRotationAxis(const Residue& residue, const std::string& axis)
{
    unsigned cntr{0};
    for (const Atom& atom : residue.atoms) {
        if (atom.name == axis) {
            return residue.atoms.at(cntr).coords;
        }
        cntr++;
    }
    return  {};
}


unsigned int detect_clashes(const std::vector<std::vector<double>>& matrix, double threshold){
    unsigned int clash_cnt = 0;
    for (const auto& row : matrix) {
        for (double value : row) {
            if (value < threshold) {
                clash_cnt++;
            }
        }
    }
    return clash_cnt;
}


double calculateRMSD(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
    if (atoms1.size() != atoms2.size()) {
        std::cerr << "Error: Sets of atomic coordinates must have the same size." << std::endl;
        return -1.0; // Return an error value
    }

    double sumSquaredDifferences = 0.0;
    std::size_t numAtoms = atoms1.size();


    for (std::size_t i = 0; i < numAtoms; ++i) {
        sumSquaredDifferences +=  sum(pow(atoms1[i].coords - atoms2[i].coords, 2));
    }

    return std::sqrt(sumSquaredDifferences / static_cast<double>(numAtoms));
}

Vector3 calculateCentroid(const Residue& res)
{
    Vector3 centroid{0., 0., 0.};
    auto num_atoms = (double)res.atoms.size();
    for (const Atom& atom: res.atoms) {
        centroid+=atom.coords;
    }
    return centroid/num_atoms;
}



bool areResiduesNeighbors(const Residue& residue1, const Residue& residue2, double threshold) {
    Vector3 cog1 = calculateCentroid(residue1);
    Vector3 cog2 = calculateCentroid(residue2);
    double distance = std::sqrt(sum(pow(cog2-cog1, 2)));
    if (distance<threshold)
    {
        return true;
    }
    return false;
}

Vector3 operator-(const Vector3& a, const Vector3& b) {
    Vector3 v;
    for(unsigned i=0; i<3; i++)
    {
        v[i]=a[i]-b[i];
    }
    return v;
}

Vector3 operator+(const Vector3& a, const Vector3& b) {
    Vector3 v{a};
    for(unsigned i=0; i<3; i++)
    {
        v[i]+=b[i];
    }
    return v;
}

void operator+=(Vector3& a, const Vector3& b)
{
    for (unsigned i{0}; i<3; i++) a[i]+=b[i];
}

void operator-=(Vector3& a, const Vector3& b)
{
    for (unsigned i{0}; i<3; i++) a[i]-=b[i];
}

Vector3 operator/(const Vector3& a, double d)
{
    Vector3 v;
    for (unsigned i{0}; i<3; i++) v[i]=a[i]/d;
    return v;
}
Vector3 pow(const Vector3 a, double p)
{
    Vector3 v;
    for(unsigned i{0}; i<3; i++)
    {
        v[i] = std::pow(a[i], p);
    }
    return v;
}

double sum(const Vector3& v)
{
    double out{0};
    for(auto& el:v)
    {
        out+=el;
    }
    return out;
}

//maybe useful later.
Vector3 crossProduct(const Vector3& vector1, const Vector3& vector2) {
    if (vector1.size() != 3 || vector2.size() != 3) {
        std::cerr << "Error: Input vectors must have size 3." << std::endl;
        return Vector3({0., 0., 0.}); // Return a zero vector
    }

    double resultX = vector1[1] * vector2[2] - vector1[2] * vector2[1];
    double resultY = vector1[2] * vector2[0] - vector1[0] * vector2[2];
    double resultZ = vector1[0] * vector2[1] - vector1[1] * vector2[0];

    return Vector3({resultX, resultY, resultZ});
}

double norm(const Vector3& vec)
{
    return std::sqrt(sum(pow(vec, 2)));
}