//Header that deals with geometric operations on atoms
#ifndef AAPERTURB_GEOMETRY_H
#define AAPERTURB_GEOMETRY_H
#include<string>
#include<array>
#include<cmath>
#include "molecules.h"


using Vector3 = std::array<double, 3>;

Vector3 operator-(const Vector3&, const Vector3&);
Vector3 operator+(const Vector3&, const Vector3&);

Vector3 pow(const Vector3, double);

void operator+=(Vector3&, const Vector3&);
void operator-=(Vector3&, const Vector3&);
bool operator==(const Vector3& a, const Vector3& b);

Vector3 operator/(const Vector3&, double);
double sum(const Vector3&);

void rotateCoordinatesAroundAxis(Vector3&, const Vector3 & , const Vector3 & axis, const double angleInDegrees);
Vector3 findRotationAxis(const Residue&, const std::string&);


//unsigned int detect_clashes(const std::vector<std::vector<double>>& matrix, double ythreshold);
Vector3 calculateCentroid(const Residue& res);

//double calculateDistance(const Atom& atom1, const Atom& atom2);
//double calculateRMSD(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2);
double norm(const Vector3&);
bool areResiduesNeighbors(const Residue& residue1, const Residue& residue2, double threshold);
//const std::map<char, std::vector<int>> findInterfaceResidues(const std::unique_ptr<std::map<char, std::vector<Residue>>> & chainmap, double cutoff);

Vector3 crossProduct(const Vector3& vector1, const Vector3& vector2);
#endif