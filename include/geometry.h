//Header that deals with geometric operations on atoms
#ifndef AAPERTURB_GEOMETRY_H
#define AAPERTURB_GEOMETRY_H
#include<string>
#include<array>
#include<cmath>
#include "molecules.h"






void rotateCoordinatesAroundAxis(std::valarray<double>&,  const std::valarray<double>& axis, const double angleInDegrees);
std::valarray<double> findRotationAxis(const Residue&, const std::string&);


unsigned int detect_clashes(const std::vector<std::vector<double>>& matrix, double ythreshold);
std::valarray<double> calculateCentroid(const Residue& res);

//double calculateDistance(const Atom& atom1, const Atom& atom2);
double calculateRMSD(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2);
double norm(const std::valarray<double>&);
bool areResiduesNeighbors(const Residue& residue1, const Residue& residue2, double threshold);
const std::map<char, std::vector<int>> findInterfaceResidues(const std::unique_ptr<std::map<char, std::vector<Residue>>> & chainMap, double cutoff);

std::valarray<double> crossProduct(const std::valarray<double>& vector1, const std::valarray<double>& vector2);
#endif