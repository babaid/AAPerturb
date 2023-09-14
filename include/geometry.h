//Header that deals with geometric operations on atoms
#ifndef AAPERTURB_GEOMETRY_H
#define AAPERTURB_GEOMETRY_H
#include<string>
#include<array>
#include<cmath>
#include "pdbparser.h"
#include "molecules.h"






void rotateCoordinatesAroundAxis(std::valarray<double>&,  const std::valarray<double>& axis, const double angleInDegrees);
std::valarray<double> findRotationAxis(const std::unique_ptr<Residue>&, const std::string&);


unsigned int detect_clashes(const std::vector<std::vector<double>>& matrix, double ythreshold);
std::valarray<double> calculateCentroid(const std::unique_ptr<Residue>& res);

std::vector<std::vector<double>> calculateDistanceMatrix(const std::map<char, std::vector<std::unique_ptr<Residue>>>& chainMap);
std::vector<std::vector<double>> calculateLocalDistanceMatrix(const std::map<char, std::vector<std::unique_ptr<Residue>>>& chainMap, const std::unique_ptr<Residue>& refres);
double calculateDistance(const Atom& atom1, const Atom& atom2);
double calculateRMSD(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2);
double norm(const std::valarray<double>&);
bool areResiduesNeighbors(const std::unique_ptr<Residue>& residue1, const std::unique_ptr<Residue>& residue2, double threshold);
const std::map<char, std::vector<int>> findInterfaceResidues(const std::map<char, std::vector<std::unique_ptr<Residue>>>& chainMap, double cutoff);

std::valarray<double> crossProduct(const std::valarray<double>& vector1, const std::valarray<double>& vector2);
#endif