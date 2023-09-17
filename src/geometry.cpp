#include<cmath>
#include<string>
#include<iostream>
#include<valarray>
#include<array>
#include "constants.h"
#include "geometry.h"
#include "molecules.h"


void rotateCoordinatesAroundAxis(std::valarray<double>& vector, const std::valarray<double>& axis, const double angleInDegrees)
{
    //std::valarray<double> axis = axis / std::sqrt(std::pow(axis[0], 2) + std::pow(axis[1], 2) + std::pow(axis[2], 2));
    if(vector.size() != 3 && axis.size() != 3 ) throw std::out_of_range("Coordinates should be 3D.");
    
    double angle = angleInDegrees * M_PI / 180.0;
    double cosAngle = std::cos(angle);
    double sinAngle = std::sin(angle);

    // Calculate the rotation matrix
    

    std::array<std::array<double, 3>, 3> rotationMatrix{{   {cosAngle + (1 - cosAngle) * std::pow(axis[0], 2), (1 - cosAngle) * axis[0] * axis[1] - sinAngle * axis[2],(1 - cosAngle) * axis[0] * axis[2] + sinAngle * axis[1]},
                                                            {(1 - cosAngle) * axis[1] * axis[0] + sinAngle * axis[2], cosAngle + (1 - cosAngle) * std::pow(axis[1], 2), (1 - cosAngle) * axis[1] * axis[2] - sinAngle * axis[0]},
                                                            {(1 - cosAngle) * axis[2] * axis[0] - sinAngle * axis[1],(1 - cosAngle) * axis[2] * axis[1] + sinAngle * axis[0], cosAngle + (1 - cosAngle) * std::pow(axis[2], 2)} }};


    // Apply the rotation matrix to the input vector
    double oldX = vector[0];
    double oldY = vector[1];
    double oldZ = vector[2];

    vector[0] = rotationMatrix[0][0] * oldX + rotationMatrix[0][1] * oldY + rotationMatrix[0][2] * oldZ;
    vector[1] = rotationMatrix[1][0] * oldX + rotationMatrix[1][1] * oldY + rotationMatrix[1][2] * oldZ;
    vector[2] = rotationMatrix[2][0] * oldX + rotationMatrix[2][1] * oldY + rotationMatrix[2][2] * oldZ;
}

std::valarray<double> findRotationAxis(const Residue& residue, const std::string& axis)
{
    for (const Atom& atom : residue.atoms) {
        if (atom.name == axis) {
            return atom.coords;
        }
    }
    return  {};
}

double calculateDistance(const Atom& atom1, const Atom& atom2) {
    return std::sqrt(std::pow(atom2.coords - atom1.coords, 2).sum());
}

std::vector<std::vector<double>> calculateDistanceMatrix(const std::unique_ptr<std::map<char, std::vector<Residue>>> & chainMap) {
    int numAtoms = 0;

    // Calculate the total number of atoms in the chainMap
    for (const auto& chainEntry : *chainMap) {
        for (auto & residue : chainEntry.second) {
            numAtoms += residue.atoms.size();
        }
    }

    std::vector<std::vector<double>> distanceMatrix(numAtoms, std::vector<double>(numAtoms, 0.0));

    int currentIndex = 0;

    // Iterate through the chainMap structure
    for (const auto& chainEntry1 : *chainMap) {
        for (auto const& residue1 : chainEntry1.second) {
            for (const Atom& atom1 : residue1.atoms) {
                int otherIndex = 0;
                for (const auto& chainEntry2 : *chainMap) {
                    for (auto const&  residue2 : chainEntry2.second) {
                        for (const Atom& atom2 : residue2.atoms) {
                            double distance = calculateDistance(atom1, atom2);
                            distanceMatrix[currentIndex][otherIndex] = distance;
                            if (currentIndex == otherIndex) {
                                distanceMatrix.at(currentIndex).at(otherIndex) = std::numeric_limits<double>::max();
                            }
                            otherIndex++;
                        }
                    }
                }
                currentIndex++;
            }
        }
    }
    return distanceMatrix;
}

std::vector<std::vector<double>> calculateLocalDistanceMatrix(const std::unique_ptr<std::map<char, std::vector<Residue>>> & chainMap, const Residue& refres)
{

    //This needs to be moved to the class Protein as it can be stored at loading of a structure.
    int numAtoms = 0;

    // Calculate the total number of atoms in the chainMap
    for (const auto& chainEntry : *chainMap) {
        for (auto const& residue : chainEntry.second) {
            numAtoms += residue.atoms.size();
        }
    }

    std::vector<std::vector<double>> distanceMatrix =  std::vector<std::vector<double>>(numAtoms, std::vector<double>(refres.atoms.size(), 0.0));

    int currentIndex = 0;

    // Iterate through the chainMap structure
    for (const auto& chainEntry1 : *chainMap) {
        for (auto const& residue1: chainEntry1.second) {
            for (const Atom& atom1: residue1.atoms) {
                int otherIndex = 0;
                for (const Atom& atom2: refres.atoms) {
                    double distance = calculateDistance(atom1, atom2);
                    distanceMatrix.at(currentIndex).at(otherIndex) = distance;
                    // Set diagonal elements to max
                    if (residue1.resSeq == refres.resSeq && atom1 == atom2) {
                        distanceMatrix.at(currentIndex).at(otherIndex) = std::numeric_limits<double>::max();
                    }
                    otherIndex++;
                }
                currentIndex++;
            }
        }
    }
    return distanceMatrix;
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
        sumSquaredDifferences +=  std::pow(atoms1[i].coords - atoms2[i].coords, 2).sum();
    }

    return std::sqrt(sumSquaredDifferences / static_cast<double>(numAtoms));
}

std::valarray<double> calculateCentroid(const Residue& res)
{
    std::valarray<double> centroid{0., 0., 0.};
    double num_atoms = res.atoms.size();
    for (const Atom& atom: res.atoms) {
        centroid+=atom.coords;
    }
    return centroid/num_atoms;
}



bool areResiduesNeighbors(const Residue& residue1, const Residue& residue2, double threshold) {
    std::valarray<double> cog1 = calculateCentroid(residue1);
    std::valarray<double> cog2 = calculateCentroid(residue2);
    double distance = std::sqrt(std::pow(cog2-cog1, 2).sum());
    if (distance<threshold)
    {
        return true;
    }
    return false;
}

const std::map<char, std::vector<int>> findInterfaceResidues(const std::unique_ptr<std::map<char, std::vector<Residue>>> & chainMap, double cutoff) {
    std::map<char, std::vector<int>>  interfaceResidues;

    for (const auto& chainEntry1 : *chainMap) {
        interfaceResidues[chainEntry1.first] = std::vector<int>();
        for (auto const& residue1 : chainEntry1.second) {
            for (const auto& chainEntry2 : *chainMap) {
                if (chainEntry1.first != chainEntry2.first) {
                    for (auto const & residue2 : chainEntry2.second) {
                        if (areResiduesNeighbors(residue1, residue2, cutoff)) {
                            // Check if residues are adjacent by comparing residue sequence numbers
                                interfaceResidues[chainEntry1.first].emplace_back(residue1.resSeq-1);
                                //interfaceResidues.push_back(residue1.resSeq);
                                break;  // No need to check further for this residue1
                        }
                    }
                }
            }
        }
    }

    return interfaceResidues;
}


std::valarray<double> crossProduct(const std::valarray<double>& vector1, const std::valarray<double>& vector2) {
    if (vector1.size() != 3 || vector2.size() != 3) {
        std::cerr << "Error: Input vectors must have size 3." << std::endl;
        return std::valarray<double>(0.0, 3); // Return a zero vector
    }

    double resultX = vector1[1] * vector2[2] - vector1[2] * vector2[1];
    double resultY = vector1[2] * vector2[0] - vector1[0] * vector2[2];
    double resultZ = vector1[0] * vector2[1] - vector1[1] * vector2[0];

    return std::valarray<double>({resultX, resultY, resultZ});
}

double norm(const std::valarray<double>& vec)
{
    return std::sqrt(std::pow(vec, 2).sum());
}