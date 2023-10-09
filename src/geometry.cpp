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