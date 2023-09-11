//Header that contains the constants
#include<string>
#include<vector>
#include<map>

namespace amino_acids{
    namespace atoms {
        const std::vector<std::string> ALA = {"N", "CA", "C", "O", "CB"};
        const std::vector<std::string> ARG = {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"};
        const std::vector<std::string> ASN = {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"};
        const std::vector<std::string> ASP = {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"};
        const std::vector<std::string> CYS = {"N", "CA", "C", "O", "CB", "SG"};
        const std::vector<std::string> GLU = {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"};
        const std::vector<std::string> GLN = {"N", "CA", "C", "O",  "CB", "CG", "CD", "OE1", "NE2"};
        const std::vector<std::string> GLY = {"N", "CA", "C", "O",};
        const std::vector<std::string> HIS = {"N", "CA", "C", "O", "CB", "CG", "ND1", "CE1", "NE2", "CD2"};
        const std::vector<std::string> ILE = {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"};
        const std::vector<std::string> LEU = {"N", "CA", "C", "O",  "CB", "CG", "CD1", "CD2"};
        const std::vector<std::string> LYS = {"N", "CA", "C", "O",  "CB", "CG", "CD", "CE", "NZ"};
        const std::vector<std::string> MET = {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"};
        const std::vector<std::string> PHE = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CE1", "CZ", "CE2", "CD2"};
        const std::vector<std::string> PRO = {"N", "CA", "C", "O", "CB", "CG", "CD"};
        const std::vector<std::string> SER = {"N", "CA", "C", "O", "CB", "OG"};
        const std::vector<std::string> THR = {"N", "CA", "C", "O", "CB", "OG1", "CG2"};
        const std::vector<std::string> TRP = {"N", "CA", "C", "O", "CB", "CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"};
        const std::vector<std::string> TYR = {"N", "CA", "C", "O", "CB", "CD1", "CE1", "CZ", "CE2", "CD2", "OH"};
        const std::vector<std::string> VAL = {"N", "CA", "C", "O", "CB", "CG1", "CG2"};

        const std::map<std::string, std::vector<std::string>> AMINO_MAP =  {
                {"ALA", ALA}, {"ARG",  ARG}, {"ASN",  ASN}, {"ASP",  ASP}, {"CYS",  CYS},
                {"GLU", GLU}, {"GLN", GLN}, {"GLY", GLY}, {"HIS", HIS}, {"ILE", ILE},
                {"LEU", LEU}, {"LYS", LYS}, {"MET", MET}, {"PHE", PHE}, {"PRO", PRO},
                {"SER", SER}, {"THR", THR}, {"TRP", TRP}, {"TYR", TYR}, {"VAL", VAL}
        };
    };

    namespace axes{
        const std::vector<std::string> ALA = {"CB"};
        const std::vector<std::string> ARG = {"CB", "CG", "CD", "NE", "CZ"};
        const std::vector<std::string> ASN = {"CB", "CG"};
        const std::vector<std::string> ASP = {"CB", "CG"};
        const std::vector<std::string> CYS = {"CB", "SG"};
        const std::vector<std::string> GLU = {"CB", "CG", "CD"};
        const std::vector<std::string> GLN = {"CB", "CG", "CD"};
        const std::vector<std::string> GLY = {}; //no axes :(
        const std::vector<std::string> HIS = {"CB", "CG"};
        const std::vector<std::string> ILE = {"CB", "CG1"};
        const std::vector<std::string> LEU = {"CB", "CG"};
        const std::vector<std::string> LYS = {"CB", "CG", "CD", "CE", "NZ"};
        const std::vector<std::string> MET = {"CB", "CG", "SD", "CE"};
        const std::vector<std::string> PHE = {"CB", "CG"};
        const std::vector<std::string> PRO = {}; //no axes :(
        const std::vector<std::string> SER = {"CB", "OG"};
        const std::vector<std::string> THR = {"CB", "OG1", "CG2"};
        const std::vector<std::string> TRP = {"CB", "CG"};
        const std::vector<std::string> TYR = {"CB", "CD1", "CE1", "CZ"};
        const std::vector<std::string> VAL = {"CB"};

        const std::map<std::string, std::vector<std::string>> AMINO_MAP =  {
                {"ALA", ALA}, {"ARG",  ARG}, {"ASN",  ASN}, {"ASP",  ASP}, {"CYS",  CYS},
                {"GLU", GLU}, {"GLN", GLN}, {"GLY", GLY}, {"HIS", HIS}, {"ILE", ILE},
                {"LEU", LEU}, {"LYS", LYS}, {"MET", MET}, {"PHE", PHE}, {"PRO", PRO},
                {"SER", SER}, {"THR", THR}, {"TRP", TRP}, {"TYR", TYR}, {"VAL", VAL}
        };
    }

};