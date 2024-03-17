#ifndef AAPERTURB_PERTURBRUN_H
#define AAPERTURB_PERTURBRUN_H

#include<string>
#include<filesystem>

namespace fs = std::filesystem;

void createdataset(const std::string, const std::string, const unsigned int, const unsigned int, double, double);
void perturbRun(fs::path, fs::path, unsigned int, double, double);



#endif //AAPERTURB_PERTURBRUN_H
