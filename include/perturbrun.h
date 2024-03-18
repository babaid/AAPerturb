#ifndef AAPERTURB_PERTURBRUN_H
#define AAPERTURB_PERTURBRUN_H

#include<string>
#include<filesystem>

namespace fs = std::filesystem;

void createdataset(const std::string&, const std::string&, size_t, size_t, double, double);
void perturbRun(fs::path&, fs::path&, size_t, double, double);



#endif //AAPERTURB_PERTURBRUN_H
