#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

int countFiles(const fs::path& dirPath, const std::string& extension) {
    int count = 0;
    try {
        if (fs::exists(dirPath) && fs::is_directory(dirPath)) {
            for (const auto& entry : fs::recursive_directory_iterator(dirPath)) {
                if (fs::is_regular_file(entry) && entry.path().extension() == extension) {
                    ++count;
                }
            }
        } else {
            std::cerr << "Error: Directory does not exist or is not a directory." << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
    return count;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <directory_path> <file_extension>" << std::endl;
        return 1;
    }

    std::string directoryPath = argv[1];
    std::string fileExtension = argv[2];

    int fileCount = countFiles(directoryPath, fileExtension);

    std::cout << "Number of files with extension " << fileExtension << " in directory "
              << directoryPath << ": " << fileCount << std::endl;

    return 0;
}
