#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;

/**
 * @brief Compares two text files line by line.
 *
 * Opens both files and reads them line by line simultaneously. If either file
 * cannot be opened, throws a runtime_error. As soon as a mismatch is found
 * (either differing lines or differing number of lines), prints the first
 * difference and returns false. Returns true if all lines match exactly.
 *
 * @param path1 Path to the first file.
 * @param path2 Path to the second file.
 * @return true if the files are identical, false otherwise.
 * @throws std::runtime_error if either file cannot be opened.
 * Marija writing
 */
bool files_are_equal(const std::string& path1, const std::string& path2) {
    std::ifstream f1(path1), f2(path2);
    if (!f1.is_open() || !f2.is_open()) {
        throw std::runtime_error("Ne mogu otvoriti jedan od fajlova.");
    }

    std::string line1, line2;
    int line_num = 1;
    while (true) {
        bool r1 = static_cast<bool>(std::getline(f1, line1));
        bool r2 = static_cast<bool>(std::getline(f2, line2));

        if (!r1 || !r2) {
            if (r1 != r2) {
                std::cout << "Fajlovi imaju različit broj linija. Prva razlika na liniji " << line_num << std::endl;
                return false;
            }
            return true;
        }

        if (line1 != line2) {
            std::cout << "Razlika na liniji " << line_num << ":\n";
            std::cout << "Fajl1: " << line1 << "\nFajl2: " << line2 << std::endl;
            return false;
        }

        ++line_num;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Upotreba: " << argv[0] << " <tar file>\n";
        return 1;
    }
        try {
            // string path1="chr1_tar.fa";
            string path1="test//" + string(argv[1]);

            string path2= "output.fa";
            bool same = files_are_equal(path1,path2);
            if (same) {
                std::cout << "Fileovi su isti.\n";
            } else {
                std::cout << "Fileovi se razlikuju.\n";
            }
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << "\n";
            return 1;
        }  
        return 0;
}
