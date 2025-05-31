#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
using namespace std;

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

int main() {
        try {
            // string path1="chr1_tar.fa";
            string path1="test\\sekvenca_tar3.txt";

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
