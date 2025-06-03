#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;

/**
 * @brief Returns CPU thread time in nanoseconds.
 * Dora writing
 */
long long getCPUTime() {
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0) {
        return static_cast<long long>(ts.tv_sec) * 1'000'000'000LL + ts.tv_nsec; // nanosekunde
    } else {
        return 0LL;
    }
}
struct Position {
    int start;
    int end;
};

/**
 * @brief Outputs a Position as “[start, end]”.
 * Dora writing
 */
ostream& operator<<(ostream& os, const Position& pos) {
    os << "[" << pos.start << ", " << pos.end << "]";
    return os;
}

/**
 * @brief Outputs a vector of Position objects in braces.
 * Dora writing
 */
ostream& operator<<(ostream& os, const vector<Position>& v) {
    os << "{ ";
    for (const auto& p : v) {
        os << p << " ";
    }
    os << "}";
    return os;
}


class SCCGD {
public:
    string meta;
    int lineLen;
    vector<Position> L_list, N_list;

    /**
     * @brief Parses the first four lines from final.txt (meta, lineLen, L_list, N_list).
     * Dora + Marija writing
     */
    void parsePreproc(const string& final) {
        ifstream in(final);
        if (!in) throw runtime_error("Ne mogu otvoriti " + final);

        string line;
        int stage = 0;
        while (getline(in, line)) {
            if (stage == 0) {
                meta = line + "\n";
                stage++;
            } else if (stage == 1) {
                lineLen = stoi(line);
                stage++;
            } else if (stage == 2 || stage == 3) {
                if (line == "x") {
                    stage++;
                    continue;
                }

                istringstream iss(line);
                string token;
                int prev = 0;
                auto& list = (stage == 2 ? L_list : N_list);

                while (iss >> token) {
                    size_t colon = token.find(':');
                    if (colon == string::npos) continue;

                    int offset = stoi(token.substr(0, colon));
                    int length = stoi(token.substr(colon + 1));
                    int start = prev + offset;
                    int end = start + length - 1;

                    list.push_back({start, end});
                    prev = end + 1;
                }

                stage++;
            } else {
                break;
            }
        }
    }

    /**
     * @brief Reads a FASTA sequence (skips header) and removes all 'N' characters.
     * Marija writing
     */
    string readSeq(const string& path) {
        ifstream in(path);
        if (!in) throw runtime_error("Ne mogu otvoriti " + path);
        string header, line, seq;
        if (!getline(in, header)) throw runtime_error("Prazna FASTA datoteka: " + path);

        while (getline(in, line)) {
            for (char c : line) {
                char C = static_cast<char>(toupper(static_cast<unsigned char>(c)));
                if (C != 'N') {
                    seq.push_back(C);
                }
            }
        }
        return seq;
    }

    /**
     * @brief Reconstructs the raw sequence using delta entries and the reference string.
     * Dora + Marija writing
     */
    string reconstruct(const string& final,
                       const string& ref) {
        ifstream in(final);
        if (!in) throw runtime_error("Ne mogu otvoriti " + final);

        string line;
        int count = 0;
        int prevEnd = 0;
        string out;

        while (getline(in, line)) {
            ++count;
            if (count <= 4) continue;
            if (line.empty()) continue;

            auto comma = line.find(',');
            if (comma != string::npos) {
                int delta = stoi(line.substr(0, comma));
                int len   = stoi(line.substr(comma + 1));
                
                int b     = prevEnd + delta;
                int e     = b + len - 1;  
                if (b < 0 || e > (int)ref.size()) {
                    throw runtime_error("  Error: segment izvan ref granica: " +
                                        to_string(delta) + "+" +
                                        to_string(len) + " > " +
                                        to_string((int)ref.size()));
                }
                out += ref.substr(b, len); 
                prevEnd = e;
            }
            else {
                out += line;
            }
        }
        return out;
    }
    /**
     * @brief Inserts 'N' regions at the specified positions.
     * Marija writing
     */   
    string insertNs(const string& raw, const vector<Position>& N_list) {
        string out;
        int rawIdx = 0;
        int pos = 0;
        size_t n = raw.size();

        for (const auto& p : N_list) {
            while (pos < p.start && rawIdx < (int)n) {
                out.push_back(raw[rawIdx++]);
                pos++;
            }

            for (int i = p.start; i <= p.end; ++i) {
                out.push_back('N');
                pos++;
            }
        }

        while (rawIdx < (int)n) {
            out.push_back(raw[rawIdx++]);
            pos++;
        }

        return out;
    }

    /**
     * @brief Applies lowercase characters within intervals from L_list.
     * Marija writing
     */
    void applyLower(const vector<Position>& L_list,
                    string& seq) {
        for (auto& p : L_list) {
            for (int i = p.start; i <= p.end && i < (int)seq.size(); ++i)
                seq[i] = tolower(seq[i]);
        }
    }

    /**
     * @brief Writes a FASTA file, wrapping lines at lineLen.
     * Dora writing
     */
    void writeFasta(const string& path,
                    const string& header,
                    const string& seq,
                    int lineLen) {
        ofstream out(path);
        if (!out) throw runtime_error("Ne mogu pisati " + path);
        out << header;
        for (int i = 0; i < (int)seq.size(); ++i) {
            if (i > 0 && i % lineLen == 0) out << "\n";
            out << seq[i];
        }
    }

    /**
     * @brief Main entry: reconstructs and writes the final FASTA.
     * Dora + Marija writing
     */
    void run(string refFile) {

        const string outFile = "output.fa";
        const string finalFile="final.txt";

        parsePreproc(finalFile);
        string refSeq = readSeq(refFile);
        string raw = reconstruct(finalFile, refSeq);

        refSeq.clear();
        refSeq.shrink_to_fit();
        string withNs = N_list.empty() ? raw : insertNs(raw, N_list);
   
        raw.clear();
        raw.shrink_to_fit();
        applyLower(L_list, withNs);

        writeFasta(outFile, meta, withNs, lineLen);
        cout << "Reconstructed FASTA: " << outFile << "\n";

        withNs.clear();
        withNs.shrink_to_fit();
    }
};

int main(int argc, char* argv[]) {
    linux long long start_time = getCPUTime();
    if (argc != 2) {
        std::cerr << "Upotreba: " << argv[0] << " <ref file name>\n";
        return 1;
        }
    
    try {
        const string refFile="test/" + string(argv[1]);
        SCCGD decompressor;
        decompressor.run(refFile);
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }
    long long end_time = getCPUTime();  
    
    cout << "Execution time " << (end_time - start_time) << endl; 
    return 0;
}
