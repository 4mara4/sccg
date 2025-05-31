#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;

struct Position {
    int start;
    int end;
};
ostream& operator<<(ostream& os, const Position& pos) {
    os << "[" << pos.start << ", " << pos.end << "]";
    return os;
}

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

    // Parsira preambulu iz final.txt
    void parsePreproc(const string& final) {
        ifstream in(final);
        if (!in) throw runtime_error("Ne mogu otvoriti " + final);

        string line;
        int stage = 0;
        while (getline(in, line)) {
            if (stage == 0) {
                meta = line + "\n";
                cout << "meta " << meta << endl;
                stage++;
            } else if (stage == 1) {
                lineLen = stoi(line);
                cout << "Line length  " << lineLen << endl;
                stage++;
            } else if (stage == 2 || stage == 3) {
                if (line.empty()) {
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


    // Čita referencu iz FASTA (preskače header)
    string readSeq(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Ne mogu otvoriti " + path);
    string header, line, seq;
    // preskoči header
    if (!getline(in, header)) 
        throw runtime_error("Prazna FASTA datoteka: " + path);

    // čitaj liniju po liniju
    while (getline(in, line)) {
        for (char c : line) {
            // uppercase
            char C = static_cast<char>(toupper(static_cast<unsigned char>(c)));
            // izbaci sve 'N'
            if (C != 'N') {
                seq.push_back(C);
            }
        }
    }
    return seq;
}

    // Rekonstruira sirovu sekvencu iz delta-podataka final.txt i reference
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
            // Preskači header, lineLen, L-list, N-list
            if (count <= 4) continue;
            if (line.empty()) continue;

            auto comma = line.find(',');
            if (comma != string::npos) {
                // delta,len
                int delta = stoi(line.substr(0, comma));
                int len   = stoi(line.substr(comma + 1));
                
                int b     = prevEnd + delta;
                int e     = b + len - 1;  // len = e-b
                /* cout << "  b " << b;
                cout << "  Len:   " << len;
                cout << "  e " << e << endl; */
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
                // literal mismatch
                out += line;
            }
        }
        return out;
    }

    // Ubacuje N-regije
   
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

    // Vraća lowercase u L-regijama
    void applyLower(const vector<Position>& L_list,
                    string& seq) {
        for (auto& p : L_list) {
            for (int i = p.start; i <= p.end && i < (int)seq.size(); ++i)
                seq[i] = tolower(seq[i]);
        }
    }

    // Ispisuje FASTA s linijskom dužinom
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

    void run() {
        //const string refFile = "chr20_ref.fa";
        const string refFile="chr18ref.fa";

        const string outFile = "output.fa";
        const string finalFile="final.txt";

        // 1) Parsiraj preambulu
        parsePreproc(finalFile);
        //cout << "l list " << L_list << endl;
        //cout << "n list " << N_list << endl;
        // 2) Učitaj referencu
        string refSeq = readSeq(refFile);

        // 3) Rekonstruiraj raw
        string raw = reconstruct(finalFile, refSeq);
        //cout << raw <<endl;

        refSeq.clear();
        refSeq.shrink_to_fit();

        // 4) Ubaci N-regije
        string withNs = N_list.empty() ? raw : insertNs(raw, N_list);
   
        raw.clear();
        raw.shrink_to_fit();

        // 5) Vraćanje lowercase
        applyLower(L_list, withNs);

        // 6) Ispiši FASTA
        writeFasta(outFile, meta, withNs, lineLen);
        cout << "Reconstructed FASTA: " << outFile << "\n";

        withNs.clear();
        withNs.shrink_to_fit();
    }
};

int main() {
    try {
        SCCGD decompressor;
        decompressor.run();
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
