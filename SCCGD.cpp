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

class SCCGD {
public:
    // Parsira preambulu iz interim.txt
    void parsePreproc(const string& interim,
                      string& meta,
                      int& lineLen,
                      vector<Position>& L_list,
                      vector<Position>& N_list) {
        ifstream in(interim);
        if (!in) throw runtime_error("Ne mogu otvoriti " + interim);
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
                if (line.empty()) {
                    stage++;
                    continue;
                }
                istringstream iss(line);
                int prev = 0, offset, length;
                auto &list = (stage == 2 ? L_list : N_list);
                while (iss >> offset >> length) {
                    list.push_back({prev + offset,
                                    prev + offset + length});
                    prev = list.back().end;
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
    string reconstruct(const string& interim,
                       const string& ref) {
        ifstream in(interim);
        if (!in) throw runtime_error("Ne mogu otvoriti " + interim);

        string line;
        int count = 0;
        int prevEnd = 0;
        string out;

        while (getline(in, line)) {
            ++count;
            // Preskači header, lineLen, L-list, N-list i prazan red (5 redova)
            if (count <= 5) continue;
            if (line.empty()) continue;

            auto comma = line.find(',');
            if (comma != string::npos) {
                // delta,len
                int delta = stoi(line.substr(0, comma));
                int len   = stoi(line.substr(comma + 1)) - delta;
                cout << "Len: " << len;
                int b     = prevEnd + delta;
                int e     = b + len;  // len = e-b
                cout << "b " << b;
                cout << "e " << e;
                if (b < 0 || e >= (int)ref.size()) {
                    throw runtime_error("Error: segment izvan ref granica: " +
                                        to_string(delta) + "+" +
                                        to_string(len) + " > " +
                                        to_string((int)ref.size()));
                }
                out += ref.substr(b, len + 1); // len+1 baza
                prevEnd = e;
            }
            else {
                // literal mismatch
                out += line;
                prevEnd = 0;
            }
        }
        return out;
    }

    // Ubacuje N-regije
    string insertNs(const string& raw,
                    const vector<Position>& N_list) {
        string out;
        int rawIdx = 0, pos = 0;
        for (auto& p : N_list) {
            while (pos < p.start && rawIdx < (int)raw.size()) {
                out.push_back(raw[rawIdx++]);
                pos++;
            }
            for (int i = p.start; i < p.end; ++i) {
                out.push_back('N');
                pos++;
            }
        }
        while (rawIdx < (int)raw.size())
            out.push_back(raw[rawIdx++]);
        return out;
    }

    // Vraća lowercase u L-regijama
    void applyLower(const vector<Position>& L_list,
                    string& seq) {
        for (auto& p : L_list) {
            for (int i = p.start; i < p.end && i < (int)seq.size(); ++i)
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
        out << "\n";
    }

    void run() {
        const string refFile = "sekvenca_ref.txt";
        const string interim = "interim.txt";
        const string outFile = "output.fa";

        string meta;
        int lineLen;
        vector<Position> L_list, N_list;

        // 1) Parsiraj preambulu
        parsePreproc(interim, meta, lineLen, L_list, N_list);

        // 2) Učitaj referencu
        string refSeq = readSeq(refFile);

        // 3) Rekonstruiraj raw
        string raw = reconstruct(interim, refSeq);

        // 4) Ubaci N-regije
        string withNs = N_list.empty() ? raw : insertNs(raw, N_list);

        // 5) Vraćanje lowercase
        applyLower(L_list, withNs);

        // 6) Ispiši FASTA
        writeFasta(outFile, meta, withNs, lineLen);
        cout << "Reconstructed FASTA: " << outFile << "\n";
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
