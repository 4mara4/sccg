#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <windows.h>  // Just testing on windows
using namespace std;
/*   ON LINUX
long long getCPUTime() {
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0) {
        return static_cast<long long>(ts.tv_sec) * 1'000'000'000LL + ts.tv_nsec; // nanosekunde
    } else {
        return 0LL;
    }
} */


struct Position {
    int startInRef;
    int endInRef;
    int startInTar;
    int endInTar;

};

/**
 * @class K_mer
 * @brief Represents a fixed-length substring (“k-mer”) of a genomic sequence and its position.
 *
 * The K_mer class encapsulates a k-mer string and the index at which this k-mer
 * begins in the full sequence. It provides simple getters and setters for both
 * the k-mer content and its start position, allowing building and querying
 * of a local hash map for fast k-mer lookups during reference-based compression.
 *
 * Dora writing
 */

class K_mer {
private:
    string k_mer;
    int start;

public:
    K_mer() {
        k_mer = "";
        start = 0;
    }
    K_mer(const string& kmer, int start_) {
        k_mer = kmer;
        start = start_;
    }

    string getKmer() const {
        return k_mer;
    }
    void setKmer(const string& kmer) {
        k_mer = kmer;
    }

    int getStart() const {
        return start;
    }
    void setStart(int start_) {
        start = start_;
    }
};

/**
 * @class SCCGC
 * @brief Core class for reference-based genome compression: handles sequence I/O and k-mer indexing.
 *
 * The SCCGC class encapsulates the main functionality required to prepare genomic data
 * for local matching against a reference:
 *   - Reading and parsing FASTA-formatted sequences for both local and global matching phases.
 *   - Tracking metadata (e.g., chromosome label, sequence length).
 *   - Building a hash-based index of k-mers to accelerate substring lookups.
 *
 * Individual matching, formatting, and compression steps are implemented in
 * separate member functions, each documented at its own definition.
 *
 * Dora writing
 */

class SCCGC {
public:
    string meta_data;
    int length;
    unordered_map<size_t, vector<K_mer>> localHash;
    unordered_map<size_t, vector<K_mer>> globalHash;
    vector<unordered_map<size_t, vector<K_mer>>> refLocalHashes;
    vector<unordered_map<size_t, vector<K_mer>>> tarLocalHashes;


    /**
     * @brief Reads a single-sequence FASTA file for local k-mer matching.
     *
     * This function opens the specified file, expects the first line to be a metadata/header,
     * and concatenates all subsequent lines into one continuous sequence string.
     * It stores the header in the class’s `meta_data` member and records the length
     * of the first data line in `length`. Designed for preparing a sequence for local matching.
     *
     * @param filename Path to the input FASTA file.
     * @return A std::string containing the full sequence (with all line breaks removed).
     * @throws std::runtime_error if the file cannot be opened or is empty.
     *
     * Dora writing
     */

    string LocReadSeq(const string& filename) { 
        ifstream sequence_file(filename);
        if (!sequence_file.is_open()) {
            throw runtime_error("Could not open file: " + filename);
        }
        stringstream ss;
        string line;
        bool is_first_line = true;
        if(getline(sequence_file, line)) {
            meta_data = line;
        }else {
            throw runtime_error("File is empty: " + filename);
        }
        while (getline(sequence_file, line)) {
            ss << line;
            if(is_first_line) {
                length=line.length();
                is_first_line = false;
            } 
        }
        sequence_file.close();
        return ss.str();
    }

    /**
     * @brief Reads a FASTA reference file for global k-mer matching.
     *
     * Opens the specified FASTA file, skips the header line, and processes each subsequent line
     * by converting all characters to uppercase and removing any 'N' bases. The cleaned sequence
     * is concatenated into a single string and returned.
     *
     * @param filename Path to the input FASTA reference file.
     * @return A std::string containing the full reference sequence with no 'N' characters.
     * @throws std::runtime_error if the file cannot be opened or if it is empty.
     *
     * Dora writing
     */
    string GloReadRefSeq(const string& filename) {
    
        ifstream sequence_file(filename);
        if (!sequence_file.is_open()) {
            throw runtime_error("Could not open file: " + filename);
        }
        stringstream ss;
        string line;
        if(!getline(sequence_file, line)) {
            throw runtime_error("File is empty: " + filename);
        }

        while (getline(sequence_file, line)) {
            for (char c : line) {
                if(!isupper(static_cast<unsigned char>(c))) {
                    c = static_cast<char>(toupper(static_cast<unsigned char>(c)));
                }
                if (c != 'N') {
                ss << c;
                }
            }
        }
        sequence_file.close();
        return ss.str();
    }

    /**
     * @brief Reads a FASTA target file for global k-mer matching and records lowercase/N regions.
     *
     * Opens the given input FASTA file, writes its header and line-length metadata to the specified output file,
     * then scans each sequence line to identify regions of lowercase bases and runs of 'N's. It logs the
     * start offsets and lengths of those regions to the output, converting all bases to uppercase
     * (excluding 'N's from the returned sequence). Finally, it returns a cleaned sequence string
     * with all 'N' characters removed.
     *
     * @param inputfilename  Path to the input FASTA target file.
     * @param outputfilename Path to the file where header, line-length, lowercase/N region data are appended.
     * @return A std::string containing the target sequence with all 'N' characters stripped out.
     * @throws std::runtime_error if either the input or output file cannot be opened.
     *
     * Dora writing
     */

    string GloReadTarSeq(const string& inputfilename, const string& outputfilename){
        ifstream infile(inputfilename);
        ofstream outfile(outputfilename, ios::app); 
        stringstream Llist, Nlist, clean_seq;

        if (!infile.is_open()) {
            throw runtime_error("Could not open file: " + inputfilename);
        }
        if (!outfile.is_open()) {
            throw runtime_error("Could not open file: " + outputfilename);
        }

        string meta_data, line;
        getline(infile, meta_data);
        outfile << meta_data << '\n';
        
        int line_length = 0, length_read = 0, length_lower = 0, lenght_N = 0, end_L = 0, end_N = 0, start_L = 0, start_N = 0;
        bool lowercase_mode = false, N_mode = false, first_line = true;

        while (getline(infile, line)) {
            line_length = line.length();

            if (first_line) {
                outfile << line_length << '\n';
                first_line = false;
            }

            for (int i = 0; i < line_length; ++i) {
                char current_ch = line[i];

                if (islower(current_ch)) {
                    length_lower++;
                    if (!lowercase_mode) {
                        lowercase_mode = true;
                        start_L = i + length_read;
                        int offset=start_L - end_L;
                        Llist << (offset) << ':';
                    }
                    current_ch = toupper(current_ch);
                } else {
                    if (lowercase_mode) {
                        Llist << length_lower << ' '; 
                        end_L = start_L + length_lower - 1;
                    }
                    lowercase_mode = false;
                    length_lower = 0;
                }

                if (current_ch == 'N') { 
                    lenght_N++;
                    if (!N_mode)  {
                        N_mode = true;
                        start_N = i + length_read;
                        int offset = start_N - end_N;
                        Nlist <<offset << ':';
                    }
                } else {
                    if (N_mode) {
                        Nlist << lenght_N << ' ';
                        end_N = start_N + lenght_N -1;
                    }
                    N_mode = false;
                    lenght_N = 0;

                    clean_seq << current_ch;
                }
            }
            length_read += line_length;
        }

        if(lowercase_mode) {
            Llist << length_lower << ' ';
        }
        if(N_mode) {
            Nlist << lenght_N << ' ';
        }
        
        outfile << Llist.str() << '\n';
        outfile << Nlist.str() << '\n';

        infile.close();
        outfile.close();

        return clean_seq.str();
    }

    /**
     * @brief Builds a local hash index of k-mers for fast substring lookup.
     *
     * Iterates over the given sequence, extracts each k-length substring (k-mer),
     * and inserts it into an unordered_map where the key is the hash of the k-mer
     * and the value is a vector of K_mer instances containing the k-mer string
     * and its start position. Consecutive 'N'-only k-mers are skipped to avoid
     * indexing non-informative regions.
     *
     * @param seq The input sequence string to index.
     * @param kmer_length The length of each k-mer to extract and index.
     *
     * Dora + Marija writing
     */
    
    void createLocalHash(const string& seq, int kmer_length) {
        localHash.clear();
        int start_ = 0;
        while (start_ <= seq.length()- kmer_length) {
            string kmer = seq.substr(start_, kmer_length);
            if(kmer == string(kmer_length, 'N')) {                             //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! reminder to check later
                    start_= start_ + kmer_length ; // skip the N sequence   
                    continue; // skip this k-mer    
            }
            K_mer kmer_instance(kmer, start_); // reading the k-mer and its starting position
            size_t key = hash<string>{}(kmer);
            localHash[key].push_back(kmer_instance);
            start_++;         
        }
    }

    /**
     * @brief Builds a global hash index of k-mers, storing full K_mer objects.
     *
     * Iterates over the given sequence, extracts each k-length substring (k-mer),
     * skips any k-mers containing 'N', and inserts a K_mer (string+start) into
     * an unordered_map keyed by the k-mer’s hash. This lets you look up not only
     * the positions but also carry around the exact k-mer text if needed.
     *
     * @param seq           The input reference sequence string to index.
     * @param kmer_length   The length of each k-mer to extract and index.
     *
     * Marija writing
     */

    void createGlobalHash(const string& seq, int kmer_length) {
        globalHash.clear();
        int limit = seq.length() - kmer_length;
        for(int i = 0; i <= limit; ++i) {
            string kmer = seq.substr(i, kmer_length);
            if(kmer.find('N') != string::npos) {
                continue;
            }
            size_t key = hash<string>{}(kmer);
            K_mer entry(kmer, i);
            globalHash[key].push_back(entry);
        }
    }

    /**
    * @brief Splits a sequence into consecutive blocks of given size.
    *
    * @param seq         The full sequence string to split.
    * @param block_size  Maximum length of each block.
    * @return vector<string>  List of sequence blocks.
    *
    * Marija writing
    */
    vector<string> createBlocks(const string& seq, int block_size) {
        vector<string> blocks;
        int total = seq.size();
        for (int start = 0; start < total; start += block_size) {
            int len = min(block_size, total - start);
            blocks.push_back(seq.substr(start, len));
        }
        return blocks;
    }

    vector<Position> localMatch(const string& reference,
                            const string& target,
                            int kmer_length) {
        vector<Position> matches;
        int index = 0, tlen = target.length();
        while (index + kmer_length <= tlen) {
            string kmer = target.substr(index, kmer_length);
            size_t key = hash<string>{}(kmer);
            if (!localHash.count(key)) { index++; continue; }
            auto &klist = localHash[key];
            int bestStart = INT_MAX, bestExtend = 0;
            for (auto &km : klist) {
                if (km.getKmer() != kmer) continue;
                int incre = 0;
                int rpos = km.getStart() + kmer_length;
                int tpos = index + kmer_length;
                while (rpos < (int)reference.length() && tpos < tlen && reference[rpos] == target[tpos]) {
                    incre++; rpos++; tpos++;
                }
                if (incre > bestExtend || (incre == bestExtend && km.getStart() < bestStart)) {
                    bestExtend = incre;
                    bestStart = km.getStart();
                }
            }
            if (bestStart == INT_MAX) { index++; continue; }
            Position p;
            p.startInTar = index;
            p.endInTar   = index + kmer_length + bestExtend - 1;
            p.startInRef = bestStart;
            p.endInRef   = bestStart + kmer_length + bestExtend - 1;
            matches.push_back(p);
            index += kmer_length + bestExtend + 1;
        }
        return matches;
    }

     vector<Position> globalMatch(const string& ref, const string& tar, int klen) {
        createGlobalHash(ref, klen);
        vector<Position> out;
        int idx=0;
        while (idx + klen <= tar.size()) {
            string k=tar.substr(idx,klen);
            size_t h=hash<string>{}(k);
            if (!globalHash.count(h)) { idx++; continue; }
            int bestS=INT_MAX, bestE=0;
            for (auto &km: globalHash[h]) if (km.getKmer()==k) {
                int ext=0, r=km.getStart()+klen, t=idx+klen;
                while (r<ref.size()&&t<tar.size()&&ref[r]==tar[t]){ext++;r++;t++;}
                if (ext>bestE||(ext==bestE&&km.getStart()<bestS)) {bestE=ext;bestS=km.getStart();}
            }
            if (bestS!=INT_MAX) {
                out.push_back({bestS, bestS+klen+bestE-1, idx, idx+klen+bestE-1});
                idx += klen+bestE;   // check later if this is correct
            } else idx++;
        }
        return out;
    }

    /*
 * format_matches: Formats and outputs aligned and unaligned segments between
 * a reference and target sequence based on a list of alignment positions.
 * 
 * Arguments:
 *  - list: vector of Position structs defining aligned regions
 *  - target: the target sequence string
 *  - sor: optional offset to shift reference positions (default 0)
 *  - fileName: optional output file to save the formatted result
 * 
 * Returns:
 *  - the last Position from the input list
 * 
 * Functionality:
 *  - Writes aligned regions as reference position ranges (adjusted by sor)
 *  - Inserts unaligned/mismatched sequence fragments from target as-is
 *  - Optionally saves the result to a file
 * Dora writing
 */
Position format_matches(const vector<Position>& list, const string& target,
                        int sor = 0, const string& fileName = "") {
    string result;
    int trouble_parts = 0, endRef = 0, endInTar = 0;

    for (size_t i = 0; i < list.size(); ++i) {
        int startInTar = list[i].startInTar;
        int startInRef = list[i].startInRef;
        int endInRef = list[i].endInRef;
        int endInTarCurr = list[i].endInTar;

        if (i == 0) {
            endInTar = endInTarCurr;
            endRef = max(endRef, endInRef);

            if (startInTar > 0) {
                string pre_allign = target.substr(0, startInTar);
                result += pre_allign + "\n";    // we write the pre-aligned part
                trouble_parts += pre_allign.length();
            }
            result += to_string(startInRef + sor) + "," + to_string(endInRef + sor) + "\n";
            continue;
        }

        int endInTarPrev = list[i - 1].endInTar;
        string mismatch = target.substr(endInTarPrev + 1, startInTar - endInTarPrev - 1);   // mismatch part between two alignments	

        if (!mismatch.empty()) {
            result += mismatch + "\n";     // if there was a mismatched part we need to write it down
            trouble_parts += mismatch.length();
        }

        result += to_string(startInRef + sor) + "," + to_string(endInRef + sor) + "\n";    // we write the aligned part
        endInTar = endInTarCurr;
        endRef = max(endRef, endInRef);
    }

    // Optional: if alignment ends early, append the remainder of target
    if (endInTar < static_cast<int>(target.length()) - 1) {
        result += target.substr(endInTar + 1) + "\n";
    }

    // Optional: write to file
    if (!fileName.empty()) {
        ofstream out(fileName);
        if (out.is_open()) {
            out << result;
            out.close();
        } else {
            cerr << "Could not open file: " << fileName << "\n";
        }
    }

    // Optional: logic for mismatch counter  !!!!!!!!!!!!!!!!!!!!!!!!!!! check later
    // if (trouble_parts > (sub_length * T1)) {
    //     mismatch++;
    // }

    return list.back(); // return the last Position
}


};
int main() {
    // long long start_time = getCPUTime();   UNCOMMENT ON LINUX
    SCCGC reader;
    const int kmer_length = 3;
    const int block_size = 16;
    try {
        string sequence = reader.LocReadSeq("sekvenca_ref.txt");  
        cout << "Meta data: " << reader.meta_data << endl;
        cout << "Local sequence read: " << sequence << endl;
        string sequence_ref = reader.GloReadRefSeq("sekvenca_ref.txt");  
        cout << "Ref sequence: " << sequence_ref<< endl;
        string sequence_tar = reader.GloReadTarSeq("sekvenca_tar.txt", "output.txt");
        cout << "Target sequence: " << sequence_tar << endl;
        cout << "block size " << block_size << endl;
        vector<string> refBlocks = reader.createBlocks(sequence_ref, block_size);
        vector<string> tarBlocks = reader.createBlocks(sequence_tar, block_size);

        cout << "=== Lokalna podudaranja ===\n";
        for (size_t i = 0; i < refBlocks.size(); ++i) {
            reader.createLocalHash(refBlocks[i], kmer_length);
            auto localMatches = reader.localMatch(refBlocks[i], tarBlocks[i], kmer_length);

            cout << "Block " << i << ": Pronadeno " << localMatches.size() << " podudaranja\n";
            for (auto &m : localMatches) {
                cout << "  Ref[" << m.startInRef << "-" << m.endInRef << "] "  
                     << "Tar[" << m.startInTar << "-" << m.endInTar << "]\n";
            }
        }
        cout << "=== Globalna podudaranja ===\n";
        reader.globalHash.clear();
        auto globalMatches = reader.globalMatch(sequence_ref, sequence_tar, kmer_length);
        cout << "Pronadeno " << globalMatches.size() << " globalnih podudaranja\n";
        for (auto &m : globalMatches) {
            cout << "Ref[" << m.startInRef << "-" << m.endInRef
                 << "] Tar[" << m.startInTar << "-" << m.endInTar << "]\n";
        }

        Position lastPosGl = reader.format_matches(globalMatches, sequence_tar, 0, "matches_output.txt");

        
    } catch (const exception& e) {
        cerr << "Greška: " << e.what() << endl;
    }

    // long long end_time = getCPUTime();    UNCOMMENT ON LINUX
    
    // cout << "Execution time " << (end_time - start_time) << endl;  UNCOMMENT ON LINUX 
    system("pause"); 
    return 0;
}
