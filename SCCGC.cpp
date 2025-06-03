#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <climits>
using namespace std;

/**
 * Function for measuring CPU time on Ubuntu Linux.
 * Dora writing
 */
long long getCPUTime() {
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0) {
        return static_cast<long long>(ts.tv_sec) * 1'000'000'000LL + ts.tv_nsec;
    } else {
        return 0LL;
    }
}


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
    double bad_segment_treshold=0.5;  
    int consecutive_bad_segment_tresh=4;
    int consec_bad_segments=0;



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

    string GloReadTarSeq(const string& inputfilename, const string& outputfilename, 
                            string& meta_data, int& line_length, string& Llist, string& Nlist){
        ifstream infile(inputfilename);
        ofstream outfile(outputfilename, ios::app); 

        if (!infile.is_open()) {
            throw runtime_error("Could not open file: " + inputfilename);
        }
        if (!outfile.is_open()) {
            throw runtime_error("Could not open file: " + outputfilename);
        }

        string line;
        getline(infile, meta_data);
        
        int length_read = 0, length_lower = 0, length_N = 0, end_L = 0, end_N = 0, start_L = 0, start_N = 0;
        bool lowercase_mode = false, N_mode = false;
        string clean_seq;

        bool first_data_line = true;
        while (getline(infile, line)) {
            if(first_data_line) {
                line_length = static_cast<int>(line.size());
                first_data_line = false;
            }
            int l = line.size();

            for (int i = 0; i < l; ++i) {
                char current_ch = line[i];

                if (islower(current_ch)) {
                    length_lower++;
                    if (!lowercase_mode) {
                        lowercase_mode = true;
                        start_L = i + length_read;
                        Llist += to_string(start_L - end_L) + ":";
                    }
                    current_ch = toupper(current_ch);
                } else if (lowercase_mode) {
                    lowercase_mode = false;
                    end_L = start_L + length_lower;
                    Llist += to_string(length_lower) + " ";
                    length_lower = 0;
                }
                if(current_ch == 'N') {
                    if(!N_mode) {
                        N_mode = true;
                        start_N = length_read + i;
                        Nlist += to_string(start_N - end_N) + ":";
                    }
                    length_N++;
                } else if (N_mode) {
                    N_mode = false;
                    end_N = start_N + length_N;
                    Nlist += to_string(length_N) + " ";
                    length_N = 0;
                }
                if(current_ch != 'N') {
                    clean_seq += current_ch;
                }
            }
            length_read += l;
        }
        if(lowercase_mode) {
            Llist += to_string(length_lower) + " ";
        }
        if(N_mode) {
            Nlist += to_string(length_N) + " ";
        }
        return clean_seq;
    }
    
    /**
     * @brief Function for writing header data to final file.
     *
     * Writes meta data, line length, list of lowercase and N positions.
     * Marija writing
     */
    void writePreamble(const string &finalFile,
                   const string &meta_data,
                   int line_length,
                   const string &Llist,
                   const string &Nlist) {
        ofstream out(finalFile, ios::trunc);
        if (!out) throw runtime_error("Cannot open " + finalFile);

        out << meta_data << "\n";
        out << line_length << "\n";
        if (Llist.empty()) {
            out << "x\n";
        } else {
            out << Llist << "\n";
        }

        if (Nlist.empty()) {
            out << "x\n";
        } else {
            out << Nlist << "\n";
        }
        out << "\n";
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
            if(kmer == string(kmer_length, 'N')) {                         
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

    /**
     * @brief Finds local alignments between a reference block and a target block.
     *
     * This function scans the target sequence one k-mer at a time. For each k-mer,
     * it looks up all matching positions in the localHash of the reference.
     * It then extends the match forward as long as the bases remain identical,
     * chooses the best (longest, or earliest) match, and records its start/end
     * coordinates in both reference and target. Matched regions are skipped over
     * so they are not re-examined, ensuring each position is aligned at most once.
     *
     * @param reference    A substring of the reference to search within.
     * @param target       A substring of the target to align.
     * @param kmer_length  Length of each k-mer used for initial matching.
     * @return A vector of Position structs defining the aligned segments.
     *
     * Dora + Marija writing
     */

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

    /**
     * @brief Finds global alignments between the full reference and target sequences.
     *
     * This function first builds a hash of all k-mers over the entire reference.
     * It then scans the target sequence one k-mer at a time, looking up each k-mer
     * in the reference hash. For each match, it extends forward as far as the bases
     * remain identical, selects the longest (or earliest on tie) extension, and records
     * its start/end coordinates in both reference and target. Matched regions are skipped
     * so that each target position is aligned only once.
     *
     * @param ref   The full reference sequence.
     * @param tar   The full target sequence to align against the reference.
     * @param klen  Length of each k-mer used for seeding matches.
     * @return      A vector of Position structs describing all alignment segments.
     *
     * Dora + Marija writing
     */

    vector<Position> globalMatch(const string& ref,
                             const string& tar,
                             int klen,
                             int limit = 100) {
   
        const int PRAG_ZA_MALI_TAR = 100000;  

        int effLimit;
        if ((int)tar.size() >= PRAG_ZA_MALI_TAR) {
            effLimit = limit; 
        } else {
            effLimit = (int)ref.size();
        }

        createGlobalHash(ref, klen);
        vector<Position> out;
        int idx = 0;
        int lastEIR = 0;

        while (idx + klen <= (int)tar.size()) {
            string k = tar.substr(idx, klen);
            size_t h = hash<string>{}(k);
            auto it = globalHash.find(h);
            if (it == globalHash.end()) {
                idx++;
                continue;
            }

            int bestS = INT_MAX;
            int bestE = 0;
            bool foundInLimit = false;
            for (auto &km : it->second) {
                if (km.getKmer() != k)
                    continue;

                if (!out.empty() && abs(km.getStart() - lastEIR) > effLimit)
                    continue;

                int ext = 0;
                int r = km.getStart() + klen;
                int t = idx + klen;
                while (r < (int)ref.size() && t < (int)tar.size() && ref[r] == tar[t]) {
                    ext++; r++; t++;
                }
                if (ext > bestE || (ext == bestE && km.getStart() < bestS)) {
                    bestE = ext;
                    bestS = km.getStart();
                    foundInLimit = true;
                }
            }
            if (!foundInLimit) {
                for (auto &km : it->second) {
                    if (km.getKmer() != k)
                        continue;

                    int ext = 0;
                    int r = km.getStart() + klen;
                    int t = idx + klen;
                    while (r < (int)ref.size() && t < (int)tar.size() && ref[r] == tar[t]) {
                        ext++; r++; t++;
                    }
                    if (ext > bestE || (ext == bestE && km.getStart() < bestS)) {
                        bestE = ext;
                        bestS = km.getStart();
                    }
                }
            }
            if (bestS != INT_MAX) {
                Position p;
                p.startInRef = bestS;
                p.endInRef   = bestS + klen + bestE - 1;
                p.startInTar = idx;
                p.endInTar   = idx + klen + bestE - 1;
                out.push_back(p);

                lastEIR = p.endInRef;
                idx += klen + bestE + 1;
            } else {
                idx++;
            }
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
    Position format_matches(const vector<Position>& list, const string& target, const string& fileName = "",bool local=true) {
        if(list.empty()) {
            if(!fileName.empty()) {
                ofstream out(fileName, ios::app);
                if (out) out << target << "\n\n";
            }
            return Position{0,0,0,0};
        }
        string result;
        int trouble_parts = 0, endRef = 0, endInTar = 0;
        bool bad_Segment=false;

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
                    result += pre_allign + "\n";
                    trouble_parts += pre_allign.length();
                }
                result += to_string(startInRef ) + "," + to_string(endInRef ) + "\n";
                continue;
            }

            int endInTarPrev = list[i - 1].endInTar;
            string mismatch = target.substr(endInTarPrev + 1, startInTar - endInTarPrev - 1);	

            if (!mismatch.empty()) {
                result += mismatch + "\n";  
                trouble_parts += mismatch.length();
            }

            result += to_string(startInRef ) + "," + to_string(endInRef ) + "\n";
            endInTar = endInTarCurr;
            endRef = max(endRef, endInRef);
        }

        if (endInTar < static_cast<int>(target.length()) - 1) {
            result += target.substr(endInTar + 1) + "\n";
        }

        if (!fileName.empty()) {
            ofstream out;
            if(local){out.open(fileName, ios::app); }
            else{out.open(fileName); }
            if (out.is_open()) {
                out << result ;
                out.close();
            } else {
                cerr << "Could not open file: " << fileName << "\n";
            }
        }

        if (trouble_parts > (target.size() * bad_segment_treshold)) {
            consec_bad_segments++;
        } else {
            consec_bad_segments=0;
        }

        return list.back();
    }

    /**
     * @brief Post-processes alignment data into delta-encoded segments.
     *
     * This function reads the input file `inFile` line by line, where aligned regions
     * and mismatches are mixed in the following format:
     *   - Lines with "start,end" specify aligned segments.
     *   - Lines without a comma represent literal mismatch sequences.
     *
     * Adjacent aligned segments (where the end of one segment + 1 == start of the next)
     * are merged into a single continuous segment. For each merged segment, it outputs:
     *   deltaB,len
     * where:
     *   - deltaB = the segment’s start position minus the end position of the previous segment
     *              (or the segment’s start if there is no previous segment),
     *   - len    = the length of the segment (end − start + 1).
     * Mismatch lines are passed through unchanged.
     *
     * @param inFile  Path to the interim file containing raw alignments (e.g., "interim.txt").
     * @param outFile Path to the final output file where delta-encoded segments and
     *                mismatch lines will be appended.
     *
     * Dora + Marija writing
     */
    void postProcess(const string& inFile, const string& outFile) {
        ifstream in(inFile);
        ofstream out(outFile, ios::app);
        vector<pair<int,int>> segs;
        string line;
        int prevEnd = 0;

        auto flushSegs = [&](){
            for (auto &s : segs) {
                int b = s.first, e = s.second;
                int deltaB =b - prevEnd + 1;
                if(prevEnd ==0 ) {deltaB--;}
                int len    = e - b + 1;
                out << deltaB << "," << len << "\n";
                prevEnd = b + len;
            }
            segs.clear();
        };

        while (getline(in, line)) {
            if (line.empty()) continue;

            if (line.find(',') != string::npos) {
                int b,e; char c;
                istringstream iss(line);
                iss >> b >> c >> e;
                if (!segs.empty() && segs.back().second + 1 == b) {
                    segs.back().second = e;
                } else {
                    segs.emplace_back(b,e);
                }
            } else {
                flushSegs();
                out << line << "\n";
            }
        }
        flushSegs();
    }
};

int main(int argc, char* argv[]) {
    long long start_time = getCPUTime(); 
    if (argc != 3) {
        std::cerr << "Upotreba: " << argv[0] << " <tar file name> <ref file name\n";
        return 1;
    }
    ofstream cleaner("local_matches_output.txt");
    cleaner.close();
    ofstream cleaner2("global_matches_output.txt");
    cleaner2.close();
    
    const string tempFile  = "interim.txt";
    const string finalFile = "final.txt";
    const string refFile="test/" + string(argv[2]);
    const string tarFile="test/" + string(argv[1]);
    ofstream(tempFile).close();    
    ofstream(finalFile).close();  


    SCCGC reader;
    const int kmer_length = 21;
    const int block_size = 30000;
    auto local = true;
    string meta;
    int line_length;
    string Llist, Nlist;
    try {
        string sequence = reader.LocReadSeq(refFile);  
        string sequence_ref = reader.GloReadRefSeq(refFile);  
        string sequence_tar = reader.GloReadTarSeq(tarFile, "output.txt", meta, line_length, Llist, Nlist);
        if(sequence_tar.length() < block_size * 5) {
            local = false;
        }
        
        vector<string> refBlocks = reader.createBlocks(sequence_ref, block_size);
        vector<string> tarBlocks = reader.createBlocks(sequence_tar, block_size);
        reader.writePreamble(tempFile, meta, line_length, Llist, Nlist);

        cout << "=== Lokalna podudaranja ===\n";
        for (size_t i = 0; i < min(refBlocks.size(), tarBlocks.size()) && local; ++i) {  
            reader.createLocalHash(refBlocks[i], kmer_length);
            vector<Position> localMatches = reader.localMatch(refBlocks[i], tarBlocks[i], kmer_length);

            int block_offset = static_cast<int>(i) * block_size;
            for (auto &m : localMatches) {
                m.startInRef += block_offset;
                m.endInRef   += block_offset;
            }
            Position lastPosLoc = reader.format_matches(localMatches, tarBlocks[i], tempFile, true);
            if(reader.consec_bad_segments >= reader.consecutive_bad_segment_tresh){
                cout << "moramo prec na globalno, ne idemo dalje lokalno"<< endl;
                local = false;
                break;
            }
            
        }

        if(local) {
            reader.postProcess(tempFile, finalFile);
        } else {
            ofstream(tempFile).close();
            reader.writePreamble("interim.txt", meta, line_length, Llist, Nlist);
            cout << "\n=== Globalna podudaranja ===\n";
            reader.globalHash.clear();
            auto globalMatches = reader.globalMatch(sequence_ref, sequence_tar, kmer_length);
            Position lastPosLoc = reader.format_matches(globalMatches, sequence_tar, tempFile);
            reader.postProcess(tempFile, finalFile); //ja
        }
        /*cout << "Meta: " << meta;
        cout << "Length: " << line_length;
        cout << "L list: " << Llist;
        cout << "N list: " << Nlist;*/
        
        
    } catch (const exception& e) {
        cerr << "Greška: " << e.what() << endl;
    }

    long long end_time = getCPUTime(); 
    
    cout << "Execution time " << (end_time - start_time) << endl;
    return 0;
}
