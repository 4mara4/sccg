#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
using namespace std;

class K_mer { // Dora writing
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
class SCCGC {   //  Dora writing 
public:
    string meta_data;
    int length;
    unordered_map<size_t, vector<K_mer>> localHash;

    string LocReadSeq(const string& filename) {  // reads the sequence for local matching
        ifstream sequence_file(filename);
        if (!sequence_file.is_open()) {
            throw runtime_error("Could not open file: " + filename);
        }
        stringstream ss;
        string line;
        bool is_first_line = true;
        if(getline(sequence_file, line)) {   // the first line is the metadata
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
    string GloReadRefSeq(const string& filename) {  // reads the reference sequence for global matching
    
        ifstream sequence_file(filename);
        if (!sequence_file.is_open()) {
            throw runtime_error("Could not open file: " + filename);
        }
        stringstream ss;
        string line;
        if(!getline(sequence_file, line)) {   // the first line is the metadata
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
    string GloReadTarSeq(const string& inputfilename, const string& outputfilename){ // reads the target sequence for global matching and finds lowercase_mode areas and N areas and returns the rest
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
        getline(infile, meta_data); // first line is metadata
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

                // checking for lowercase_mode
                if (islower(current_ch)) {
                    length_lower++;
                    if (!lowercase_mode) { // the start of a new lowercase_mode sequence
                        lowercase_mode = true;
                        start_L = i + length_read;
                        int offset=start_L - end_L;  // offset from the end of the last lowercase_mode sequence
                        Llist << (offset) << ':';
                    }
                    current_ch = toupper(current_ch);  // now that we have noted where it is lowercase_mode, we convert it to uppercase
                } else {
                    if (lowercase_mode) {
                        Llist << length_lower << ' '; 
                        end_L = start_L + length_lower - 1;
                    }
                    lowercase_mode = false;
                    length_lower = 0;
                }

                // checking for N
                if (current_ch == 'N') { 
                    lenght_N++;
                    if (!N_mode)  {  // the start of a new N sequence
                        N_mode = true;
                        start_N = i + length_read;
                        int offset = start_N - end_N;  // offset from the end of the last N sequence
                        Nlist <<offset << ':';
                    }
                } else {
                    if (N_mode) {
                        Nlist << lenght_N << ' ';
                        end_N = start_N + lenght_N -1;
                    }
                    N_mode = false;
                    lenght_N = 0;

                    clean_seq << current_ch;  // the N characters are left out of the clean sequence
                }
            }
            length_read += line_length;
        }

        if(lowercase_mode) {  // if we were left in lowercase mode at the end of the file
            Llist << length_lower << ' ';
        }
        if(N_mode) {  // if we were left in N mode at the end of the file
            Nlist << lenght_N << ' ';
        }
        
        outfile << Llist.str() << '\n';
        outfile << Nlist.str() << '\n';

        infile.close();
        outfile.close();

        return clean_seq.str();
    }
    
    void createLocalHash(const string& seq, int kmer_length) {
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
};
int main() {
    SCCGC reader;
    try {
        string sequence = reader.LocReadSeq("sekvenca_ref.txt");  
        cout << "Meta data: " << reader.meta_data << endl;
        cout << "Local sequence read: " << sequence << endl;
        string sequence_ref = reader.GloReadRefSeq("sekvenca_ref.txt");  
        cout << "Ref sequence: " << sequence_ref<< endl;
        string sequence_tar = reader.GloReadTarSeq("sekvenca_tar.txt", "output.txt");
        cout << "Target sequence: " << sequence_tar << endl;
        reader.createLocalHash(sequence_ref, 3);  // Example k-mer length of 3
        for (const auto& pair : reader.localHash) {
            size_t key = pair.first;
            const vector<K_mer>& kmers = pair.second;

            cout << "Key: " << key << " -> kmers: ";
            for (const K_mer& kmer_obj : kmers) {
                cout << "[" << kmer_obj.getKmer() << ", start=" << kmer_obj.getStart() << "] ";
            }
            cout << endl;
        }
        
    } catch (const exception& e) {
        cerr << "Greška: " << e.what() << endl;
    }
    system("pause"); 
    return 0;
}
