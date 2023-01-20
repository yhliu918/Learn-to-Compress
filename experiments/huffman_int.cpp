/**
 * This is an implementation of Huffman coding.
 *
 * The core algorithm is taken from the CLR book (Introduction of Algorithms),
 * Chapter 16.3, and directly used to implement the 'build_tree()' routine.
 *
 * After the tree is built, a code table that maps a character to a binary
 * code is built from the tree, and used for encoding text. Decoding is done
 * by traversing the Huffman tree, as prescribed by the algorithm.
 *
 * Binary codes are represented by std::vector<bool>, which is a specialized
 * vector that optimizes space.
 */

#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <string>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include "../headers/common.h"
typedef uint32_t leco_type;
using namespace std;

// A Huffman Tree Node
struct HuffmanTree {
    char c; // character in an alphabet
    int cfreq; // frequency of c.
    struct HuffmanTree *left;
    struct HuffmanTree *right;
    HuffmanTree(char c, int cfreq, struct HuffmanTree *left=NULL,
                struct HuffmanTree *right=NULL) :
        c(c), cfreq(cfreq), left(left), right(right) {
    }
    ~HuffmanTree() {
        delete left, delete right;
    }
    // Compare two tree nodes
    class Compare {
    public:
        bool operator()(HuffmanTree *a, HuffmanTree *b) {
            return a->cfreq > b->cfreq;
        }
    };
};

/**
 * Builds a Huffman Tree from an input of alphabet C, where C is a vector
 * of (character, frequency) pairs.
 */
HuffmanTree *build_tree(vector< pair<char, unsigned> > &alph) {
    // First build a min-heap
    // Build leaf nodes first
    priority_queue<HuffmanTree *, vector<HuffmanTree *>, HuffmanTree::Compare > alph_heap;
    for (vector< pair<char, unsigned> >::iterator it = alph.begin();
         it != alph.end(); ++it) {
        HuffmanTree *leaf = new HuffmanTree(it->first, it->second);
        alph_heap.push(leaf);
    }

    // HuffmanTree algorithm: Merge two lowest weight leaf nodes until
    // only one node is left (root).
    HuffmanTree *root = NULL;
    while (alph_heap.size() > 1) {
        HuffmanTree *l, *r;
        l = alph_heap.top();
        alph_heap.pop();
        r = alph_heap.top();
        alph_heap.pop();
        root = new HuffmanTree(0, l->cfreq + r->cfreq, l, r);
        alph_heap.push(root);
    }

    return root;
}

/**
 * Prints the tree nodes in breadth-first order
 */
void print_tree(HuffmanTree *t) {
    deque< pair<HuffmanTree *, int> > q;

    q.push_back(make_pair(t, 0));
    int curlevel = -1;
    while (!q.empty()) {
        HuffmanTree *parent = q.front().first;
        int level = q.front().second;
        q.pop_front();
        if (curlevel != level) {
            curlevel = level;
            cout << "Level " << curlevel << endl;
        }
        cout << parent->cfreq << " " << parent->c << endl;
        if (parent->left)
            q.push_back(make_pair(parent->left, level + 1));
        if (parent->right)
            q.push_back(make_pair(parent->right, level + 1));
    }
}

typedef vector<bool> code_t;
typedef map<char, code_t> codetable;
/**
 * Makes a lookup table (std::map) of (c -> code) from a HuffmanTree, where
 * code is an unsigned long representing the binary code.
 */
map<char, code_t> build_lookup_table(HuffmanTree *htree, int* size) {
    int totalsize = 0;
    codetable lookup;
    deque< pair<HuffmanTree *, code_t> > q;

    q.push_back(make_pair(htree, code_t()));
    while (!q.empty()) {
        HuffmanTree *node, *lc, *rc;
        code_t code;
        node = q.front().first;
        code = q.front().second;
        q.pop_front();
        lc = node->left;
        rc = node->right;
        if (lc) {
            // HuffmanTree is always full (either no children or two children)
            // Left child is appended a 0 and right child a 1.
            code_t code_cp(code);
            q.push_back(make_pair(lc, (code.push_back(0), code)));
            q.push_back(make_pair(rc, (code_cp.push_back(1), code_cp)));
        } else {
            // Leaf node: contains the character
            lookup.insert(make_pair(node->c, code));
            totalsize+= (8+code.size());
            // cout << "(" << node->c << ", ";
            // for (unsigned i = 0; i < code.size(); i++) {
            //     cout << code[i];
            // }
            // cout << ")" << endl;
        }
    }
    *size = totalsize;
    return lookup;
}


/**
 * Encodes an input string. returns a byte vector.
 */
code_t encode(string input, codetable &lookup) {
    code_t result;

    for (string::iterator it = input.begin(); it != input.end(); ++it) {
            code_t b = lookup[*it];
            result.insert(result.end(), b.begin(), b.end());
    }

    return result;
}

/**
 * Look up the next valid code in @biter using @htree and returns the
 * resulting string. Note the iterator @biter is advanced by the actual
 * length of the next valid code, which varies.
 */
char code_lookup(code_t::iterator &biter, const code_t::iterator &biter_end,
                 const HuffmanTree *htree) {
    const HuffmanTree *node = htree;

    while (true) {
        if (!node->left) {
            // Huffman tree is full: always contains both children or none.
            break;
        }
        if (biter == biter_end) {
            throw std::out_of_range("No more bits");
        }
        if (*biter) {
            node = node->right;
        } else {
            node =node->left;
        }
        ++biter;
    }

    return node->c;
}

/**
 * Decodes a compressed string represented by a bit vector (vector<char>)
 * @compressed, using a HuffmanTree @htree.
 * Returns the original string.
 */
string decode(code_t &compressed, const HuffmanTree *htree) {
    string result;

    code_t::iterator biter = compressed.begin();
    while (true) {
        try {
            result += code_lookup(biter, compressed.end(), htree);
        } catch (const std::out_of_range &oor) {
            // Iterator exhausted.
            break;
        }
    }

    return result;
}

/**
 * Tests
 */
// Make frequency table from a string.
vector< pair<char, unsigned> > make_freq_table(string inp) {
    map<char, unsigned> cfmap;
    vector< pair<char, unsigned> >cfvec;

    for (unsigned i = 0; i < inp.size(); i++) {
        if (cfmap.find(inp[i]) == cfmap.end()) {
            cfmap.insert(make_pair(inp[i], 1));
        }
        cfmap[inp[i]] += 1;
    }

    for (map<char, unsigned>::iterator it = cfmap.begin();
         it != cfmap.end(); ++it) {
        cfvec.push_back(make_pair(it->first, it->second));
    }

    return cfvec;
}

string bitvec_to_string(code_t &bitvec) {
    string result;
    size_t nbits;

    nbits = bitvec.size() & 7;

    // Write the number of "hanging bits" at the first byte
    result += static_cast<char>(nbits); // at most 7

    char byte = 0;
    for (unsigned i = 0; i < bitvec.size(); i++) {
        unsigned boff = i & 7;
        byte |= bitvec[i] << boff;
        if (boff == 7) {
            // Write a byte
            result += byte;
            byte = 0;
        }
    }
    if (nbits) {
        result += byte;
    }

    return result;
}

code_t string_to_bitvec(string packed) {
    code_t result;

    assert(packed.size());
    if (packed.size() == 1) {
        assert(packed[0] == 0);
        return result;
    }
    unsigned nbits = packed[0];
    for (string::iterator it = packed.begin() + 1; it != packed.end(); ++it) {
        for (unsigned i = 0; i < 8; i++) {
            result.push_back((*it >> i) & 1);
        }
    }
    // fix the last byte
    if (nbits) {
        for (unsigned i = 0; i < (8 - nbits); i++) {
            result.pop_back();
        }
    }

    return result;
}

#include <cstdio>
void hexdump(const unsigned char *bytes, int nbytes) {
    int i, j;

    for (i = 0; i < nbytes; i += 16) {
        printf("%06x: ", i);
        for (j = 0; j < 16; j++) {
            if (i + j < nbytes) {
                printf("%02x ", bytes[i + j]);
            } else {
                printf("   ");
            }
        }
        printf(" ");
        for (j = 0; j < 16; j++) 
            if (i + j < nbytes)
                printf("%c", isprint(bytes[i+j]) ? bytes[i+j] : '.');
        printf("\n");
    }
}

int random(int m)
{
    return rand() % m;
}
template <typename T>
static std::vector<T> load_data_binary(const std::string& filename,
    bool print = true) {
    std::vector<T> data;

    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    // Read size.
    uint64_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    data.resize(size);
    // Read values.
    in.read(reinterpret_cast<char*>(data.data()), size * sizeof(T));
    in.close();

    return data;
}


template <typename T>
static std::vector<T> load_data(const std::string& filename) {
    std::vector<T> data;
    std::ifstream srcFile(filename, std::ios::in);
    if (!srcFile) {
        std::cout << "error opening source file." << std::endl;
        return data;
    }

    while (srcFile.good()) {
        T next;
        srcFile >> next;
        if (!srcFile.good()) { break; }
        data.emplace_back(next);

    }
    srcFile.close();

    return data;
}
double getNow() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int main(int argc, const char *argv[]) {
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int model_size = atoi(argv[3]);
    bool binary = atoi(argv[4]);
    std::vector<leco_type> data;
    source_file = "../data/"+source_file;

    if(binary){
        data = load_data_binary<leco_type>(source_file);
    }
    else{
        data = load_data<leco_type>(source_file);
    }
    int N = data.size();
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    }
    int delta = 32;
    std::vector<std::string> block_start_vec;
    std::vector<int> start_index;
    std::vector<int> seglen;
    std::vector<HuffmanTree*> htrees;

    int totalsize = 0;
    double start = getNow();
    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        std::string s((char *)(data.data() + i * block_size),block_length*sizeof(leco_type) );
        vector< pair<char, unsigned> > cfvec = make_freq_table(s);
        int htree_size = 0;
        HuffmanTree *htree = build_tree(cfvec);
        //print_tree(htree);
        codetable ctbl = build_lookup_table(htree, &htree_size);
        code_t t = encode(s, ctbl);
        // cout << "original:" << endl << s << endl;
        // cout << "encoded (compression ratio: "
        //      << (t.size() + 7) / 8 << "/" << s.size() << " or "
        //      << ((float)(t.size() / 8) / s.size()) << "):" << endl;
        std::string packed = bitvec_to_string(t);
        htrees.push_back(htree);
        block_start_vec.push_back(packed);
        totalsize += packed.size() + (htree_size/8);
    }
    double end = getNow();
    double compress_time = end - start;

    double compress_throughput = N * sizeof(leco_type) / (compress_time * 1000000000);

    double origin_size = (sizeof(leco_type) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize)*100.0 / origin_size;

    bool flag = true;
    std::vector<leco_type> recover(data.size());
    double totaltime = 0.0;
    // std::cout << "decompress all!" << std::endl;
    int repeat = 1;

    start = getNow();
    for (int iter = 0; iter < repeat; iter++)
    {
        for (int i = 0; i < blocks; i++)
        {
            int block_length = block_size;
            if (i == blocks - 1)
            {
                block_length = N - (blocks - 1) * block_size;
            }
            // std::string uncomp_str;
            // const char* reover= reinterpret_cast<const char*>(recover + block_size*i);
            // snappy::Uncompress(block_start_vec[i].c_str(), seglen[i], &uncomp_str);
            // Decode
            code_t t1 = string_to_bitvec(block_start_vec[i]);
            string s1 = decode(t1, htrees[i]);

        }
        
    }
    end = getNow();

    totaltime += (end - start);
    double da_ns = totaltime / (repeat * data.size()) * 1000000000;
    // std::cout << da_ns << std::endl;
    // std::cout << "random access decompress!" << std::endl;
    
    std::vector<uint32_t> ra_pos;
    repeat = 1;
    int total_ra = 1000;
    for(int i=0;i<total_ra;i++){
        ra_pos.emplace_back(random(N));
        // ra_pos.push_back(i);
    }
    std::vector<uint32_t> buffer(data.size());

    double randomaccesstime = 0.0;
    start = getNow();
    uint32_t mark = 0;
    for (auto index : ra_pos)
    {
        code_t t1 = string_to_bitvec(block_start_vec[index/block_size]);
        string s1 = decode(t1, htrees[index/block_size]);
        leco_type tmpvalue;
        memcpy(&tmpvalue, s1.c_str() + (index % block_size) * sizeof(leco_type), sizeof(leco_type));
        // std::cout<<tmpvalue<<std::endl;
        mark += tmpvalue;
        // if (data[index] != tmpvalue)
        // {

        //     std::cout << "num: " << index << "true is: " << data[index] << " predict is: " << tmpvalue << std::endl;
        //     flag = false;
        //     std::cout << "something wrong! decompress failed" << std::endl;
        // }
        // if (!flag)
        // {
        //     break;
        // }

    }
    
    end = getNow();
    randomaccesstime += (end - start);
    std::ofstream outfile("fix_log", std::ios::app);
    outfile<<mark<<std::endl;

    double ra_ns = randomaccesstime / total_ra * 1000000000;

    std::cout<<"Huffman "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<cr_model<<" "<<cr_wo_model<<" "<<da_ns<<" "<<ra_ns<<" "<<compress_throughput<<std::endl;

   
}
