#include <vector>

#include "common.h"
#include "../libs/bit-buffer/src/bit_buffer.hpp"
#include "Tree.hpp"
#include "io/file_io.hpp"

/* Binary one-hot encoding of nucleotides including ambiguities, following the IUPAC code
 * https://www.bioinformatics.org/sms/iupac.html
 */ 
enum class IUPAC_DNA : uint8_t {
    A = 0b1000,
    C = 0b0100,
    T = 0b0010,
    G = 0b0001,
    R = A | G,
    Y = C | T,
    S = G | C,
    W = A | T,
    K = G | T,
    M = A | C,
    B = C | G | T,
    D = A | C | T,
    V = A | C | G,
    N = A | C | T | G,
    GAP = 0b0000,
    INVALID = 0b11111111,
    EMPTY = 0b0000
};

IUPAC_DNA operator|(const IUPAC_DNA lhs, const IUPAC_DNA rhs);
IUPAC_DNA operator&(const IUPAC_DNA lhs, const IUPAC_DNA rhs);
bool is_single_state(const IUPAC_DNA state);
char iupac_dna2char(IUPAC_DNA code);
IUPAC_DNA char2iupac_dna(char code);
bool valid(const IUPAC_DNA state);

class CompressedSequence {
    public:
    CompressedSequence(const std::string& sequence);
    std::string to_string();
    char operator[](size_t idx);
    size_t length() const;
    size_t compressed_size_in_bytes() const;

    private:
    bit_buffer _sequence;
    size_t _sequence_length = 0;

    void from_string(const std::string& sequence);
    size_t position(size_t pos) { return pos * 4; }
    static const unsigned char STATE_ENCODING_LENGTH;
};

class PllTree : public BasicTree {
    public:
    typedef pll_unode_t * Node;
    enum class TRAVERSAL_MODE {
        PREORDER,
        POSTORDER,
        INORDER
    };

    PllTree(const pll_utree_t * tree);
    PllTree(const std::string& newick_string);

    // The following method is inefficient, copying the tree more times
    // than it needs to. But I'm using it for testing only.
    static PllTree load_from_file(const std::string& file_name);

    Node root() const;
    bool is_root(Node node) const;
    void traverse(std::function<bool(Node)> callback_on_node_visit, TRAVERSAL_MODE mode);
    
    // These are put here instead of in spearate node classes to avoid recreating
    // Node objects all the time for performance reasons.
    static bool is_leaf(Node node);
    static bool is_inner_node(Node node);
    static Node left_child(Node node);
    static Node right_child(Node node);
    static Node parent(Node node);
    static size_t node_id(Node node);
    static std::string label(Node node);

    private:
    pll_utree_t * const _tree;
    pll_unode_t * const _tree_root;
};

/* 
 * Algorithms and datastructure for tree based MSA compression.
 * The best compression is archived using the most parsimonous phylogenetic tree.
 * Should be deallocated once compression/decompression is complete to free memory.
 * 
 * Note: This contains code for tree traversal and other stuff that should from an
 * encapsulation standpoint be put into the Tree class. It has been decided, however,
 * that, as compression is only useful in special cases, all it's code should be
 * separated to no increase the main code's complexity.
 */
class MSATreeCompression {
    public:
    class NodeStates {
        public:
        NodeStates(size_t num_nodes);
        IUPAC_DNA operator[](PllTree::Node node) const;
        IUPAC_DNA& operator[](PllTree::Node node);
        size_t size() { return _states.size(); }
        static IUPAC_DNA choose_random_state(const IUPAC_DNA states);
        IUPAC_DNA choose_random_state(PllTree::Node node);
        void reset();
        bool check_all_single_state();

        private:
        std::vector<IUPAC_DNA> _states;
    };

    MSATreeCompression(const Tree& tree);
    void compress(const MSA& msa);
    bool working_buffers_allocated();
    void free_working_buffers();
    const size_t num_sequences;

    private:
    PllTree _tree;
    std::shared_ptr<CompressedSequence> _root_sequence = nullptr;
    std::shared_ptr<std::vector<std::shared_ptr<CompressedSequence>>> _node_sequences = nullptr;
    std::shared_ptr<bit_buffer> _change_encoding = nullptr;
    const size_t _root_sequence_idx = 0;

    void build_ancestral_states(const MSA& msa, size_t position, MSATreeCompression::NodeStates& node_states);
};

/*
 * Container for in-memory compressed MSA, handles compression, decompression and memory management.
 */
class CompressedMSA {
    public:
    CompressedMSA(const MSA& msa, const Tree& tree);
    const MSA& decompress();

    private:
    MSATreeCompression _tree_compression;
    void compress();
};