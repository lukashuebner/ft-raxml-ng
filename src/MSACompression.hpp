#include <vector>
#include <memory>

#include "common.h"
#include "../libs/bit-buffer/src/bit_buffer.hpp"
#include "Tree.hpp"
#include "io/file_io.hpp"

/* 
 * Binary one-hot encoding of nucleotides including ambiguities, following the IUPAC code
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
    D = A | G | T,
    H = A | C | T,
    V = A | C | G,
    N = A | C | T | G,
    GAP = 0b0000,
    INVALID = 0b11111111,
    EMPTY = 0b0000,
};
const uint8_t IUPAC_DNA_CODE_LENGTH = 4;

/* 
 * Expanded version of IUPAC DNA Coding where each ambiguous state has its own one-hot encoding bit.
 * This enables us to encode for example M | H, as we need to build the ancestral states if ambiguities
 * are allowed in the input sequences.
 */ 
enum class IUPAC_DNA_16BIT : uint32_t {
    A =     0b1000000000000000,
    C =     0b0100000000000000,
    T =     0b0010000000000000,
    G =     0b0001000000000000,
    R =     0b0000100000000000,
    Y =     0b0000010000000000,
    S =     0b0000001000000000,
    W =     0b0000000100000000,
    K =     0b0000000010000000,
    M =     0b0000000001000000,
    B =     0b0000000000100000,
    D =     0b0000000000010000,
    H =     0b0000000000001000,
    V =     0b0000000000000100,
    N =     0b0000000000000010,
    GAP =   0b0000000000000001, // It's a state after all
    INVALID = 0xFFFFFFFF,
    EMPTY = 0x00000000
};
const uint8_t IUPAC_DNA_16BIT_CODE_LENGTH = 32;

IUPAC_DNA operator|(const IUPAC_DNA lhs, const IUPAC_DNA rhs);
IUPAC_DNA operator&(const IUPAC_DNA lhs, const IUPAC_DNA rhs);
IUPAC_DNA operator^(const IUPAC_DNA lhs, const IUPAC_DNA rhs);
std::ostream& operator<<(std::ostream& os, const IUPAC_DNA& rhs);
std::ostream& operator<<(std::ostream& os, const IUPAC_DNA_16BIT& rhs);

IUPAC_DNA_16BIT operator|(const IUPAC_DNA_16BIT lhs, const IUPAC_DNA_16BIT rhs);
IUPAC_DNA_16BIT operator&(const IUPAC_DNA_16BIT lhs, const IUPAC_DNA_16BIT rhs);
IUPAC_DNA_16BIT operator^(const IUPAC_DNA_16BIT lhs, const IUPAC_DNA_16BIT rhs);

IUPAC_DNA iupac_dna_to_4bit (const IUPAC_DNA_16BIT & state);

bool is_single_state(const IUPAC_DNA state);
bool is_single_state(const IUPAC_DNA_16BIT state);

char iupac_dna2char(IUPAC_DNA code);
std::string iupac_dna2string_16bit(IUPAC_DNA_16BIT code);

IUPAC_DNA char2iupac_dna(char code);
IUPAC_DNA_16BIT char2iupac_dna_16bit(char code);

bool valid(const IUPAC_DNA state);
bool valid(const IUPAC_DNA_16BIT state);


class CompressedSequence {
    public:
    CompressedSequence();
    CompressedSequence(const std::string& sequence);

    operator std::string() { return move(to_string()); }
    std::string to_string();

    CompressedSequence slice(size_t start, size_t length);
    
    char operator[](size_t idx);
    IUPAC_DNA get_code(size_t idx);
    
    size_t length() const;
    size_t compressed_size_in_bytes() const;

    void append(IUPAC_DNA code);
    void append(char c);
    void append(CompressedSequence & sequence);

    private:
    bit_buffer _sequence;
    size_t _sequence_length = 0;

    void from_string(const std::string& sequence);
    size_t position(size_t pos) const { return pos * 4; }
    static const unsigned char STATE_ENCODING_LENGTH;
};

std::ostream& operator<<(std::ostream& os, CompressedSequence & rhs);

class PllTree : public BasicTree {
    public:
    enum class TRAVERSAL_MODE {
        PREORDER,
        POSTORDER,
        INORDER
    };
    
    // Wrapper for a pll_rnode_t * and respective C functions
    // A user should never create a Node object themselfs, let the PllTree object handle this for you.
    // The life of a Node object ends when the corresponding PllTree is destructed.
    class Node {
        public:
        Node(const pll_rnode_t * pll_node, PllTree * tree);

        bool is_leaf() const;
        bool is_inner_node() const;
        bool is_root() const;
        
        Node * left_child();
        Node * right_child();
        Node * parent();
        size_t node_id() const;
        std::string label() const;

        friend bool operator==(const Node & lhs, const Node & rhs);
        friend bool operator!=(const Node & lhs, const Node & rhs);

        private:
        const pll_rnode_t * const _pll_node;
        PllTree * const _tree = nullptr;
    };

    PllTree(const pll_rtree_t * tree);
    PllTree(const std::string& newick_string);
    ~PllTree();

    Node * root();
    void traverse(std::function<bool(Node *)> callback_on_node_visit, TRAVERSAL_MODE mode);
    
    // A rooted tree has one more inner node than a unrooted tree with the same amount of nodes.
    virtual size_t num_inner() const { return _num_tips - 1; };
    // The same goes for the number of branches
    virtual size_t num_branches() const { return _num_tips ? _num_tips + _num_tips - 2 : 0; };

    private:
    // This function will cache Node objects and thus avoid recreating them each time they're requested.
    Node * create_Node_from_pll_node(const pll_rnode_t * pll_node);

    const pll_rtree_t * const _tree;
    std::vector<Node *> _node_ptrs;
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
        // Constructors
        // num_nodes = 0 is not allowed
        NodeStates(size_t num_nodes);

        // The following two methods use the node id of the passed Node object as index
        // It is assumed that the node id is in [0,..,num_nodes], do not use clv_index
        IUPAC_DNA_16BIT operator[](const PllTree::Node & node) const;
        IUPAC_DNA_16BIT& operator[](const PllTree::Node & node);
        IUPAC_DNA_16BIT operator[](const size_t node_id) const;
        IUPAC_DNA_16BIT& operator[](const size_t node_id);

        // Returns the number of Node's states in this NodeStates object (some of which may not be used)
        size_t size() const { return _states.size(); }

        // Choose a single state from the given ambiguous state
        // If a non-ambiguous state (including GAP) is given, it will be returned
        // Throws an exception if passed INVALID
        // Does not guarantee even pseudo-randomness
        static IUPAC_DNA choose_random_state(const IUPAC_DNA states);
        static IUPAC_DNA_16BIT choose_random_state(const IUPAC_DNA_16BIT states);

        // Sets the Nodes state to a random state which is covered by it's current (ambiguous)
        // state description. I.e. if the node's state is W = A | T, it's state will be set
        // to either A or T but not C org G.
        // Throws an exception if the Node's state is INVALID
        // Does not guarantee even pseudo-randomness
        IUPAC_DNA_16BIT choose_random_state(const PllTree::Node & node);

        // Sets all Node's states back to INVALID
        void reset();
        
        // Returns true iff each Node's state in this object a single state
        // Throws if any state is INVALID or GAP
        // EMTPY states will cause this function to return false
        bool check_all_single_state() const;

        // Support iteration using range-based for loops
        typedef IUPAC_DNA_16BIT * iterator;
        typedef const IUPAC_DNA_16BIT * const_iterator;
        
        iterator begin() { return &_states[0]; }
        const_iterator begin() const { return &_states[0]; }
        iterator end() { return &_states[_states.size()]; }
        const_iterator end() const { return &_states[_states.size()]; }

        private:
        std::vector<IUPAC_DNA_16BIT> _states;
    };

    class ChangeEncoding {
        public:
        ChangeEncoding(size_t max_edge_id) : _max_edge_id(max_edge_id) {}

        // Return the encoding in the least significant bits of the returned bytes
        // The remaining higher bits are set to 0
        virtual uint64_t encoding() const = 0;

        // Returns the number of bits in the encoding which are significant
        uint8_t bit_length() const;
        static uint8_t bit_length(size_t max_edge_id);

        // Trivial getters
        virtual size_t edge_id() const = 0;
        size_t max_edge_id() const { return _max_edge_id; }

        private:
        const size_t _max_edge_id;

        uint8_t bits_for_edge_id() const;
        static uint8_t bits_for_edge_id(size_t max_edge_id);
    };

    // Base Class for nucleotide encoding and decoding
    class ChangeEncoder : public ChangeEncoding {
        public:
        // Used to encode a change, that is construct a ChangeEncoding given the unencoded 
        // edge id, from state and to state.
        // The max_edge_id is used to determine the number of bits used for encoding edge_ids.
        ChangeEncoder(size_t edge_id, size_t max_edge_id, IUPAC_DNA from_state, IUPAC_DNA to_state);
        
        // Return the encoding in the least significant bits of the returned bytes
        // The remaining higher bits are set to 0
        uint64_t encoding() const;

        // Trivial getters
        virtual size_t edge_id() const { return _edge_id; }
        IUPAC_DNA from_state() const { return _from_state; }
        IUPAC_DNA to_state() const { return _to_state; }

        private:
        const size_t _edge_id;
        const IUPAC_DNA _from_state;
        const IUPAC_DNA _to_state;
    };

    // Encodes a nucleotide as a tuple of edge_id and nucleotide change.
    // We use the node id as the edge id of the edge that is leading to the node in
    // the (rooted and directed) compression tree.
    // Can handle ambiguous nucleotide states and gaps.
    class ChangeDecoder : public ChangeEncoding {
        public:
        // Used to decode a change, that is construct a ChangeEncoding given the encoded
        // bit-sequence.
        // The max_edge_id is used to determine the number of bits used for encoding edge_ids.
        ChangeDecoder(uint64_t encoding, size_t max_edge_id);

        // Apply change to the given state and resturns the result
        // apply_change ° apply_change = Identity
        IUPAC_DNA apply_change(const IUPAC_DNA state) const;
        
        virtual uint64_t encoding() const { return _encoding; };
        virtual size_t edge_id() const;
        
        private:
        const uint64_t _encoding;
    };

    MSATreeCompression(const Tree& tree);
    MSATreeCompression(const std::string newick_string_rooted);
    
    // Takes a reference to a MSA which matches the given tree and compresses it using
    // tree-based MSA compression (the root sequence and the changes along the edges
    // are stored).
    // Returns the compressed sites in bits (for now, the length of the change encoding
    // plus the length of the root sequence encoding).
    size_t compress(const MSA& msa);
    
    // TODO support partitioned MSAs
    // Builds and populates a MSA with the stored data
    std::shared_ptr<MSA> decompress(RangeList& range_list);

    // Returns the number of sequences stored in the MSA, that is
    // the numver of tips of the compression tree.
    size_t num_sequences() const { return _num_sequences; }

    uint64_t build_ancestral_states(const MSA& msa, size_t position, MSATreeCompression::NodeStates& node_states);

    private:
    PllTree _tree;
    const size_t _num_sequences;
    std::shared_ptr<CompressedSequence> _root_sequence = nullptr;
    std::shared_ptr<bit_buffer> _encoding_of_changes = nullptr;
    std::shared_ptr<std::vector<size_t>> _site_idx_in_encoding = nullptr;

    MSATreeCompression();
};