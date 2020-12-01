#include <sstream>
#include <exception>
#include <algorithm>

#include "MSACompression.hpp"

using namespace std;

IUPAC_DNA operator|(const IUPAC_DNA lhs, const IUPAC_DNA rhs) {
    return static_cast<IUPAC_DNA> (
        static_cast<std::underlying_type<IUPAC_DNA>::type>(lhs) |
        static_cast<std::underlying_type<IUPAC_DNA>::type>(rhs)
    );
}

IUPAC_DNA operator&(const IUPAC_DNA lhs, const IUPAC_DNA rhs) {
    return static_cast<IUPAC_DNA> (
        static_cast<std::underlying_type<IUPAC_DNA>::type>(lhs) &
        static_cast<std::underlying_type<IUPAC_DNA>::type>(rhs)
    );
}

IUPAC_DNA operator^(const IUPAC_DNA lhs, const IUPAC_DNA rhs) {
    return static_cast<IUPAC_DNA> (
        static_cast<std::underlying_type<IUPAC_DNA>::type>(lhs) ^
        static_cast<std::underlying_type<IUPAC_DNA>::type>(rhs)
    );
}

IUPAC_DNA_16BIT operator|(const IUPAC_DNA_16BIT lhs, const IUPAC_DNA_16BIT rhs) {
    return static_cast<IUPAC_DNA_16BIT> (
        static_cast<std::underlying_type<IUPAC_DNA_16BIT>::type>(lhs) |
        static_cast<std::underlying_type<IUPAC_DNA_16BIT>::type>(rhs)
    );
}

IUPAC_DNA_16BIT operator&(const IUPAC_DNA_16BIT lhs, const IUPAC_DNA_16BIT rhs) {
    return static_cast<IUPAC_DNA_16BIT> (
        static_cast<std::underlying_type<IUPAC_DNA_16BIT>::type>(lhs) &
        static_cast<std::underlying_type<IUPAC_DNA_16BIT>::type>(rhs)
    );
}

IUPAC_DNA_16BIT operator^(const IUPAC_DNA_16BIT lhs, const IUPAC_DNA_16BIT rhs) {
    return static_cast<IUPAC_DNA_16BIT> (
        static_cast<std::underlying_type<IUPAC_DNA_16BIT>::type>(lhs) ^
        static_cast<std::underlying_type<IUPAC_DNA_16BIT>::type>(rhs)
    );
}

ostream& operator<<(ostream& os, const IUPAC_DNA& rhs) {
    os << iupac_dna2char(rhs);
    return os;
}

std::ostream& operator<<(std::ostream& os, const IUPAC_DNA_16BIT& rhs) {
    os << iupac_dna2string_16bit(rhs);
    return os;
}

std::ostream& operator<<(std::ostream& os, CompressedSequence & rhs) {
    os << rhs.to_string();
    return os;
}

IUPAC_DNA iupac_dna_to_4bit (const IUPAC_DNA_16BIT & state) {
       switch(state) {
        case IUPAC_DNA_16BIT::A:
            return IUPAC_DNA::A;
        case IUPAC_DNA_16BIT::C:
            return IUPAC_DNA::C;
        case IUPAC_DNA_16BIT::T:
            return IUPAC_DNA::T;
        case IUPAC_DNA_16BIT::G:
            return IUPAC_DNA::G;
        case IUPAC_DNA_16BIT::R:
            return IUPAC_DNA::R;
        case IUPAC_DNA_16BIT::Y:
            return IUPAC_DNA::Y;
        case IUPAC_DNA_16BIT::S:
            return IUPAC_DNA::S;
        case IUPAC_DNA_16BIT::W:
            return IUPAC_DNA::W;
        case IUPAC_DNA_16BIT::K:
            return IUPAC_DNA::K;
        case IUPAC_DNA_16BIT::M:
            return IUPAC_DNA::M;
        case IUPAC_DNA_16BIT::B:
            return IUPAC_DNA::B;
        case IUPAC_DNA_16BIT::D:
            return IUPAC_DNA::D;
        case IUPAC_DNA_16BIT::H:
            return IUPAC_DNA::H;
        case IUPAC_DNA_16BIT::V:
            return IUPAC_DNA::V;
        case IUPAC_DNA_16BIT::N:
            return IUPAC_DNA::N;
        case IUPAC_DNA_16BIT::GAP:
            return IUPAC_DNA::GAP;
        case IUPAC_DNA_16BIT::INVALID:
        default:
            throw new runtime_error("Invalid code");
            return IUPAC_DNA::INVALID; // Silence "control reaches end of non-void function"
    } 
}

char iupac_dna2char(IUPAC_DNA code) {
    switch(code) {
        case IUPAC_DNA::A:
            return 'A';
        case IUPAC_DNA::C:
            return 'C';
        case IUPAC_DNA::T:
            return 'T';
        case IUPAC_DNA::G:
            return 'G';
        case IUPAC_DNA::R:
            return 'R';
        case IUPAC_DNA::Y:
            return 'Y';
        case IUPAC_DNA::S:
            return 'S';
        case IUPAC_DNA::W:
            return 'W';
        case IUPAC_DNA::K:
            return 'K';
        case IUPAC_DNA::M:
            return 'M';
        case IUPAC_DNA::B:
            return 'B';
        case IUPAC_DNA::D:
            return 'D';
        case IUPAC_DNA::H:
            return 'H';
        case IUPAC_DNA::V:
            return 'V';
        case IUPAC_DNA::N:
            return 'N';
        case IUPAC_DNA::GAP:
            return '-';
        case IUPAC_DNA::INVALID:
        default:
            throw new runtime_error("Invalid code");
            return ' '; // Silence "control reaches end of non-void function"
    }
}

string iupac_dna2string_16bit(IUPAC_DNA_16BIT code) {
    auto code_bits = static_cast<underlying_type<IUPAC_DNA_16BIT>::type>(code);

    if (!valid(code)) {
        return string("(Invalid)");
    }

    string str("(");
    while (code_bits) {
        // ( x & ~(x-1) ) will *return* the lowest bit set
        underlying_type<IUPAC_DNA_16BIT>::type lowest_bit = code_bits & ~(code_bits - 1);

        // Append the character representing this bit
        str.push_back(
            iupac_dna2char(
                iupac_dna_to_4bit(
                    static_cast<IUPAC_DNA_16BIT>(lowest_bit)
        )));

        // Clear the lowest bit set
        code_bits ^= lowest_bit;

        // Append a , if this was not the last possibility of the ambiguity
        if (code_bits) {
            str.push_back(',');
        }
    }
    str.push_back(')');

    return str;
}

IUPAC_DNA char2iupac_dna(char code) {
    switch(code) {
        case 'A':
            return IUPAC_DNA::A;
        case 'C':
            return IUPAC_DNA::C;
        case 'T':
            return IUPAC_DNA::T;
        case 'G':
            return IUPAC_DNA::G;
        case 'R':
            return IUPAC_DNA::R;
        case 'Y':
            return IUPAC_DNA::Y;
        case 'S':
            return IUPAC_DNA::S;
        case 'W':
            return IUPAC_DNA::W;
        case 'K':
            return IUPAC_DNA::K;
        case 'M':
            return IUPAC_DNA::M;
        case 'B':
            return IUPAC_DNA::B;
        case 'D':
            return IUPAC_DNA::D;
        case 'H':
            return IUPAC_DNA::H;
        case 'V':
            return IUPAC_DNA::V;
        case 'N':
            return IUPAC_DNA::N;
        case '-':
        case '?':
            return IUPAC_DNA::GAP;
        default:
            throw runtime_error("Invalid nucleotide code: " + code);
    }
}

IUPAC_DNA_16BIT char2iupac_dna_16bit(char code) {
    switch(code) {
        case 'A':
            return IUPAC_DNA_16BIT::A;
        case 'C':
            return IUPAC_DNA_16BIT::C;
        case 'T':
            return IUPAC_DNA_16BIT::T;
        case 'G':
            return IUPAC_DNA_16BIT::G;
        case 'R':
            return IUPAC_DNA_16BIT::R;
        case 'Y':
            return IUPAC_DNA_16BIT::Y;
        case 'S':
            return IUPAC_DNA_16BIT::S;
        case 'W':
            return IUPAC_DNA_16BIT::W;
        case 'K':
            return IUPAC_DNA_16BIT::K;
        case 'M':
            return IUPAC_DNA_16BIT::M;
        case 'B':
            return IUPAC_DNA_16BIT::B;
        case 'D':
            return IUPAC_DNA_16BIT::D;
        case 'H':
            return IUPAC_DNA_16BIT::H;
        case 'V':
            return IUPAC_DNA_16BIT::V;
        case 'N':
            return IUPAC_DNA_16BIT::N;
        case '-':
        case '?':
            return IUPAC_DNA_16BIT::GAP;
        default:
            throw runtime_error("Invalid nucleotide code: " + code);
    }
}

bool is_single_state(const IUPAC_DNA state) {
    assert(state != IUPAC_DNA::INVALID);
    auto s = static_cast<underlying_type<IUPAC_DNA>::type>(state);
    // s & (s - 1) will clear the least significant bit set
    return state != IUPAC_DNA::EMPTY && (s & (s - 1)) == 0;
}

bool is_single_state(const IUPAC_DNA_16BIT state) {
    assert(state != IUPAC_DNA_16BIT::INVALID);
    auto s = static_cast<underlying_type<IUPAC_DNA_16BIT>::type>(state);
    // s & (s - 1) will clear the least significant bit set
    return (s & (s - 1)) == 0;
}

bool is_single_state_16bit(const IUPAC_DNA_16BIT state) {
    return is_single_state(state);
}

bool valid(const IUPAC_DNA state) {
    return state != IUPAC_DNA::INVALID &&
           (static_cast<underlying_type<IUPAC_DNA>::type>(state) & 0b11110000) == 0;
}

bool valid(const IUPAC_DNA_16BIT state) {
    return state != IUPAC_DNA_16BIT::INVALID &&
           (static_cast<underlying_type<IUPAC_DNA_16BIT>::type>(state) & 0xFFFF0000) == 0;
}

CompressedSequence::CompressedSequence() :
    _sequence_length(0)
{}

CompressedSequence::CompressedSequence(const std::string& sequence) {
    from_string(sequence);
}

std::string CompressedSequence::to_string() {
    stringstream sstream;
    for (size_t idx = 0; idx < _sequence_length; idx++) {
        sstream << (*this)[idx];
    }

    return move(sstream.str());
}

IUPAC_DNA CompressedSequence::get_code(size_t idx) {
    if (idx >= _sequence_length) {
        throw runtime_error("Sequence index out of bounds.");
    }
    assert(ceil((position(idx) + STATE_ENCODING_LENGTH) / 8) <= compressed_size_in_bytes());
    
    return (IUPAC_DNA)_sequence.read_bits(position(idx), STATE_ENCODING_LENGTH);
}

char CompressedSequence::operator[](size_t idx) {
    return iupac_dna2char(get_code(idx));
}

void CompressedSequence::append(IUPAC_DNA code) {
    if (code == IUPAC_DNA::INVALID) {
        throw runtime_error("Trying to append invalid IUPAC code");
    }
    assert((static_cast<underlying_type<IUPAC_DNA>::type>(code) & 0b0000) == 0);

    _sequence.write_bits(code, STATE_ENCODING_LENGTH);
    _sequence_length++;
}

void CompressedSequence::append(CompressedSequence & sequence ) {
    for (size_t idx = 0; idx < sequence.length(); idx++) {
        this->append(sequence.get_code(idx));
    }
}

void CompressedSequence::append(char c) {
    IUPAC_DNA code = char2iupac_dna(c);
    append(code);
}

void CompressedSequence::from_string(const std::string& sequence) {
    if (_sequence_length > 0) {
        throw runtime_error("Sequence not empty");
    }

    for (char c: sequence) {
        append(c);
    }
    assert(_sequence_length == sequence.length());
    assert(compressed_size_in_bytes() == (sequence.length() + 1) / (8 / STATE_ENCODING_LENGTH));
}

const unsigned char CompressedSequence::STATE_ENCODING_LENGTH = IUPAC_DNA_CODE_LENGTH;
size_t CompressedSequence::length() const {
    return _sequence_length;
}

size_t CompressedSequence::compressed_size_in_bytes() const {
    return _sequence.size();
}

MSATreeCompression::MSATreeCompression(const string newick_string_rooted) :
   _tree(newick_string_rooted),
   _num_sequences(_tree.num_tips())
{
    if (_num_sequences < 2) {
        throw runtime_error("Can't compress a tree with less than two sequences.");
    }

    _encoding_of_changes = make_shared<bit_buffer>();
    _site_idx_in_encoding = make_shared<vector<size_t>>();
    assert(_encoding_of_changes);
    assert(_site_idx_in_encoding);
}

MSATreeCompression::MSATreeCompression(const Tree& tree) :
    // TODO This might be slow: benchmark/profile
    // We don't care where the root is as long as it's the same for compression and decompression
    MSATreeCompression(to_newick_string_rooted(tree))
{}

size_t MSATreeCompression::compress(const MSA& msa) {
    if (msa.size() != _num_sequences) {
        throw runtime_error("The number of sequences in the MSA and the number of leaves in the tree do not match.");
    }

    // Initialize variables
    _root_sequence = make_shared<CompressedSequence>();
    _encoding_of_changes = make_shared<bit_buffer>();

    // Compress tree msa column by msa column
    NodeStates node_states(_tree.num_nodes());
    size_t bits_written = 0;
    uint64_t changes_across_all_sites = 0;
    for (size_t position = 0; position < msa.length(); position++) {
        // Save current site's index in the output stream to enable random acces while decompressing
        _site_idx_in_encoding->push_back(bits_written);

        // Build ancestral states of inner nodes 
        uint64_t nChanges = build_ancestral_states(msa, position, node_states);
        changes_across_all_sites += nChanges;

        // Traverse along the tree, storing only the changes along the tree
        // and the root sequence
        _tree.traverse([this, &node_states, &bits_written, &nChanges](PllTree::Node * node) -> bool {
            IUPAC_DNA_16BIT my_state = node_states[node->node_id()];
            if (node->is_root()) {
                // Store character at root
                _root_sequence->append(iupac_dna_to_4bit(my_state));
            } else { // Inner node or leaf
                IUPAC_DNA_16BIT parent_state = node_states[node->parent()->node_id()];
                if (my_state != parent_state) {                        
                    // Store change
                    assert(node->node_id() <= _tree.num_nodes() - 1);
                    ChangeEncoder change_encoder(
                        node->node_id(),
                        _tree.num_nodes() - 1,
                        iupac_dna_to_4bit(parent_state),
                        iupac_dna_to_4bit(my_state)
                    );
                    _encoding_of_changes->write_bits(change_encoder.encoding(), change_encoder.bit_length());
                    bits_written += change_encoder.bit_length();
                    nChanges--;
                }
            }
            return true;
        }, PllTree::TRAVERSAL_MODE::PREORDER);

        assert(nChanges == 0);
        // Reset node states
        node_states.reset();
    }
    
    uint64_t changes_written = bits_written / ChangeEncoder::bit_length(_tree.num_nodes() - 1);
    assert(changes_across_all_sites == changes_written);

    // Guard element, like for adjacency arrays
    _site_idx_in_encoding->push_back(bits_written);

    //cout << "Encoded " << bits_written << " bits (" << bits_written / ChangeEncoder::bit_length(_tree.num_nodes() - 1) << " changes)" << endl;
    size_t size_of_encoding = ceil(bits_written / 8.); // in bytes
    size_of_encoding += _root_sequence->compressed_size_in_bytes();
    //cout << "the root sequence took up another " << _root_sequence->compressed_size_in_bytes() << " bytes" << endl;
    //cout << "The tree is stored using TODO bytes" << endl;
    // TODO
    return size_of_encoding;
}

shared_ptr<MSA> MSATreeCompression::decompress(RangeList& range_list) {
    // A vector to store the sequences. "Compressed" here means using 4 bits per state,
    // not compressed as in tree-based MSA compression.
    // We index this array by node_id of the corresponding tips in the phylogenetic tree.
    vector<CompressedSequence> sequences(_tree.num_nodes());
    fill(sequences.begin(), sequences.end(), CompressedSequence());
    sequences[_tree.root()->node_id()] = string();

    // Go over the list of ranges of sites and load each site from the encoding
    size_t num_decoded_sites = 0;
    for (auto range: range_list) {
        CompressedSequence segment_of_root_sequence =_root_sequence->slice(range.start, range.length);
        sequences[_tree.root()->node_id()].append(segment_of_root_sequence);

        for (size_t current_site = range.start; current_site < range.start + range.length; current_site++) {
            size_t max_node_id = _tree.num_nodes() - 1; // We use the node id as the edge id of the edge leading to this node in
                                                        // the directed and rooted phylogenetic tree used for compression.
            
            size_t position_in_encoding = (*_site_idx_in_encoding)[current_site];
            auto read_next_change = [this, max_node_id, &position_in_encoding, current_site]() -> unique_ptr<ChangeDecoder> {
                if ((*_site_idx_in_encoding)[current_site + 1] != position_in_encoding) { // There are more changes for this site
                    uint64_t encoding_of_change = _encoding_of_changes->read_bits(position_in_encoding, ChangeDecoder::bit_length(max_node_id));
                    position_in_encoding += ChangeDecoder::bit_length(max_node_id);

                    auto change = unique_ptr<ChangeDecoder>(new ChangeDecoder(encoding_of_change, max_node_id));
                    assert(change->edge_id() <= max_node_id);
                    return change;
                } else {
                    return nullptr; // No more changes for this site
                }
            };

            unique_ptr<ChangeDecoder> next_change = read_next_change();
            _tree.traverse([&sequences, &next_change, &num_decoded_sites, &current_site, &read_next_change](PllTree::Node * node) -> bool {
                if (node->is_root()) {
                    // We already set the root sequence
                    return true;
                } else { // Tip or inner node
                    // Assert, that the parent of this node has already been processed
                    assert(sequences[node->parent()->node_id()].length() >= num_decoded_sites + 1);

                    IUPAC_DNA this_nodes_state = sequences[node->parent()->node_id()].get_code(num_decoded_sites);
                    if (next_change && next_change->edge_id() == node->node_id()) { // see above
                        //cout << node->label() << ": at pos " << current_site << " changed " << this_nodes_state;
                        this_nodes_state = next_change->apply_change(this_nodes_state);
                        //cout << " to " << this_nodes_state << endl;
                        next_change = read_next_change();
                    } else {
                        //cout << node->label() << ": at pos " << current_site << " kept " << this_nodes_state << endl;
                    }

                    assert(sequences[node->node_id()].length() == num_decoded_sites);
                    sequences[node->node_id()].append(this_nodes_state);

                    return true;
                }
            }, PllTree::TRAVERSAL_MODE::PREORDER);

            num_decoded_sites++;
        }
    }

   auto msa = make_shared<MSA>();
   _tree.traverse([msa, &sequences](PllTree::Node * node) -> bool {
       if (node->is_leaf()) {
           msa->append(sequences[node->node_id()], node->label());
       }
       return true;
   }, PllTree::TRAVERSAL_MODE::INORDER);

   return msa;
}

CompressedSequence CompressedSequence::slice(size_t start, size_t length) {
    CompressedSequence seq;

    if (start + length > _sequence_length) {
        throw runtime_error("The requested slice is out of the sequence's bounds.");
    }

    for (size_t idx = start; idx < start + length; idx++) {
        seq.append((*this)[idx]);
    }

    return seq;
}

MSATreeCompression::NodeStates::NodeStates(size_t num_nodes) {
    if (num_nodes == 0) {
        throw runtime_error("NodeState objects with zero nodes are not allowed.");
    }
    _states.resize(num_nodes);
    fill(_states.begin(), _states.end(), IUPAC_DNA_16BIT::INVALID);
    assert(_states.size() == num_nodes);
}


IUPAC_DNA_16BIT MSATreeCompression::NodeStates::operator[](const size_t node_id) const {
    if(node_id >= _states.size()) {
        throw runtime_error("node_id out of range");
    }
    return _states[node_id];
}

IUPAC_DNA_16BIT& MSATreeCompression::NodeStates::operator[](const size_t node_id) {
    if(node_id >= _states.size()) {
        throw runtime_error("node_id out of range");
    }
    return _states[node_id];
}

IUPAC_DNA_16BIT MSATreeCompression::NodeStates::operator[](const PllTree::Node & node) const {
    return (*this)[node.node_id()];
}

IUPAC_DNA_16BIT& MSATreeCompression::NodeStates::operator[](const PllTree::Node & node) {
    return (*this)[node.node_id()];
}

IUPAC_DNA MSATreeCompression::NodeStates::choose_random_state(const IUPAC_DNA states) {
    if (states == IUPAC_DNA::INVALID) {
        throw runtime_error("Cannot choose random state form INVALID or EMPTY");
    }

    auto s = static_cast<underlying_type<IUPAC_DNA>::type>(states);
    
    // s & (s - 1) will clear the least significant bit set
    // ^-ing this with s will result in only the least significant bit set
    s = s ^ (s & (s - 1));

    return static_cast<IUPAC_DNA>(s);
}

IUPAC_DNA_16BIT MSATreeCompression::NodeStates::choose_random_state(const IUPAC_DNA_16BIT states) {
    if (states == IUPAC_DNA_16BIT::INVALID) {
        throw runtime_error("Cannot choose random state form INVALID or EMPTY");
    }

    auto s = static_cast<underlying_type<IUPAC_DNA_16BIT>::type>(states);
    
    // s & (s - 1) will clear the least significant bit set
    // ^-ing this with s will result in only the least significant bit set
    s = s ^ (s & (s - 1));

    return static_cast<IUPAC_DNA_16BIT>(s);
}

IUPAC_DNA_16BIT MSATreeCompression::NodeStates::choose_random_state(const PllTree::Node & node) {
    return ((*this)[node] = choose_random_state((*this)[node]));
}

void MSATreeCompression::NodeStates::reset() {
    fill(_states.begin(), _states.end(), IUPAC_DNA_16BIT::INVALID);
}

bool MSATreeCompression::NodeStates::check_all_single_state() const {
    return all_of(_states.begin(), _states.end(), is_single_state_16bit);
}

PllTree::PllTree(const pll_rtree_t * tree) :
    BasicTree(tree->tip_count),
    _tree(tree),
    _node_ptrs(num_nodes() + 1, nullptr)
{
    if (!tree) {
        throw runtime_error("Invalid tree pointer.");
    }
    if (tree->root == nullptr) {
        throw runtime_error("Tree root is not set.");
    }
    assert(num_nodes() > 0);
    assert(num_nodes() + 1 == _node_ptrs.size());
}

PllTree::PllTree(const string& newick_str) :
    PllTree(
        ([newick_str]() {
            if (newick_str.empty()) {
                throw runtime_error("Empty Newick string");
            }
            pll_rtree_t * tree =  pll_rtree_parse_newick_string(newick_str.c_str());
            libpll_check_error("ERROR reading Newick string");
            return tree;
        })()
    )
{ }

PllTree::~PllTree() {
    for (PllTree::Node * node_ptr: _node_ptrs) {
        if (node_ptr) {
            delete node_ptr;
            node_ptr = nullptr;
        }
    }
}

PllTree::Node * PllTree::root() {
    return create_Node_from_pll_node(_tree->root);
}

PllTree::Node::Node(const pll_rnode_t * pll_node, PllTree * tree) :
    _pll_node(pll_node),
    _tree(tree)
{
    if (!pll_node || !tree) {
        throw runtime_error("Either pll_node or tree or both is/are (a) nullptr(s)");
    }
}

bool operator==(const PllTree::Node & lhs, const PllTree::Node & rhs) {
    return lhs._pll_node == rhs._pll_node;
}

bool operator!=(const PllTree::Node & lhs, const PllTree::Node & rhs) {
    return !(lhs == rhs);
}

PllTree::Node * PllTree::create_Node_from_pll_node(const pll_rnode_t * pll_node) {
    size_t node_id = pll_node->node_index;
    assert(node_id < _node_ptrs.size());
    if (!_node_ptrs[node_id]) {
       _node_ptrs[node_id] = new Node(pll_node, this);
    }
    assert(_node_ptrs[node_id]);
    return _node_ptrs[node_id];
}

string PllTree::Node::label() const {
    assert(_pll_node);
    if(_pll_node->label == nullptr) {
        throw runtime_error("Label not set on node " + to_string(_pll_node->clv_index));
    }
    return move(string(_pll_node->label));
}

PllTree::Node * PllTree::Node::left_child() {
    assert(_pll_node);
    if(_pll_node->left) { 
        return _tree->create_Node_from_pll_node(_pll_node->left);
    }
    else {
        return nullptr;
    }
}

PllTree::Node * PllTree::Node::right_child() {
    assert(_pll_node);
    if(_pll_node->right) { 
        return _tree->create_Node_from_pll_node(_pll_node->right);
    }
    else {
        return nullptr;
    }
}

PllTree::Node * PllTree::Node::parent() {
    assert(_pll_node);
    if(_pll_node->parent) { 
        return _tree->create_Node_from_pll_node(_pll_node->parent);
    }
    else {
        return nullptr;
    }
}

bool PllTree::Node::is_leaf() const {
    assert(_pll_node);
    return pllmod_rtree_is_tip(_pll_node);
}

bool PllTree::Node::is_root() const {
    assert(_pll_node);
    return _pll_node == _tree->root()->_pll_node;
}

bool PllTree::Node::is_inner_node() const {
    return !is_leaf();
}

size_t PllTree::Node::node_id() const {
    assert(_pll_node);
    return _pll_node->node_index;
}

uint64_t MSATreeCompression::build_ancestral_states(const MSA& msa, size_t position, MSATreeCompression::NodeStates& node_states) {
    assert(node_states.size() == 2 * msa.size() - 1);
    assert(node_states.size() == _tree.num_nodes());
    
    // Phase 1: Build possible ancestral states
    // cout << "### Phase 1 starts" << endl;
    _tree.traverse([this, &node_states, &msa, position](PllTree::Node * node) -> bool {
        if (node->is_leaf()) {
            assert(node_states[*node] == IUPAC_DNA_16BIT::INVALID);
            node_states[*node] = char2iupac_dna_16bit(msa[node->label()][position]);
        } else {
            IUPAC_DNA_16BIT code_left = node_states[*(node->left_child())];
            IUPAC_DNA_16BIT code_right = node_states[*(node->right_child())];
            IUPAC_DNA_16BIT code_intersection = code_left & code_right;
            assert(valid(code_left));
            assert(valid(code_right));
            assert(valid(code_intersection));
            
            if (code_intersection != IUPAC_DNA_16BIT::EMPTY) {
                node_states[*node] = code_intersection;
            } else {
                IUPAC_DNA_16BIT code_union = code_left | code_right;
                assert(code_union != IUPAC_DNA_16BIT::EMPTY && valid(code_union));
                node_states[*node] = code_union;
            }

            // cout << "I'm node " << node->node_id() 
            //      << " left of me is " << code_left
            //      << " right of me is " << code_right
            //      << " (intersection: " << code_intersection
            //      << ") -> I choose " << node_states[*node]
            //      << endl;
            
        }
        return true;
    }, PllTree::TRAVERSAL_MODE::POSTORDER);

    // Phase 2: Select ancestral states
    // cout << "### Phase 2 starts" << endl;
    uint64_t nChanges = 0;
    _tree.traverse([this, &node_states, &nChanges](PllTree::Node * node) -> bool {
        if (node->is_root()) {
            node_states[*node] = node_states.choose_random_state(node_states[*node]);
            //cout << "I'm the root node, I have states " << node_states[*node]
            //     << " I choose state " << node_states[*node]<< endl;
            assert(is_single_state(node_states[*node]));
        } else {
            IUPAC_DNA_16BIT parent_state = node_states[*(node->parent())];
            assert(is_single_state(parent_state));

            // cout << "I'm node " << node->node_id() 
            //      << " my states are " << node_states[*node]  
            //      << " my parent's state is " << parent_state << " ";
            // As leaves have a single state, and the parent is also set to a single
            // state, this will only evaluate to true at a leaf node if our state and
            // our parent's state is equal -> leaves will always keep their state.
            if ((parent_state & node_states[*node]) != IUPAC_DNA_16BIT::EMPTY) {
                assert(!node->is_leaf() || node_states[*node] == parent_state);
                node_states[*node] = parent_state;
            } else {
                assert(node_states[*node] != IUPAC_DNA_16BIT::EMPTY);
                assert(valid(node_states[*node]));
                assert(!node->is_leaf() || node_states[*node] == node_states.choose_random_state(*node));
                node_states.choose_random_state(*node);
                nChanges++;
            }
            // cout << "-> I choose " << node_states[*node]
            //      << endl;
            assert(is_single_state(node_states[*node]));
        }
        return true;
    }, PllTree::TRAVERSAL_MODE::PREORDER);

    // _tree.traverse([&node_states](PllTree::Node * node) -> bool {
    //     if (node->is_leaf()) {
    //         cout << node->label() << ": " << node_states[*node] << endl;
    //     }
    //     return true;
    // }, PllTree::TRAVERSAL_MODE::INORDER);
    return nChanges;
}

void PllTree::traverse(function<bool(PllTree::Node *)> callback_on_node_visit, TRAVERSAL_MODE mode) {
    // We're doing a bit of the compilers work here by converting the (probably caputuring) lambda function converted
    // to a std::function object to a function pointer.
    // The idea is to pass the function object as the context to the pllmod call and the following function to
    // call the operator() of this function object as the callback for pllmod.
    int (*pllmod_callback) (pll_rnode_t *, void * context) = [](pll_rnode_t * pll_node, void * context) -> int {
        // Reinterpret the void * context as a function object
        auto context_pair = reinterpret_cast<pair<PllTree*,function<bool(PllTree::Node *)>> *>(context);
        PllTree * plltree_obj = context_pair->first;
        function<bool(PllTree::Node *)> func = context_pair->second;
        // Call the operator() of the function object
        return func(plltree_obj->create_Node_from_pll_node(pll_node));
    };

    // Choose traversal mode
    int (*pre_traversal) (pll_rnode_t *, void *) { nullptr };
    int (*in_traversal) (pll_rnode_t *, void *) { nullptr };
    int (*post_traversal) (pll_rnode_t *, void *) { nullptr };

    switch (mode) {
        case TRAVERSAL_MODE::PREORDER:
            pre_traversal = pllmod_callback;
            break;
        case TRAVERSAL_MODE::INORDER:
            in_traversal = pllmod_callback;
            break;
        case TRAVERSAL_MODE::POSTORDER:
            post_traversal = pllmod_callback;
    }

    // Use pll to perform the traversal
    std::pair<PllTree*,function<bool(PllTree::Node *)>> context_pair;
    context_pair.first = this;
    context_pair.second = callback_on_node_visit;
    auto ret = pllmod_rtree_traverse_apply(
        (pll_rnode_t*) _tree->root, // Yeah, I know, casting away constness is evil, but it's done inside libpll, too
        pre_traversal, in_traversal, post_traversal,
        reinterpret_cast<void*>(&context_pair) // The function object is passed as context
    );
    if (!ret) {
        throw runtime_error("Error traversing tree");
    }
    libpll_check_error("Error traversing tree");
}

CompressedMSA::CompressedMSA(const MSA& msa, const Tree& tree) :
    _tree_compression(tree)
{
    _tree_compression.compress(msa);
}

const MSA& CompressedMSA::decompress() { return *(new MSA()); };

uint8_t MSATreeCompression::ChangeEncoding::bit_length() const {
    return bits_for_edge_id() + IUPAC_DNA_CODE_LENGTH;
}

uint8_t MSATreeCompression::ChangeEncoding::bits_for_edge_id() const {
    return ceil(log2(_max_edge_id));
}

uint8_t MSATreeCompression::ChangeEncoding::bit_length(size_t max_edge_id) {
    return bits_for_edge_id(max_edge_id) + IUPAC_DNA_CODE_LENGTH;
}

uint8_t MSATreeCompression::ChangeEncoding::bits_for_edge_id(size_t max_edge_id) {
    return ceil(log2(max_edge_id));
}

MSATreeCompression::ChangeEncoder::ChangeEncoder(size_t edge_id, size_t max_edge_id, IUPAC_DNA from_state, IUPAC_DNA to_state) :
    ChangeEncoding(max_edge_id),
    _edge_id(edge_id),
    _from_state(from_state),
    _to_state(to_state)
{
    if (edge_id > max_edge_id) {
        throw runtime_error("The given edge id is greater than the maximum edge id.");
    }
    if (from_state == IUPAC_DNA::INVALID || to_state == IUPAC_DNA::INVALID) {
        throw runtime_error("Invalid node state");
    }
}

uint64_t MSATreeCompression::ChangeEncoder::encoding() const {
    uint64_t edge_id_shifted = _edge_id << 4;
    uint64_t nucleotide_change = static_cast<std::underlying_type<IUPAC_DNA>::type>(_from_state ^ _to_state);

    assert((edge_id_shifted & 0b0000) == 0);
    assert(nucleotide_change < 0b10000);
    
    return edge_id_shifted | nucleotide_change;
}

MSATreeCompression::ChangeDecoder::ChangeDecoder(size_t encoding, size_t max_edge_id) :
    ChangeEncoding(max_edge_id),
    _encoding(encoding)
{
    if (_encoding >> bit_length() != 0) {
        throw runtime_error("Bits which are not needed for the encoding are set. Is the maximum edge id correct?");
    }
    if (edge_id() > max_edge_id) {
        throw runtime_error("The edge id is greater than the maximum edge id.");
    }
}
    
IUPAC_DNA MSATreeCompression::ChangeDecoder::apply_change(const IUPAC_DNA state) const {
    uint64_t changemask = _encoding & 0b1111;
    assert(changemask < 0b10000);

    return static_cast<IUPAC_DNA>(
        static_cast<std::underlying_type<IUPAC_DNA>::type>(state) ^
        changemask
    );
};

size_t MSATreeCompression::ChangeDecoder::edge_id() const {
   return _encoding >> 4; 
}