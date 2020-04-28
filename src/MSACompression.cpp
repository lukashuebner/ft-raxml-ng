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
        case 'V':
            return IUPAC_DNA::V;
        case 'N':
            return IUPAC_DNA::N;
        case '-':
            return IUPAC_DNA::GAP;
        default:
            throw runtime_error("Invalid nucleotide code: " + code);
    }
}

bool is_single_state(const IUPAC_DNA state) {
    assert(state != IUPAC_DNA::INVALID);
    auto s = static_cast<underlying_type<IUPAC_DNA>::type>(state);
    // s & (s - 1) will clear the least significand bit set
    return state != IUPAC_DNA::EMPTY && (s & (s - 1)) == 0;
}

bool valid(const IUPAC_DNA state) {
    return state != IUPAC_DNA::INVALID &&
           (static_cast<underlying_type<IUPAC_DNA>::type>(state) & 0b11110000) == 0;
}

CompressedSequence::CompressedSequence(const std::string& sequence) {
    if (sequence.empty()) {
        throw runtime_error("An empty sequence is probably not what you want.");
    }
    from_string(sequence);
}

std::string CompressedSequence::to_string() {
    stringstream sstream;
    for (size_t idx = 0; idx < _sequence_length; idx++) {
        sstream << (*this)[idx];
    }

    return sstream.str();
}

char CompressedSequence::operator[](size_t idx) {
    if (idx >= _sequence_length) {
        throw runtime_error("Sequence index out of bounds.");
    }
    assert(position(idx) + STATE_ENCODING_LENGTH < compressed_size_in_bytes() * 8);
    return iupac_dna2char((IUPAC_DNA)_sequence.read_bits(position(idx), STATE_ENCODING_LENGTH));
}

void CompressedSequence::from_string(const std::string& sequence) {
    for (char c: sequence) {
        IUPAC_DNA code = char2iupac_dna(c);
        assert(code != IUPAC_DNA::INVALID);
        assert((static_cast<underlying_type<IUPAC_DNA>::type>(code) & 0b0000) == 0);
        _sequence.write_bits(code, STATE_ENCODING_LENGTH);
        _sequence_length++;
    }
    assert(_sequence_length == sequence.length());
    assert(compressed_size_in_bytes() == (sequence.length() + 1) / (8 / STATE_ENCODING_LENGTH));
}

const unsigned char CompressedSequence::STATE_ENCODING_LENGTH = 4;
size_t CompressedSequence::length() const {
    return _sequence_length;
}

size_t CompressedSequence::compressed_size_in_bytes() const {
    return _sequence.size();
}

MSATreeCompression::MSATreeCompression(const Tree& tree) :
    num_sequences(tree.num_tips()),
    _tree(&(tree.pll_utree()))
{
    if (num_sequences < 2) {
        throw runtime_error("Can't compress a tree with less than two sequences.");
    }
    _change_encoding = make_shared<bit_buffer>();
    assert(_change_encoding);
}

void MSATreeCompression::compress(const MSA& msa) {
    if (msa.num_patterns() != num_sequences) {
        throw runtime_error("The number of sequences in the MSA and the number of leaves in the tree do not match.");
    }

    // Use new objects in case anyone has still references to it.
    // Repeated compression is probably not that much of an use case.
    if (working_buffers_allocated()) {
        free_working_buffers();
    }

    // Allocate memory for sequences of inner nodes and leaves
    _node_sequences = make_shared<vector<shared_ptr<CompressedSequence>>>();
    assert(_node_sequences);
    _node_sequences->reserve(msa.num_patterns());
    
    // Compress and store leaf sequences
    for (auto & sequence: msa) {
        _node_sequences->push_back(make_shared<CompressedSequence>(sequence));
    }
    assert(_node_sequences->size() == msa.num_patterns());

    _root_sequence = _node_sequences->at(_root_sequence_idx);
    assert(_root_sequence);

    // Compress tree msa column by msa column
    NodeStates node_states(_tree.num_nodes());
    for (size_t position = 0; position < msa.length(); position++) {
        // Build ancestral states of inner nodes 
        build_ancestral_states(msa, position, node_states);

        // Traverse along the tree, storing only the changes

    }
}

MSATreeCompression::NodeStates::NodeStates(size_t num_nodes) {
    _states.reserve(num_nodes);
    #ifndef NDEBUG
    fill(_states.begin(), _states.end(), IUPAC_DNA::INVALID);
    #endif 
}

IUPAC_DNA MSATreeCompression::NodeStates::operator[](PllTree::Node node) const {
    size_t node_id = PllTree::node_id(node);
    assert(node_id < _states.size());
    return _states[node_id];
}

IUPAC_DNA& MSATreeCompression::NodeStates::operator[](PllTree::Node node) {
    size_t node_id = PllTree::node_id(node);
    assert(node_id < _states.size());
    return _states[node_id];
}

IUPAC_DNA MSATreeCompression::NodeStates::choose_random_state(const IUPAC_DNA states) {
    auto s = static_cast<underlying_type<IUPAC_DNA>::type>(states);
    
    // s & (s - 1) will clear the least significant bit set
    // ^-ing this with s will result in only the least significant bit set
    s = s ^ (s & (s - 1));

    return static_cast<IUPAC_DNA>(s);
}

IUPAC_DNA MSATreeCompression::NodeStates::choose_random_state(PllTree::Node node) {
    return ((*this)[node] = choose_random_state((*this)[node]));
}

void MSATreeCompression::NodeStates::reset() {
    #ifndef NDEBUG
    fill(_states.begin(), _states.end(), IUPAC_DNA::INVALID);
    #endif 
}

bool MSATreeCompression::NodeStates::check_all_single_state() {
    return all_of(_states.begin(), _states.end(), is_single_state);
}

PllTree::PllTree(const pll_utree_t * tree) :
    BasicTree(tree->tip_count),
    _tree(pll_utree_wraptree(pll_utree_graph_clone(tree->vroot), tree->tip_count)),
    _tree_root(_tree->vroot)
{
    if (!tree) {
        throw runtime_error("Invalid tree pointer.");
    }
    if (!_tree_root) {
        throw runtime_error("(Virtual) tree root is not set.");
    }
}

PllTree::Node PllTree::root() const {
    return _tree_root;
}

string PllTree::label(Node node) {
    assert(node);
    if(node->label == nullptr) {
        throw runtime_error("Label not set on node " + to_string(node->clv_index));
    }
    return string(node->label);
}

PllTree PllTree::load_from_file(const string& file_name) {
    Tree tree;
    NewickStream ns(file_name, std::ios::in);
    if (!ns.good()) {
        throw runtime_error("Cannot open " + file_name);
    }

    ns >> tree;

    return PllTree(&(tree.pll_utree()));
}

PllTree::Node PllTree::left_child(PllTree::Node node) {
    assert(node && node->next && node->next->next);
    assert(node->next->next->next == node); // bifurcating tree
    return node->next->next;
}

PllTree::Node PllTree::right_child(PllTree::Node node) {
    assert(node && node->next && node->next->next);
    assert(node->next->next->next == node); // bifurcating tree
    return node->next;
}

PllTree::Node PllTree::parent(PllTree::Node node) {
    assert(node && node->back);
    return node->back;
}

bool PllTree::is_leaf(PllTree::Node node) {
    assert(node);
    return pllmod_utree_is_tip(node);
}

bool PllTree::is_root(PllTree::Node node) const {
    assert(node);
    return node == _tree_root;
}

bool PllTree::is_inner_node(PllTree::Node node) {
    assert(node);
    return !is_leaf(node);
}

size_t PllTree::node_id(PllTree::Node node) {
    assert(node);
    return node->clv_index;
}

void MSATreeCompression::build_ancestral_states(const MSA& msa, size_t position, MSATreeCompression::NodeStates& node_states) {
    assert(node_states.size() == msa.num_patterns());
    assert(msa.num_patterns() == _tree.num_nodes());
    
    // Phase 1: Build possible ancestral states
    _tree.traverse([this, &node_states, &msa, position](pll_unode_t * node) -> bool {
        if (_tree.is_leaf(node)) {
            assert(node_states[node] == IUPAC_DNA::INVALID);
            node_states[node] = char2iupac_dna(msa[PllTree::node_id(node)][position]);
        } else {
            IUPAC_DNA code_left = node_states[_tree.left_child(node)];
            IUPAC_DNA code_right = node_states[_tree.right_child(node)];
            IUPAC_DNA code_intersection = code_left | code_right;
            assert(code_left != IUPAC_DNA::INVALID);
            assert(code_right != IUPAC_DNA::INVALID);
            assert(code_intersection != IUPAC_DNA::INVALID);
            
            if (code_intersection != IUPAC_DNA::EMPTY) {
                node_states[node] = code_intersection;
            } else {
                IUPAC_DNA code_union = code_left & code_right;
                assert(code_union != IUPAC_DNA::EMPTY && code_union != IUPAC_DNA::INVALID);
                node_states[node] = code_union;
            }
        }
        return true;
    }, PllTree::TRAVERSAL_MODE::PREORDER);

    // Phase 2: Select ancestral states
    _tree.traverse([this, &node_states](pll_unode_t * node) -> bool {
        if (_tree.is_root(node)) {
            node_states.choose_random_state(node_states[node]);
        } else {
            IUPAC_DNA parent_state = node_states[_tree.parent(node)];
            assert(is_single_state(parent_state));

            if ((parent_state & node_states[node]) != IUPAC_DNA::EMPTY) {
                node_states[node] = parent_state;
            } else {
                node_states.choose_random_state(node);
            }
        }
        return true;
    }, PllTree::TRAVERSAL_MODE::POSTORDER);
}

void PllTree::traverse(function<bool(PllTree::Node)> callback_on_node_visit, TRAVERSAL_MODE mode) {
    // We're doing a bit of the compilers work here by converting the (probably caputuring) lambda function converted
    // to a std::function object to a function pointer.
    // The idea is to pass the function object as the context to the pllmod call and the following function to
    // call the operator() of this function object as the callback for pllmod.
    int (*pllmod_callback) (pll_unode_t *, void * context) = [](pll_unode_t * node, void * context) -> int {
        // Reinterpret the void * context as a function object
        function<bool(pll_unode_t *)> * func = reinterpret_cast<function<bool(pll_unode_t *)> *>(context);
        // Call the operator() of the function object
        return (*func)(node);
    };

    // Choose traversal mode
    int (*pre_traversal) (pll_unode_t *, void *) { nullptr };
    int (*in_traversal) (pll_unode_t *, void *) { nullptr };
    int (*post_traversal) (pll_unode_t *, void *) { nullptr };

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
    assert(_tree->vroot == _tree_root);
    auto ret = pllmod_utree_traverse_apply(
        (pll_unode_t*) _tree_root, // Yeah, I know, casting away constness is evil, but it's done inside libpll, too
        pre_traversal, in_traversal, post_traversal,
        reinterpret_cast<void*>(&callback_on_node_visit) // The function object is passed as context
    );
    assert(ret);
}

CompressedMSA::CompressedMSA(const MSA& msa, const Tree& tree) :
    _tree_compression(tree)
{
    _tree_compression.compress(msa);
    _tree_compression.free_working_buffers();
}

const MSA& CompressedMSA::decompress() { return *(new MSA()); };

void MSATreeCompression::free_working_buffers() {
    _node_sequences = nullptr;
}

bool MSATreeCompression::working_buffers_allocated() {
    return _node_sequences != nullptr;
}
