#include <string>
#include <iostream>
#include <unistd.h>

#include "RaxmlTest.hpp"
#include "../../src/MSACompression.hpp"

using namespace std;

TEST(IUPAC_DNATest, Operators) {
    ASSERT_EQ(IUPAC_DNA::A | IUPAC_DNA::C | IUPAC_DNA::G, IUPAC_DNA::V);
    ASSERT_EQ(IUPAC_DNA::V | IUPAC_DNA::T, IUPAC_DNA::N);
    ASSERT_EQ(IUPAC_DNA::R & IUPAC_DNA::G, IUPAC_DNA::G);
    ASSERT_EQ(IUPAC_DNA::Y & IUPAC_DNA::S, IUPAC_DNA::C);
    ASSERT_EQ(IUPAC_DNA::S & IUPAC_DNA::W, IUPAC_DNA::EMPTY);
}

TEST(IUPAC_DNATest, Helpers) {
    const vector<IUPAC_DNA> single_types = { IUPAC_DNA::A, IUPAC_DNA::C, IUPAC_DNA::T, IUPAC_DNA::G };
    const vector<IUPAC_DNA> multi_types =  { IUPAC_DNA::B, IUPAC_DNA::D, IUPAC_DNA::K, IUPAC_DNA::M, IUPAC_DNA::N,
                                             IUPAC_DNA::R, IUPAC_DNA::S, IUPAC_DNA::V, IUPAC_DNA::W, IUPAC_DNA::Y };
    const vector<IUPAC_DNA> special_types = { IUPAC_DNA::GAP, IUPAC_DNA::EMPTY };
    vector<IUPAC_DNA> special_multi_types(special_types);
    special_multi_types.insert(special_multi_types.end(), multi_types.begin(), multi_types.end());
    vector<IUPAC_DNA> single_multi_types(single_types);
    single_multi_types.insert(single_multi_types.end(), multi_types.begin(), multi_types.end());
    vector<IUPAC_DNA> all_types(special_multi_types);
    all_types.insert(all_types.end(), single_types.begin(), single_types.end());

    // is_single_state
    for (auto state: single_types) {
        ASSERT_TRUE(is_single_state(state));
    }

    for (auto state: special_multi_types) {
        ASSERT_FALSE(is_single_state(state));
    }

    // char <-> code conversion
    ASSERT_ANY_THROW(char2iupac_dna('x'));
    ASSERT_ANY_THROW(char2iupac_dna('\n'));
    ASSERT_ANY_THROW(iupac_dna2char(IUPAC_DNA::INVALID));

    for (auto state: single_multi_types) {
        ASSERT_EQ(char2iupac_dna(iupac_dna2char(state)), state);
    }
    ASSERT_EQ(char2iupac_dna('-'), IUPAC_DNA::GAP);
    ASSERT_EQ(char2iupac_dna('A'), IUPAC_DNA::A);
    ASSERT_EQ(char2iupac_dna('R'), IUPAC_DNA::R);
    ASSERT_EQ(char2iupac_dna('N'), IUPAC_DNA::N);
    ASSERT_EQ(iupac_dna2char(IUPAC_DNA::GAP), '-');
    ASSERT_EQ(iupac_dna2char(IUPAC_DNA::W), 'W');
    ASSERT_EQ(iupac_dna2char(IUPAC_DNA::S), 'S');
    ASSERT_EQ(iupac_dna2char(IUPAC_DNA::M), 'M');

    // validity
    for (auto state: all_types) {
        ASSERT_TRUE(valid(state));
    }
    ASSERT_FALSE(valid(IUPAC_DNA::INVALID));
    ASSERT_FALSE(valid(static_cast<IUPAC_DNA>(0b0100000)));
    ASSERT_FALSE(valid(static_cast<IUPAC_DNA>(-1)));
}

TEST(CompressedSequenceTest, EncodeDecode) {
    string char_seq("ACTGRYSWKMBDVN-");
    auto compressed_seq = CompressedSequence(char_seq);
    
    ASSERT_EQ(compressed_seq.to_string(), char_seq);
    ASSERT_EQ(compressed_seq.length(), char_seq.length());
    ASSERT_EQ(compressed_seq.compressed_size_in_bytes(), (compressed_seq.length() + 1) / 2);
    ASSERT_EQ(CompressedSequence("ACT").compressed_size_in_bytes(), 2);
    ASSERT_EQ(CompressedSequence("ACV-").compressed_size_in_bytes(), 2);
    
    for (size_t idx = 0; idx < char_seq.length(); idx++) {
        ASSERT_EQ(compressed_seq[idx], char_seq[idx]);
    }
}

TEST(CompressedSequenceTest, GarbageIn) {
    ASSERT_ANY_THROW(CompressedSequence("XAC"));
    ASSERT_ANY_THROW(CompressedSequence(""));
    ASSERT_ANY_THROW(CompressedSequence("XAC"));

    CompressedSequence sequence("ACTG");
    ASSERT_ANY_THROW(sequence[4]);
    ASSERT_ANY_THROW(sequence[-1]);
}

class PllTreeTest : public ::testing::Test {
    protected:
    shared_ptr<PllTree> simple_tree;
    shared_ptr<PllTree> another_tree;
    shared_ptr<PllTree> elaborate_tree;

    virtual void SetUp() {
        cout << string(get_current_dir_name()) << endl;
        simple_tree = make_shared<PllTree>(PllTree::load_from_file("test/data/simple.newick"));
        // another_tree = make_shared<PllTree>(PllTree::load_from_file("test/data/another.newick"));
        // elaborate_tree = make_shared<PllTree>(PllTree::load_from_file("test/data/elaborate.newick"));
    }
};

TEST_F(PllTreeTest, SimpleTree) {
    // root, left_child, right_child, parent
    PllTree::Node root, left_child, right_child;
    ASSERT_NO_THROW(root = simple_tree->root());
    ASSERT_NO_THROW(left_child = simple_tree->left_child(root));
    ASSERT_NO_THROW(right_child = simple_tree->right_child(root));
    ASSERT_NE(root, nullptr);
    ASSERT_NE(left_child, nullptr);
    ASSERT_NE(right_child, nullptr);

    ASSERT_EQ(simple_tree->label(root), string("root"));
    ASSERT_EQ(simple_tree->label(left_child), string("left_child"));
    ASSERT_EQ(simple_tree->label(right_child), string("right_child"));

    ASSERT_EQ(simple_tree->parent(left_child), root);
    ASSERT_EQ(simple_tree->parent(right_child), root);

    // is_leaf, is_root, is_inner_node
    ASSERT_TRUE(simple_tree->is_root(root));
    ASSERT_TRUE(simple_tree->is_inner_node(root));
    ASSERT_FALSE(simple_tree->is_leaf(root));

    ASSERT_FALSE(simple_tree->is_root(left_child));
    ASSERT_FALSE(simple_tree->is_inner_node(left_child));
    ASSERT_TRUE(simple_tree->is_leaf(left_child));

    ASSERT_FALSE(simple_tree->is_inner_node(right_child));
    ASSERT_TRUE(simple_tree->is_leaf(right_child));
    ASSERT_TRUE(simple_tree->is_leaf(right_child));

    // node_id in range
}

TEST_F(PllTreeTest, TreeMSAInterplay) {
    // Are all node_ids of leaves valid MSA indices?

    // Does the label of all leaves match between MSA and node structure?
}

TEST_F(PllTreeTest, GarbageIn) {
    PllTree::Node root, left_child, right_child;
    root = simple_tree->root();
    left_child = simple_tree->left_child(root);
    right_child = simple_tree->right_child(root);
    assert(root && left_child && right_child);

    // Loading too small tree
    ASSERT_ANY_THROW(PllTree::load_from_file("test/data/single_taxon.newick"));

    // Getting label on node with no label set
    auto no_label_tree = PllTree::load_from_file("test/data/no_labels.newick");
    ASSERT_ANY_THROW(no_label_tree.label(no_label_tree.left_child(no_label_tree.root())));
    ASSERT_ANY_THROW(no_label_tree.label(no_label_tree.right_child(no_label_tree.root())));
    ASSERT_ANY_THROW(no_label_tree.label(no_label_tree.root()));

    // parent on root
    ASSERT_EQ(simple_tree->parent(root), nullptr);

    // left_child and right_child on leaves
    ASSERT_EQ(simple_tree->right_child(right_child), nullptr);
    ASSERT_EQ(simple_tree->left_child(right_child), nullptr);
    ASSERT_EQ(simple_tree->right_child(left_child), nullptr);
    ASSERT_EQ(simple_tree->left_child(left_child), nullptr);
}

TEST_F(PllTreeTest, Traversal) {
    // Check if the order in which the nodes are visited are correct

}

TEST(NodeStatesTest, Basics) {
    // Getting and setting using operator[]

    // size

    // choosing random states

    // check if all states are single states
}

TEST(NodeStateTest, InvalidOperations) {
    // Setting invalid states using operator[]

    // Choosing random state from invalid state
}

TEST(MSATreeCompressionTest, Basics) {
    // Is num_sequences set correctly?

    // Are the working buffers freed upon request?
}

TEST(MSATreeCompressionTest, GarbageIn) {
    // Invalid tree

    // Invalid MSA
}

TEST(MSATreeCompressionTest, CompressionAndDecompression) {
    // Is compression the inverse of decompression?

    // Is the compression' size the expected one?
}

TEST(MSATreeCompressionTest, ComputingAncestralStates) {
    // Are the ancestral states choosen correctly?
}

TEST(CompressedMSATest, Basics) {
    // Is compression the inverse of decompression?
}

TEST(CompressedMSATest, GarbageIn) {
    // Invalid tree

    // Invalid MSA
}