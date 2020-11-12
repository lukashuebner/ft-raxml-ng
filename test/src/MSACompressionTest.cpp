#include <string>
#include <iostream>
#include <unistd.h>

#include "RaxmlTest.hpp"
#include "../../src/MSACompression.hpp"

using namespace std;
using namespace testing;

// Any function which uses ASSERT_* must have a void return type. Any function we want to pass to
// PllTree.traverse(...) must have a boolean return type.
// This decorator accepts a function returning void and returns a function doing exactly the same
// but returning a boolean.
function<bool(PllTree::Node * node)> decorate_for_success(function<void(PllTree::Node * node)> func) {
    return [&func](PllTree::Node * node) {
        func(node);
        return true;
    };
}

TEST(IUPAC_DNATest, Operators) {
    // IUPAC_DNA
    ASSERT_EQ(IUPAC_DNA::A | IUPAC_DNA::C | IUPAC_DNA::G, IUPAC_DNA::V);
    ASSERT_EQ(IUPAC_DNA::V | IUPAC_DNA::T, IUPAC_DNA::N);
    ASSERT_EQ(IUPAC_DNA::R & IUPAC_DNA::G, IUPAC_DNA::G);
    ASSERT_EQ(IUPAC_DNA::Y & IUPAC_DNA::S, IUPAC_DNA::C);
    ASSERT_EQ(IUPAC_DNA::S & IUPAC_DNA::W, IUPAC_DNA::EMPTY);

    // IUPAC_DNA_16BIT
    // TODO
    // ASSERT_EQ(IUPAC_DNA_16BIT::A | IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::G, 0b1101000000000000);
    // ASSERT_EQ(IUPAC_DNA_16BIT::V | IUPAC_DNA_16BIT::T, IUPAC_DNA::N);
    // ASSERT_EQ(IUPAC_DNA_16BIT::R & IUPAC_DNA_16BIT::G, IUPAC_DNA::G);
    // ASSERT_EQ(IUPAC_DNA_16BIT::Y & IUPAC_DNA_16BIT::S, IUPAC_DNA::C);
    // ASSERT_EQ(IUPAC_DNA_16BIT::S & IUPAC_DNA_16BIT::W, IUPAC_DNA::EMPTY);
}

TEST(IUPAC_DNATest, Helpers) {
    // TODO Test IUPAC_DNA_16BIT
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

    auto all_types_16bit = {
        IUPAC_DNA_16BIT::A, IUPAC_DNA_16BIT::C, IUPAC_DNA_16BIT::G, IUPAC_DNA_16BIT::T,
        IUPAC_DNA_16BIT::R, IUPAC_DNA_16BIT::Y, IUPAC_DNA_16BIT::S, IUPAC_DNA_16BIT::W,
        IUPAC_DNA_16BIT::K, IUPAC_DNA_16BIT::M, IUPAC_DNA_16BIT::B, IUPAC_DNA_16BIT::D,
        IUPAC_DNA_16BIT::H, IUPAC_DNA_16BIT::V, IUPAC_DNA_16BIT::N, IUPAC_DNA_16BIT::GAP        
    };

    // is_single_state
    for (auto state: single_types) {
        ASSERT_TRUE(is_single_state(state));
    }

    for (auto state: special_multi_types) {
        ASSERT_FALSE(is_single_state(state));
    }

    // char <-> code conversion - IUPAC_DNA
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

    // char <-> code conversion - IUPAC_DNA_16BIT
    ASSERT_ANY_THROW(char2iupac_dna_16bit('x'));
    ASSERT_ANY_THROW(char2iupac_dna_16bit('\n'));
    
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::INVALID), "(Invalid)");

    for (auto state: all_types_16bit) {
        ASSERT_EQ(char2iupac_dna_16bit(iupac_dna2string_16bit(state)[1]), state);
    }
    ASSERT_EQ(char2iupac_dna_16bit('-'), IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(char2iupac_dna_16bit('A'), IUPAC_DNA_16BIT::A);
    ASSERT_EQ(char2iupac_dna_16bit('R'), IUPAC_DNA_16BIT::R);
    ASSERT_EQ(char2iupac_dna_16bit('N'), IUPAC_DNA_16BIT::N);
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::GAP), "(-)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::W), "(W)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::S), "(S)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::M), "(M)");

    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::W | IUPAC_DNA_16BIT::A), "(W,A)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::V), "(V,C)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::N | IUPAC_DNA_16BIT::K), "(N,K)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::GAP | IUPAC_DNA_16BIT::A), "(-,A)");
    ASSERT_EQ(iupac_dna2string_16bit(IUPAC_DNA_16BIT::GAP | IUPAC_DNA_16BIT::GAP), "(-)");
    
    ASSERT_EQ(
        iupac_dna2string_16bit(
            IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::V | IUPAC_DNA_16BIT::GAP
        ),
        "(-,V,C)"
    );

    ASSERT_EQ(
        iupac_dna2string_16bit(
            IUPAC_DNA_16BIT::K | IUPAC_DNA_16BIT::T | IUPAC_DNA_16BIT::H | IUPAC_DNA_16BIT::R
        ),
        "(H,K,R,T)"
    );

    ASSERT_EQ(
        iupac_dna2string_16bit(
            IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::V | IUPAC_DNA_16BIT::G | IUPAC_DNA_16BIT::Y | IUPAC_DNA_16BIT::B
        ),
        "(V,B,Y,G,C)"
    );

    // validity - IUPAC_DNA
    for (auto state: all_types) {
        ASSERT_TRUE(valid(state));
    }
    ASSERT_FALSE(valid(IUPAC_DNA::INVALID));
    ASSERT_FALSE(valid(static_cast<IUPAC_DNA>(0b0100000)));
    ASSERT_FALSE(valid(static_cast<IUPAC_DNA>(-1)));

    // validity - IUPAC_DNA_16BIT
    for (auto state: all_types_16bit) {
        ASSERT_TRUE(valid(state));
    }
    ASSERT_FALSE(valid(IUPAC_DNA_16BIT::INVALID));
    ASSERT_FALSE(valid(static_cast<IUPAC_DNA_16BIT>(0xFF0F00FF)));
    ASSERT_FALSE(valid(static_cast<IUPAC_DNA_16BIT>(-1)));
}

TEST(CompressedSequenceTest, EncodeDecode) {
    string char_seq("ACTGRYSWKMBDVN-");
    auto compressed_seq = CompressedSequence(char_seq);
    
    ASSERT_EQ(compressed_seq.to_string(), char_seq);
    ASSERT_EQ(compressed_seq.length(), char_seq.length());
    ASSERT_EQ(compressed_seq.compressed_size_in_bytes(), (compressed_seq.length() + 1) / 2);

    ASSERT_EQ(CompressedSequence("ACT").compressed_size_in_bytes(), 2);
    ASSERT_EQ(CompressedSequence("ACV-").compressed_size_in_bytes(), 2);
    ASSERT_EQ(CompressedSequence("G").compressed_size_in_bytes(), 1);
    ASSERT_EQ(CompressedSequence("MN-TAC").compressed_size_in_bytes(), 3);
    ASSERT_EQ(CompressedSequence("---TG-CCA").compressed_size_in_bytes(), 5);
    
    for (size_t idx = 0; idx < char_seq.length(); idx++) {
        ASSERT_EQ(compressed_seq[idx], char_seq[idx]);
    }
}

TEST(CompressedSequenceTest, ToString) {
    auto compressed_sequence = CompressedSequence("ACTGRYSWKMBDVN-");
    ASSERT_EQ(compressed_sequence.to_string(), "ACTGRYSWKMBDVN-");
}

TEST(CompressedSequenceTest, Append) {
    // Building a compressed sequence char-by-char (used to build the root sequence)
    CompressedSequence root_seq;
    root_seq.append('A');
    root_seq.append('C');
    root_seq.append(IUPAC_DNA::T);
    root_seq.append(IUPAC_DNA::G);

    ASSERT_EQ(root_seq.to_string(), "ACTG");
    ASSERT_EQ(root_seq.length(), 4);
    ASSERT_EQ(root_seq.compressed_size_in_bytes(), 2);

    CompressedSequence s1("AC"), s2("GT");
    s1.append(s2);
    ASSERT_EQ(s1.to_string(), "ACGT");
    ASSERT_EQ(s2.to_string(), "GT");
}

TEST(CompressedSequenceTest, Slice) {
    CompressedSequence seq("ACTGRYSWKMBDVN-");  
    ASSERT_EQ(seq.slice(0, 1).to_string(), "A");
    ASSERT_EQ(seq.slice(0, 3).to_string(), "ACT");
    ASSERT_EQ(seq.slice(2, 3).to_string(), "TGR");
    ASSERT_EQ(seq.slice(0, 15).to_string(), "ACTGRYSWKMBDVN-");
}

TEST(CompressedSequenceTest, GarbageIn) {
    // Invalid characters
    ASSERT_ANY_THROW(CompressedSequence("XAC"));
    
    // Out-of-bounds access
    CompressedSequence sequence("ACTG");
    ASSERT_ANY_THROW(sequence[4]);
    ASSERT_ANY_THROW(sequence[-1]);
    ASSERT_ANY_THROW(sequence.slice(5, 1));
    ASSERT_ANY_THROW(sequence.slice(2, 10));
}

class PllTreeTest : public ::testing::Test {
    protected:
    shared_ptr<PllTree> simple_tree;
    shared_ptr<PllTree> another_tree;
    shared_ptr<PllTree> elaborate_tree;

    virtual void SetUp() {
        simple_tree = make_shared<PllTree>("((bat,cow)innerA,(elk,fox)innerB)root;");
        another_tree = make_shared<PllTree>("(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;");
        elaborate_tree = make_shared<PllTree>("((ll,(lrl,lrr)lr)l,(rl,((rrll,rrlr)rrl,(rrrl,rrrr)rrr)rr)r)root;");
    }
};

TEST_F(PllTreeTest, SimpleTree) {
    // root, left_child, right_child, parent
    PllTree::Node * root, * innerA, * innerB, * bat, * cow;
    ASSERT_NO_THROW(root = simple_tree->root());
    ASSERT_NO_THROW(innerA = root->left_child());
    ASSERT_NO_THROW(innerB = root->right_child());
    ASSERT_NO_THROW(bat = innerA->left_child());
    ASSERT_NO_THROW(cow = innerA->right_child());

    ASSERT_NE(root, nullptr);
    ASSERT_NE(innerA, nullptr);
    ASSERT_NE(innerB, nullptr);
    ASSERT_NE(bat, nullptr);
    ASSERT_NE(cow, nullptr);

    ASSERT_EQ(root->label(), string("root"));
    ASSERT_EQ(innerA->label(), string("innerA"));
    ASSERT_EQ(innerB->label(), string("innerB"));
    ASSERT_EQ(bat->label(), string("bat"));
    ASSERT_EQ(cow->label(), string("cow"));

    ASSERT_EQ(innerA->parent(), root);
    ASSERT_EQ(innerB->parent(), root);
    ASSERT_EQ(bat->parent(), innerA);
    ASSERT_EQ(bat->parent(), innerA);

    // is_leaf, is_root, is_inner_node
    ASSERT_TRUE(root->is_root());
    ASSERT_TRUE(root->is_inner_node());
    ASSERT_FALSE(root->is_leaf());

    ASSERT_FALSE(innerA->is_root());
    ASSERT_TRUE(innerA->is_inner_node());
    ASSERT_FALSE(innerA->is_leaf());

    ASSERT_FALSE(innerB->is_root());
    ASSERT_TRUE(innerB->is_inner_node());
    ASSERT_FALSE(innerB->is_leaf());

    ASSERT_FALSE(bat->is_root());
    ASSERT_FALSE(bat->is_inner_node());
    ASSERT_TRUE(bat->is_leaf());

    // node_id in range
    for (auto node: {root, innerA, innerB, cow, bat}) {
        size_t node_id = node->node_id();
        ASSERT_LE(node_id, 6);
    }
}

TEST_F(PllTreeTest, AnotherTree) {
    PllTree::Node * one, * two, * three, * four;
    
    one = another_tree->root()->left_child()->left_child()->left_child();
    two = one->parent()->right_child();
    three = two->parent()->parent()->right_child()->left_child();
    four = three->parent()->right_child();

    // Is each node the one we're trying to get?
    ASSERT_EQ(one->label(), "One");
    ASSERT_EQ(two->label(), "Two");
    ASSERT_EQ(three->label(), "Three");
    ASSERT_EQ(four->label(), "Four");
    
    // Do the parents match up?
    ASSERT_EQ(one->parent(), two->parent());
    ASSERT_EQ(three->parent(), four->parent());
    ASSERT_EQ(one->parent()->parent(), four->parent()->parent());

    // Some leaf and inner node tests
    ASSERT_TRUE(another_tree->root()->is_root());
    ASSERT_TRUE(another_tree->root()->is_inner_node());
    ASSERT_FALSE(another_tree->root()->is_leaf());

    ASSERT_FALSE(one->parent()->is_root());
    ASSERT_TRUE(one->parent()->is_inner_node());
    ASSERT_FALSE(one->parent()->is_leaf());

    ASSERT_FALSE(one->is_root());
    ASSERT_FALSE(one->is_inner_node());
    ASSERT_TRUE(one->is_leaf());

    // Some tree properties
    ASSERT_TRUE(another_tree->binary());
    ASSERT_FALSE(another_tree->empty());
    ASSERT_EQ(another_tree->num_inner(), 4);
    ASSERT_EQ(another_tree->num_branches(), 8);
    ASSERT_EQ(another_tree->num_tips(), 5);
    ASSERT_EQ(another_tree->num_nodes(), 9);
}

TEST_F(PllTreeTest, ElaborateTree) {
    PllTree::Node *root, *l, *ll, *lr, *lrl, *lrr, *r, *rl, *rr, *rrl, *rrll, *rrlr, *rrr, *rrrl, *rrrr;

    root = elaborate_tree->root();
    l = root->left_child();
    r = root->right_child();
    ll = l->left_child();
    lr = l->right_child();
    rl = r->left_child();
    rr = r->right_child();
    lrl = lr->left_child();
    lrr = lr->right_child();
    rrl = rr->left_child();
    rrr = rr->right_child();
    rrll = rrl->left_child();
    rrlr = rrl->right_child();
    rrrl = rrr->left_child();
    rrrr = rrr->right_child();

    ASSERT_EQ(root->label(), "root");
    ASSERT_EQ(r->label(), "r");
    ASSERT_EQ(l->label(), "l");
    ASSERT_EQ(rr->label(), "rr");
    ASSERT_EQ(rl->label(), "rl");
    ASSERT_EQ(lr->label(), "lr");
    ASSERT_EQ(ll->label(), "ll");
    ASSERT_EQ(rrl->label(), "rrl");
    ASSERT_EQ(rrr->label(), "rrr");
    ASSERT_EQ(lrl->label(), "lrl");
    ASSERT_EQ(lrr->label(), "lrr");
    ASSERT_EQ(rrrl->label(), "rrrl");
    ASSERT_EQ(rrrr->label(), "rrrr");
    ASSERT_EQ(rrll->label(), "rrll");
    ASSERT_EQ(rrlr->label(), "rrlr");

    // Ensure that our test cases below (MSATreeCompression) do not assume different node ids
    ASSERT_EQ(root->node_id(), 14);
    ASSERT_EQ(r->node_id(), 13);
    ASSERT_EQ(l->node_id(), 9);
    ASSERT_EQ(rr->node_id(), 12);
    ASSERT_EQ(rl->node_id(), 3);
    ASSERT_EQ(lr->node_id(), 8);
    ASSERT_EQ(ll->node_id(), 0);
    ASSERT_EQ(rrl->node_id(), 10);
    ASSERT_EQ(rrr->node_id(), 11);
    ASSERT_EQ(lrl->node_id(), 1);
    ASSERT_EQ(lrr->node_id(), 2);
    ASSERT_EQ(rrrl->node_id(), 6);
    ASSERT_EQ(rrrr->node_id(), 7);
    ASSERT_EQ(rrll->node_id(), 4);
    ASSERT_EQ(rrlr->node_id(), 5);
}

TEST_F(PllTreeTest, GarbageIn) {
    PllTree::Node * root, * inner, * leaf;
    root = simple_tree->root();
    inner = root->left_child();
    leaf = inner->right_child();
    assert(root && inner && leaf);

    // label() when label not set
    ASSERT_ANY_THROW(another_tree->root()->label());

    // parent on root
    ASSERT_EQ(root->parent(), nullptr);

    // left_child and right_child on leafs
    ASSERT_EQ(leaf->right_child(), nullptr);
    ASSERT_EQ(leaf->left_child(), nullptr);
}

TEST_F(PllTreeTest, Traversal) {
    // Check if the order in which the nodes are visited is correct

    // libpll only supports strictly bifurcating trees
    //auto snake_tree = make_shared<PllTree>("(((One)Two)Tree)Four;");

    vector<string> simple_preorder = {"root", "innerA", "bat", "cow", "innerB", "elk", "fox"};
    vector<string> simple_inorder = {"bat", "innerA", "cow", "root", "elk", "innerB", "fox"};
    vector<string> simple_postorder = {"bat", "cow", "innerA", "elk", "fox", "innerB", "root"};
    vector<string> snake_preorder = {"Four", "Tree", "Two", "One"};
    vector<string> elaborate_preorder = {"root", "l", "ll", "lr", "lrl", "lrr", "r", "rl", "rr", "rrl", "rrll", "rrlr", "rrr", "rrrl", "rrrr"};

    size_t idx;
    
    idx = 0;
    simple_tree->traverse(decorate_for_success([&idx, &simple_preorder](PllTree::Node * node) -> void {
        ASSERT_LT(idx, simple_preorder.size());
        ASSERT_EQ(simple_preorder[idx++], node->label());
    }), PllTree::TRAVERSAL_MODE::PREORDER);

    idx = 0;
    simple_tree->traverse(decorate_for_success([&idx, &simple_inorder](PllTree::Node * node) -> void {
        ASSERT_LT(idx, simple_inorder.size());
        //cout << simple_inorder[idx] << " ?= " << node->label() << endl;
        ASSERT_EQ(simple_inorder[idx++], node->label());
    }), PllTree::TRAVERSAL_MODE::INORDER);

    idx = 0;
    simple_tree->traverse(decorate_for_success([&idx, &simple_postorder](PllTree::Node * node) -> void {
        ASSERT_LT(idx, simple_postorder.size());
        ASSERT_EQ(simple_postorder[idx++], node->label());
    }), PllTree::TRAVERSAL_MODE::POSTORDER);

    /* see above
    idx = 0;
    snake_tree->traverse(decorate_for_success([&idx, &snake_preorder](PllTree::Node * node) -> void {
        ASSERT_LT(idx, snake_preorder.size());
        cout << "snake " << snake_preorder[idx] << endl;
        ASSERT_EQ(snake_preorder[idx++], node->label());
    }), PllTree::TRAVERSAL_MODE::PREORDER);
    */

    idx = 0;
    elaborate_tree->traverse(decorate_for_success([&idx, &elaborate_preorder](PllTree::Node * node) -> void {
        ASSERT_LT(idx, elaborate_preorder.size());
        ASSERT_EQ(elaborate_preorder[idx++], node->label());
    }), PllTree::TRAVERSAL_MODE::PREORDER);
}

TEST(NodeStatesTest, Basics) {
    MSATreeCompression::NodeStates node_states(10);

    // Does the iterator work?
    size_t i = 0;
    for(auto state: node_states) {
        RAXML_UNUSED(state);
        i++;
    }
    ASSERT_EQ(i, node_states.size());

    // Everything initialized to invalid?
    for(auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::INVALID);
    }

    // Getting and setting using operator[]
    for(size_t set_idx = 0; set_idx < node_states.size(); set_idx++) {
        node_states[set_idx] = IUPAC_DNA_16BIT::A;
        for(size_t chk_idx = 0; chk_idx < node_states.size(); chk_idx++) {                
            if(chk_idx <= set_idx) {
                ASSERT_EQ(node_states[chk_idx], IUPAC_DNA_16BIT::A);
            } else {
                ASSERT_EQ(node_states[chk_idx], IUPAC_DNA_16BIT::INVALID);
            }
        }
    }

    pll_rnode_t pll_node;
    pll_node.node_index = 1;
    PllTree::Node node(&pll_node, reinterpret_cast<PllTree*>(0xFF)); // We do not access the tree
    ASSERT_EQ(node.node_id(), 1);

    ASSERT_EQ(&node_states[node.node_id()], &node_states[node]);
    ASSERT_EQ(node_states[node], IUPAC_DNA_16BIT::A);
    node_states[node] = IUPAC_DNA_16BIT::K;
    ASSERT_EQ(node_states[node], IUPAC_DNA_16BIT::K);

    // size
    ASSERT_EQ(node_states.size(), 10);

    // choosing random states - IUPAC 4 bit
    ASSERT_EQ(node_states.choose_random_state(IUPAC_DNA::A), IUPAC_DNA::A);
    ASSERT_EQ(node_states.choose_random_state(IUPAC_DNA::C), IUPAC_DNA::C);
    ASSERT_EQ(node_states.choose_random_state(IUPAC_DNA::GAP), IUPAC_DNA::GAP);

    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA::A)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA::G)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA::M)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA::K)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA::V)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA::N)));

    vector<IUPAC_DNA> k_states = {IUPAC_DNA::G, IUPAC_DNA::T};
    ASSERT_THAT(k_states, Contains(node_states.choose_random_state(IUPAC_DNA::K)));

    vector<IUPAC_DNA> b_states = {IUPAC_DNA::G, IUPAC_DNA::T, IUPAC_DNA::C};
    ASSERT_THAT(b_states, Contains(node_states.choose_random_state(IUPAC_DNA::B)));

    vector<IUPAC_DNA> n_states = {IUPAC_DNA::G, IUPAC_DNA::T, IUPAC_DNA::C, IUPAC_DNA::A};
    ASSERT_THAT(n_states, Contains(node_states.choose_random_state(IUPAC_DNA::N)));

    vector<IUPAC_DNA> m_states = {IUPAC_DNA::A, IUPAC_DNA::C};
    ASSERT_THAT(m_states, Contains(node_states.choose_random_state(IUPAC_DNA::M)));

    vector<IUPAC_DNA> v_states = {IUPAC_DNA::G, IUPAC_DNA::C, IUPAC_DNA::A};
    ASSERT_THAT(v_states, Contains(node_states.choose_random_state(IUPAC_DNA::V)));

    // choosing random states - IUPAC 16 bit
    ASSERT_EQ(node_states.choose_random_state(IUPAC_DNA_16BIT::A), IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states.choose_random_state(IUPAC_DNA_16BIT::C), IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states.choose_random_state(IUPAC_DNA_16BIT::GAP), IUPAC_DNA_16BIT::GAP);

    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA_16BIT::A)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA_16BIT::G)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA_16BIT::M)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA_16BIT::K)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA_16BIT::V)));
    ASSERT_TRUE(is_single_state(node_states.choose_random_state(IUPAC_DNA_16BIT::N)));

    vector<IUPAC_DNA_16BIT> k_states_16bit = {IUPAC_DNA_16BIT::G, IUPAC_DNA_16BIT::T};
    ASSERT_THAT(k_states_16bit, Contains(node_states.choose_random_state(
        IUPAC_DNA_16BIT::G | IUPAC_DNA_16BIT::T)));

    vector<IUPAC_DNA_16BIT> b_states_16bit = {IUPAC_DNA_16BIT::G, IUPAC_DNA_16BIT::T, IUPAC_DNA_16BIT::C};
    ASSERT_THAT(b_states_16bit, Contains(node_states.choose_random_state(
        IUPAC_DNA_16BIT::G | IUPAC_DNA_16BIT::T | IUPAC_DNA_16BIT::C)));

    vector<IUPAC_DNA_16BIT> n_states_16bit = {IUPAC_DNA_16BIT::G, IUPAC_DNA_16BIT::T, IUPAC_DNA_16BIT::C, IUPAC_DNA_16BIT::A};
    ASSERT_THAT(n_states_16bit, Contains(node_states.choose_random_state(
        IUPAC_DNA_16BIT::T | IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::A)));

    vector<IUPAC_DNA_16BIT> m_states_16bit = {IUPAC_DNA_16BIT::A, IUPAC_DNA_16BIT::C};
    ASSERT_THAT(m_states_16bit, Contains(node_states.choose_random_state(
        IUPAC_DNA_16BIT::A | IUPAC_DNA_16BIT::C)));

    vector<IUPAC_DNA_16BIT> v_states_16bit = {IUPAC_DNA_16BIT::G, IUPAC_DNA_16BIT::C, IUPAC_DNA_16BIT::A};
    ASSERT_THAT(v_states_16bit, Contains(node_states.choose_random_state(
        IUPAC_DNA_16BIT::G | IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::A)));

    // choosing random state - node addressed

    node_states[node.node_id()] = IUPAC_DNA_16BIT::G | IUPAC_DNA_16BIT::T | IUPAC_DNA_16BIT::C;
    ASSERT_THAT(b_states_16bit, Contains(node_states.choose_random_state(node)));

    node_states[node.node_id()] = IUPAC_DNA_16BIT::A | IUPAC_DNA_16BIT::C;
    ASSERT_THAT(m_states_16bit, Contains(node_states.choose_random_state(node)));

    node_states[node.node_id()] = IUPAC_DNA_16BIT::G | IUPAC_DNA_16BIT::C | IUPAC_DNA_16BIT::A;
    ASSERT_THAT(v_states_16bit, Contains(node_states.choose_random_state(node)));

    node_states[node.node_id()] = IUPAC_DNA_16BIT::C;
    ASSERT_EQ(node_states[node], IUPAC_DNA_16BIT::C);

    node_states[node.node_id()] = IUPAC_DNA_16BIT::D | IUPAC_DNA_16BIT::G;
    IUPAC_DNA_16BIT choosen_state = node_states.choose_random_state(node);
    ASSERT_THAT(choosen_state, AnyOf(Eq(IUPAC_DNA_16BIT::D), Eq(IUPAC_DNA_16BIT::G)));

    node_states[node.node_id()] = IUPAC_DNA_16BIT::R | IUPAC_DNA_16BIT::T | IUPAC_DNA_16BIT::N;
    choosen_state = node_states.choose_random_state(node);
    ASSERT_THAT(choosen_state, AnyOf(Eq(IUPAC_DNA_16BIT::R), Eq(IUPAC_DNA_16BIT::T), Eq(IUPAC_DNA_16BIT::N)));

    // check if all states are single states - IUPAC_DNA_16BIT_16BIT
    for(auto state: {IUPAC_DNA_16BIT::A, IUPAC_DNA_16BIT::B, IUPAC_DNA_16BIT::M, IUPAC_DNA_16BIT::N}) {
        for(IUPAC_DNA_16BIT & node: node_states) {
            node = state;
        }
        assert(node_states[3] == state);
        ASSERT_TRUE(node_states.check_all_single_state());
    }

    node_states[4] = IUPAC_DNA_16BIT::K | IUPAC_DNA_16BIT::A;
    ASSERT_FALSE(node_states.check_all_single_state());

    // Does reset set all states back to INVALID?
    node_states.reset();
    for(auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::INVALID);
    }
}

TEST(NodeStateTest, InvalidOperations) {
    MSATreeCompression::NodeStates node_states(10);

    pll_rnode_t pll_node;
    pll_node.node_index = 1;
    PllTree::Node node(&pll_node, reinterpret_cast<PllTree*>(0xFF)); // We do not access the tree
    ASSERT_EQ(node.node_id(), 1);

    // Zero length initialization
    ASSERT_ANY_THROW(MSATreeCompression::NodeStates(0));

    // Using an invalid node id
    ASSERT_ANY_THROW(node_states[-1]);
    ASSERT_ANY_THROW(node_states[10]);
    ASSERT_ANY_THROW(node_states[11]);
    ASSERT_ANY_THROW(node_states[-1] = IUPAC_DNA_16BIT::G);
    ASSERT_ANY_THROW(node_states[10] = IUPAC_DNA_16BIT::G);
    ASSERT_ANY_THROW(node_states[11] = IUPAC_DNA_16BIT::G);

    pll_node.node_index = 11;
    ASSERT_EQ(node.node_id(), 11);
    ASSERT_EQ(node_states.size(), 10);
    ASSERT_ANY_THROW(node_states[node]);
    ASSERT_ANY_THROW(node_states[node] = IUPAC_DNA_16BIT::G);

    // Choosing random state from invalid state
    ASSERT_ANY_THROW(node_states.choose_random_state(IUPAC_DNA_16BIT::INVALID));
}

TEST(ChangeEncodingTest, Basics) {
    // Encoding
    const uint64_t max_edge_id = 1337;
    const uint8_t bits_for_edge_id = ceil(log2(1337));
    assert(bits_for_edge_id == 11);

    MSATreeCompression::ChangeEncoder change1(   0, max_edge_id, IUPAC_DNA::A, IUPAC_DNA::C);
    MSATreeCompression::ChangeEncoder change2(   1, max_edge_id, IUPAC_DNA::A, IUPAC_DNA::G);
    MSATreeCompression::ChangeEncoder change3( 202, max_edge_id, IUPAC_DNA::A, IUPAC_DNA::M);
    MSATreeCompression::ChangeEncoder change4(  32, max_edge_id, IUPAC_DNA::M, IUPAC_DNA::N);
    MSATreeCompression::ChangeEncoder change5(   4, max_edge_id, IUPAC_DNA::M, IUPAC_DNA::Y);
    MSATreeCompression::ChangeEncoder change6( 635, max_edge_id, IUPAC_DNA::W, IUPAC_DNA::K);
    MSATreeCompression::ChangeEncoder change7(1337, 2 * max_edge_id, IUPAC_DNA::K, IUPAC_DNA::D);
    MSATreeCompression::ChangeEncoder change8(  42, 2 * max_edge_id, IUPAC_DNA::V, IUPAC_DNA::T);

    // Trivial getters
    ASSERT_EQ(change1.edge_id(), 0);    
    ASSERT_EQ(change2.edge_id(), 1);    
    ASSERT_EQ(change3.edge_id(), 202);    

    ASSERT_EQ(change1.from_state(), IUPAC_DNA::A);
    ASSERT_EQ(change5.from_state(), IUPAC_DNA::M);
    ASSERT_EQ(change7.from_state(), IUPAC_DNA::K);
    
    ASSERT_EQ(change1.to_state(), IUPAC_DNA::C);
    ASSERT_EQ(change4.to_state(), IUPAC_DNA::N);
    ASSERT_EQ(change8.to_state(), IUPAC_DNA::T);

    // Check if all encoding have the expected length
    ASSERT_EQ(change1.bit_length(), bits_for_edge_id + 4);
    ASSERT_EQ(change2.bit_length(), bits_for_edge_id + 4);
    ASSERT_EQ(change3.bit_length(), bits_for_edge_id + 4);
    ASSERT_EQ(change4.bit_length(), bits_for_edge_id + 4);
    ASSERT_EQ(change5.bit_length(), bits_for_edge_id + 4);
    ASSERT_EQ(change6.bit_length(), bits_for_edge_id + 4);
    ASSERT_EQ(change7.bit_length(), bits_for_edge_id + 1 + 4);
    ASSERT_EQ(change8.bit_length(), bits_for_edge_id + 1 + 4);

    // Check if all encodings equal to what we expect
    ASSERT_EQ(change1.encoding(), 0b000000000001100);
    ASSERT_EQ(change2.encoding(), 0b000000000011001);
    ASSERT_EQ(change3.encoding(), 0b000110010100100);
    ASSERT_EQ(change4.encoding(), 0b000001000000011);
    ASSERT_EQ(change5.encoding(), 0b000000001001010);
    ASSERT_EQ(change6.encoding(), 0b010011110111001);
    ASSERT_EQ(change7.encoding(), 0b0101001110011000);
    ASSERT_EQ(change8.encoding(), 0b0000001010101111);

    // Is the Encoding ° Decdoding = Identity ?
    MSATreeCompression::ChangeDecoder decode_change1(change1.encoding(), max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change2(change2.encoding(), max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change3(change3.encoding(), max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change4(change4.encoding(), max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change5(change5.encoding(), max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change6(change6.encoding(), max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change7(change7.encoding(), 2 * max_edge_id);
    MSATreeCompression::ChangeDecoder decode_change8(change8.encoding(), 2 * max_edge_id);

    ASSERT_EQ(decode_change1.bit_length(), change1.bit_length());
    ASSERT_EQ(decode_change2.bit_length(), change2.bit_length());
    ASSERT_EQ(decode_change3.bit_length(), change3.bit_length());
    ASSERT_EQ(decode_change4.bit_length(), change4.bit_length());
    ASSERT_EQ(decode_change5.bit_length(), change5.bit_length());
    ASSERT_EQ(decode_change6.bit_length(), change6.bit_length());
    ASSERT_EQ(decode_change7.bit_length(), change7.bit_length());
    ASSERT_EQ(decode_change8.bit_length(), change8.bit_length());

    ASSERT_EQ(decode_change1.apply_change(change1.to_state()), change1.from_state());
    ASSERT_EQ(decode_change2.apply_change(change2.to_state()), change2.from_state());
    ASSERT_EQ(decode_change3.apply_change(change3.to_state()), change3.from_state());
    ASSERT_EQ(decode_change4.apply_change(change4.to_state()), change4.from_state());
    ASSERT_EQ(decode_change5.apply_change(change5.to_state()), change5.from_state());
    ASSERT_EQ(decode_change6.apply_change(change6.to_state()), change6.from_state());
    ASSERT_EQ(decode_change7.apply_change(change7.to_state()), change7.from_state());
    ASSERT_EQ(decode_change8.apply_change(change8.to_state()), change8.from_state());
    
    ASSERT_EQ(decode_change1.apply_change(change1.from_state()), change1.to_state());
    ASSERT_EQ(decode_change2.apply_change(change2.from_state()), change2.to_state());
    ASSERT_EQ(decode_change3.apply_change(change3.from_state()), change3.to_state());
    ASSERT_EQ(decode_change4.apply_change(change4.from_state()), change4.to_state());
    ASSERT_EQ(decode_change5.apply_change(change5.from_state()), change5.to_state());
    ASSERT_EQ(decode_change6.apply_change(change6.from_state()), change6.to_state());
    ASSERT_EQ(decode_change7.apply_change(change7.from_state()), change7.to_state());
    ASSERT_EQ(decode_change8.apply_change(change8.from_state()), change8.to_state());

    ASSERT_EQ(decode_change1.edge_id(), change1.edge_id());
    ASSERT_EQ(decode_change2.edge_id(), change2.edge_id());
    ASSERT_EQ(decode_change3.edge_id(), change3.edge_id());
    ASSERT_EQ(decode_change4.edge_id(), change4.edge_id());
    ASSERT_EQ(decode_change5.edge_id(), change5.edge_id());
    ASSERT_EQ(decode_change6.edge_id(), change6.edge_id());
    ASSERT_EQ(decode_change7.edge_id(), change7.edge_id());
    ASSERT_EQ(decode_change8.edge_id(), change8.edge_id());
}

TEST(ChangeEncodingTest, GarbageIn) {
    // Invalid constructor arguments
    ASSERT_ANY_THROW(MSATreeCompression::ChangeEncoder(-1, 10, IUPAC_DNA::A, IUPAC_DNA::A));
    ASSERT_ANY_THROW(MSATreeCompression::ChangeEncoder(11, 10, IUPAC_DNA::A, IUPAC_DNA::A));
    ASSERT_ANY_THROW(MSATreeCompression::ChangeEncoder( 1, 10, IUPAC_DNA::INVALID, IUPAC_DNA::A));
    ASSERT_ANY_THROW(MSATreeCompression::ChangeEncoder( 1, 10, IUPAC_DNA::A, IUPAC_DNA::INVALID));
    ASSERT_ANY_THROW(MSATreeCompression::ChangeDecoder(0b1110000, 2)); // Higher bits set
    ASSERT_ANY_THROW(MSATreeCompression::ChangeDecoder(0b0110000, 2));  // edge id too large
}

class MSATreeCompressionTest : public ::testing::Test {
    protected:
    string simple_tree_str = "((bat,cow)innerA,(elk,fox)innerB)root;";
    string another_tree_str = "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;";
    string elaborate_tree_str = "((ll,(lrl,lrr)lr)l,(rl,((rrll,rrlr)rrl,(rrrl,rrrr)rrr)rr)r)root;";

    virtual void SetUp() {}
};

TEST_F(MSATreeCompressionTest, Basics) {
    // Is num_sequences set correctly?
    MSATreeCompression simple_tree_compression(simple_tree_str);
    MSATreeCompression another_tree_compression(another_tree_str);
    MSATreeCompression elaborate_tree_compression(elaborate_tree_str);

    ASSERT_EQ(simple_tree_compression.num_sequences(), 4);
    ASSERT_EQ(another_tree_compression.num_sequences(), 5);
    ASSERT_EQ(elaborate_tree_compression.num_sequences(), 8);

       // None of the above changed the number of sequences in the tree
    ASSERT_EQ(simple_tree_compression.num_sequences(), 4);
    ASSERT_EQ(another_tree_compression.num_sequences(), 5);
    ASSERT_EQ(elaborate_tree_compression.num_sequences(), 8);
}

TEST_F(MSATreeCompressionTest, GarbageIn) {
    // Tree and MSA have unequal number of tips
    MSATreeCompression simple_tree_compression(simple_tree_str);

    MSA msa1, msa2;
    msa1.append("AAACCGCT", "bat");
    msa1.append("ACAGTTCT", "cow");
    msa1.append("AAAGATCA", "elk");
    ASSERT_ANY_THROW(simple_tree_compression.compress(msa1));    

    msa2.append("AAACCGCT", "bat");
    msa2.append("ACAGTTCT", "cow");
    msa2.append("AAAGATCA", "elk");
    msa2.append("AAATCTCA", "fox"); 
    msa2.append("AAATCTCA", "reindeer"); 
    ASSERT_ANY_THROW(simple_tree_compression.compress(msa2));
}

TEST_F(MSATreeCompressionTest, CompressionAndDecompression) {
    // decompress ° compress == Identity?
    // Is the compression' size the expected one?

    // Simple tree
    MSATreeCompression simple_tree_compression(simple_tree_str);

    MSA simple_msa;
    simple_msa.append("AAACCGCT", "bat");
    simple_msa.append("ACAGTTCT", "cow");
    simple_msa.append("AAAGATCA", "elk");
    simple_msa.append("AAATCTCA", "fox"); 

    RangeList simple_all_sites = { SiteRange(0, 8) };
    RangeList simple_partial = { SiteRange(1, 3) };
    RangeList simple_multi_range = { SiteRange(0, 1), SiteRange(1, 4), SiteRange(5, 3)};

    size_t compressed_size_in_bytes = simple_tree_compression.compress(simple_msa);
    auto simple_decompressed_msa = simple_tree_compression.decompress(simple_all_sites);
    ASSERT_EQ(simple_decompressed_msa->at("bat"), simple_msa.at("bat"));
    ASSERT_EQ(simple_decompressed_msa->at("cow"), simple_msa.at("cow"));
    ASSERT_EQ(simple_decompressed_msa->at("elk"), simple_msa.at("elk"));
    ASSERT_EQ(simple_decompressed_msa->at("fox"), simple_msa.at("fox"));

    // 7 changes -> 49 bits -> ~7 bytes
    // + 8 bytes for root sequence
    // = 11 bytes
    ASSERT_EQ(compressed_size_in_bytes, 11);

    auto simple_partial_decompressed_msa = simple_tree_compression.decompress(simple_partial);
    ASSERT_EQ(simple_partial_decompressed_msa->at("bat"), "AAC");
    ASSERT_EQ(simple_partial_decompressed_msa->at("cow"), "CAG");
    ASSERT_EQ(simple_partial_decompressed_msa->at("elk"), "AAG");
    ASSERT_EQ(simple_partial_decompressed_msa->at("fox"), "AAT");

    auto simple_multi_range_decompressed_msa = simple_tree_compression.decompress(simple_multi_range);
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("bat"), simple_msa.at("bat"));
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("cow"), simple_msa.at("cow"));
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("elk"), simple_msa.at("elk"));
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("fox"), simple_msa.at("fox"));

    // Another tree
    MSATreeCompression another_tree_compression(another_tree_str);

    MSA another_msa;
    another_msa.append("TACGTT", "One");
    another_msa.append("TAAGGT", "Two");
    another_msa.append("TATGGT", "Three");
    another_msa.append("TACTGT", "Four");
    another_msa.append("TACGTT", "Five");

    RangeList another_sites = { SiteRange(1, 4) };

    compressed_size_in_bytes = another_tree_compression.compress(another_msa);
    auto another_decompressed_msa = another_tree_compression.decompress(another_sites);
    ASSERT_EQ(another_decompressed_msa->at("One"), another_msa.at("One").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Two"), another_msa.at("Two").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Three"), another_msa.at("Three").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Four"), another_msa.at("Four").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Five"), another_msa.at("Five").substr(1, 4));

    // 5 changes -> 35 bits -> ~5 bytes
    // 3 + bytes for root sequence
    // = 8 bytes
    ASSERT_EQ(compressed_size_in_bytes, 8);

    // Elaborate tree
    MSATreeCompression elaborate_tree_compression(elaborate_tree_str);

    MSA elaborate_msa;
    elaborate_msa.append("CATGATAT", "ll");
    elaborate_msa.append("CATGATTT", "lrl");
    elaborate_msa.append("CCTGACTA", "lrr");
    elaborate_msa.append("CTGAACGA", "rl");
    elaborate_msa.append("CTGAACCG", "rrll");
    elaborate_msa.append("CCCCAAAA", "rrlr");
    elaborate_msa.append("CGGAAGTA", "rrrl");
    elaborate_msa.append("CATCACGT", "rrrr");

    RangeList elaborate_all_sites = { SiteRange(0, 8) };

    compressed_size_in_bytes = elaborate_tree_compression.compress(elaborate_msa);
    auto elaborate_decompressed_msa = elaborate_tree_compression.decompress(elaborate_all_sites);
    ASSERT_EQ(elaborate_decompressed_msa->at("ll"), elaborate_msa.at("ll"));
    ASSERT_EQ(elaborate_decompressed_msa->at("lrl"), elaborate_msa.at("lrl"));
    ASSERT_EQ(elaborate_decompressed_msa->at("lrr"), elaborate_msa.at("lrr"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rl"), elaborate_msa.at("rl"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrll"), elaborate_msa.at("rrll"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrlr"), elaborate_msa.at("rrlr"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrrl"), elaborate_msa.at("rrrl"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrrr"), elaborate_msa.at("rrrr"));

    // 24 changes -> 192 bits -> 24 bytes
    // + 4 byte for the root sequence
    // = 28 bytes
    ASSERT_EQ(compressed_size_in_bytes, 28);
}

TEST_F(MSATreeCompressionTest, CompressionAndDecompression16Bit) {
    // decompress ° compress == Identity?
    // Is the compression' size the expected one?

    // Simple tree
    MSATreeCompression simple_tree_compression(simple_tree_str);

    MSA simple_msa;
    simple_msa.append("---CCKCT", "bat");
    simple_msa.append("-C-KTTCT", "cow");
    simple_msa.append("---K-TC-", "elk");
    simple_msa.append("---TCTC-", "fox"); 

    RangeList simple_all_sites = { SiteRange(0, 8) };
    RangeList simple_partial = { SiteRange(1, 3) };
    RangeList simple_multi_range = { SiteRange(0, 1), SiteRange(1, 4), SiteRange(5, 3)};

    size_t compressed_size_in_bytes = simple_tree_compression.compress(simple_msa);
    auto simple_decompressed_msa = simple_tree_compression.decompress(simple_all_sites);
    ASSERT_EQ(simple_decompressed_msa->at("bat"), simple_msa.at("bat"));
    ASSERT_EQ(simple_decompressed_msa->at("cow"), simple_msa.at("cow"));
    ASSERT_EQ(simple_decompressed_msa->at("elk"), simple_msa.at("elk"));
    ASSERT_EQ(simple_decompressed_msa->at("fox"), simple_msa.at("fox"));

    // 7 changes -> 49 bits -> ~7 bytes
    // + 8 bytes for root sequence
    // = 11 bytes
    ASSERT_EQ(compressed_size_in_bytes, 11);

    auto simple_partial_decompressed_msa = simple_tree_compression.decompress(simple_partial);
    ASSERT_EQ(simple_partial_decompressed_msa->at("bat"), "--C");
    ASSERT_EQ(simple_partial_decompressed_msa->at("cow"), "C-K");
    ASSERT_EQ(simple_partial_decompressed_msa->at("elk"), "--K");
    ASSERT_EQ(simple_partial_decompressed_msa->at("fox"), "--T");

    auto simple_multi_range_decompressed_msa = simple_tree_compression.decompress(simple_multi_range);
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("bat"), simple_msa.at("bat"));
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("cow"), simple_msa.at("cow"));
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("elk"), simple_msa.at("elk"));
    ASSERT_EQ(simple_multi_range_decompressed_msa->at("fox"), simple_msa.at("fox"));

    // Another tree
    MSATreeCompression another_tree_compression(another_tree_str);

    MSA another_msa;
    another_msa.append("NAVGNN", "One");
    another_msa.append("NAAGGN", "Two");
    another_msa.append("NANGGN", "Three");
    another_msa.append("NAVNGN", "Four");
    another_msa.append("NAVGNN", "Five");

    RangeList another_sites = { SiteRange(1, 4) };

    compressed_size_in_bytes = another_tree_compression.compress(another_msa);
    auto another_decompressed_msa = another_tree_compression.decompress(another_sites);
    ASSERT_EQ(another_decompressed_msa->at("One"), another_msa.at("One").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Two"), another_msa.at("Two").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Three"), another_msa.at("Three").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Four"), another_msa.at("Four").substr(1, 4));
    ASSERT_EQ(another_decompressed_msa->at("Five"), another_msa.at("Five").substr(1, 4));

    // 5 changes -> 35 bits -> ~5 bytes
    // 3 + bytes for root sequence
    // = 8 bytes
    ASSERT_EQ(compressed_size_in_bytes, 8);

    // Elaborate tree
    MSATreeCompression elaborate_tree_compression(elaborate_tree_str);

    MSA elaborate_msa;
    elaborate_msa.append("C-TD-T-T", "ll");
    elaborate_msa.append("C-TD-TTT", "lrl");
    elaborate_msa.append("CCTD-CT-", "lrr");
    elaborate_msa.append("CTD--CD-", "rl");
    elaborate_msa.append("CTD--CCD", "rrll");
    elaborate_msa.append("CCCC----", "rrlr");
    elaborate_msa.append("CDD--DT-", "rrrl");
    elaborate_msa.append("C-TC-CDT", "rrrr");

    RangeList elaborate_all_sites = { SiteRange(0, 8) };

    compressed_size_in_bytes = elaborate_tree_compression.compress(elaborate_msa);
    auto elaborate_decompressed_msa = elaborate_tree_compression.decompress(elaborate_all_sites);
    ASSERT_EQ(elaborate_decompressed_msa->at("ll"), elaborate_msa.at("ll"));
    ASSERT_EQ(elaborate_decompressed_msa->at("lrl"), elaborate_msa.at("lrl"));
    ASSERT_EQ(elaborate_decompressed_msa->at("lrr"), elaborate_msa.at("lrr"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rl"), elaborate_msa.at("rl"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrll"), elaborate_msa.at("rrll"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrlr"), elaborate_msa.at("rrlr"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrrl"), elaborate_msa.at("rrrl"));
    ASSERT_EQ(elaborate_decompressed_msa->at("rrrr"), elaborate_msa.at("rrrr"));

    // 24 changes -> 192 bits -> 24 bytes
    // + 4 byte for the root sequence
    // = 28 bytes
    ASSERT_EQ(compressed_size_in_bytes, 28);
}

TEST_F(MSATreeCompressionTest, ComputingAncestralStatesSimple) {
    // Are the ancestral states choosen correctly?
    MSA msa;
    msa.append("AACC", "bat");
    msa.append("ACGT", "cow");
    msa.append("AAGA", "elk");
    msa.append("AATC", "fox"); 

    MSATreeCompression::NodeStates node_states(7);
    MSATreeCompression msa_tree_compression(simple_tree_str);

    msa_tree_compression.build_ancestral_states(msa, 0, node_states);
    for (auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::A);
    }
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 1, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::A);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 2, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::G);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 3, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::C);
}

TEST_F(MSATreeCompressionTest, ComputingAncestralStatesSimple16Bit) {
    // Are the ancestral states choosen correctly?
    MSA msa;
    msa.append("VVMM", "bat");
    msa.append("VMGT", "cow");
    msa.append("VVGV", "elk");
    msa.append("VVTM", "fox"); 

    MSATreeCompression::NodeStates node_states(7);
    MSATreeCompression msa_tree_compression(simple_tree_str);

    msa_tree_compression.build_ancestral_states(msa, 0, node_states);
    for (auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::V);
    }
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 1, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::V);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::M);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::V);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::V);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::V);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::V);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::V);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 2, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::M);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::G);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 3, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::M);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::V);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::M);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::M);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::M);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::M);
}

TEST_F(MSATreeCompressionTest, ComputingAncestralStatesAnother) {
    // Are the ancestral states choosen correctly?
    MSA msa;
    msa.append("ACAA", "One");
    msa.append("ATAA", "Two");
    msa.append("AACC", "Three");
    msa.append("ACTT", "Four");
    msa.append("ATAA", "Five");

    MSATreeCompression::NodeStates node_states(9);
    MSATreeCompression msa_tree_compression(another_tree_str);

    msa_tree_compression.build_ancestral_states(msa, 0, node_states);
    for (auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::A);
    }
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 1, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::T);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 2, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::A);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 3, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::A);
}

TEST_F(MSATreeCompressionTest, ComputingAncestralStatesAnother16Bit) {
    // Are the ancestral states choosen correctly?
    MSA msa;
    msa.append("KCKK", "One");
    msa.append("K-KK", "Two");
    msa.append("KKCC", "Three");
    msa.append("KC--", "Four");
    msa.append("K-KK", "Five");

    MSATreeCompression::NodeStates node_states(9);
    MSATreeCompression msa_tree_compression(another_tree_str);

    msa_tree_compression.build_ancestral_states(msa, 0, node_states);
    for (auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::K);
    }
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 1, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::GAP);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 2, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::K);
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 3, node_states);
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::K);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::K);
}

TEST_F(MSATreeCompressionTest, ComputingAncestralStatesElaborate) {
    // Are the ancestral states choosen correctly?
    MSA msa;
    msa.append("CA", "ll");
    msa.append("CA", "lrl");
    msa.append("CC", "lrr");
    msa.append("CT", "rl");
    msa.append("CT", "rrll");
    msa.append("CC", "rrlr");
    msa.append("CG", "rrrl");
    msa.append("CA", "rrrr");

    MSATreeCompression::NodeStates node_states(15);
    MSATreeCompression msa_tree_compression(elaborate_tree_str);

    msa_tree_compression.build_ancestral_states(msa, 0, node_states);
    for (auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::C);
    }
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 1, node_states);
    // See ElaborateTree test on how the node ids are assigned to nodes
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::C);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[9], IUPAC_DNA_16BIT::A);
    ASSERT_EQ(node_states[10], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[11], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[12], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[13], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[14], IUPAC_DNA_16BIT::T);
    node_states.reset();
}

TEST_F(MSATreeCompressionTest, ComputingAncestralStatesElaborate16Bit) {
    // Are the ancestral states choosen correctly?
    MSA msa;
    msa.append("N-", "ll");
    msa.append("N-", "lrl");
    msa.append("NN", "lrr");
    msa.append("NT", "rl");
    msa.append("NT", "rrll");
    msa.append("NN", "rrlr");
    msa.append("NG", "rrrl");
    msa.append("N-", "rrrr");

    MSATreeCompression::NodeStates node_states(15);
    MSATreeCompression msa_tree_compression(elaborate_tree_str);

    msa_tree_compression.build_ancestral_states(msa, 0, node_states);
    for (auto state: node_states) {
        ASSERT_EQ(state, IUPAC_DNA_16BIT::N);
    }
    node_states.reset();

    msa_tree_compression.build_ancestral_states(msa, 1, node_states);
    // See ElaborateTree test on how the node ids are assigned to nodes
    ASSERT_EQ(node_states[0], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[1], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[2], IUPAC_DNA_16BIT::N);
    ASSERT_EQ(node_states[3], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[4], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[5], IUPAC_DNA_16BIT::N);
    ASSERT_EQ(node_states[6], IUPAC_DNA_16BIT::G);
    ASSERT_EQ(node_states[7], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[8], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[9], IUPAC_DNA_16BIT::GAP);
    ASSERT_EQ(node_states[10], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[11], IUPAC_DNA_16BIT::GAP); // GAP is choosen over G (leftmost Bit)
    ASSERT_EQ(node_states[12], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[13], IUPAC_DNA_16BIT::T);
    ASSERT_EQ(node_states[14], IUPAC_DNA_16BIT::GAP); // GAP is choosen over T (leftmost Bit)
    node_states.reset();
}

void load_msa_and_build_parsimony_tree(string filename, PartitionedMSA& part_msa, Tree& tree) {
    // Load MSa file into memory
    MSA load_msa = msa_load_from_file(filename, FileFormat::autodetect);

    // Build the parsimony tree
    unsigned int score;

    unsigned int attrs = sysutil_simd_autodetect();
    attrs |= PLL_ATTRIB_PATTERN_TIP;

    part_msa.emplace_part_info("p1", DataType::dna, "GTR");
    part_msa.full_msa(std::move(load_msa));
    part_msa.split_msa();
    
    tree = Tree::buildParsimony(part_msa, 0, attrs, &score);
}

TEST(MSACompressionFromFile, CompressingPrim) {
    PartitionedMSA part_msa;
    Tree tree;
    
    load_msa_and_build_parsimony_tree("test/data/prim.phy", part_msa, tree);
    
    const MSA& msa = part_msa.full_msa();

    // Compress
    MSATreeCompression compressed_msa(tree);
    size_t compressed_size_in_bytes = compressed_msa.compress(msa);
    ASSERT_EQ(compressed_size_in_bytes, 1758); // fixed at what it was on first successfull compression
                                               // not hand-calculated as the values above

    // Decompress and check if decompression ° compression = identity 
    RangeList all_sites = { SiteRange(0, msa.num_sites()) };
    auto decompressed_msa = compressed_msa.decompress(all_sites);
    assert(msa.labels().size() == 12);
    for (string taxon: msa.labels()) {
        ASSERT_EQ(decompressed_msa->at(taxon), msa.at(taxon));
    }
}

TEST(MSACompressionFromFile, CompressingFusob) {
    PartitionedMSA part_msa;
    Tree tree;
    
    // TODO support RNA data and RNA+DNA data
    load_msa_and_build_parsimony_tree("test/data/fusob_dnaified.phy", part_msa, tree);
    
    const MSA& msa = part_msa.full_msa();

    // Compress
    MSATreeCompression compressed_msa(tree);
    size_t compressed_size_in_bytes = compressed_msa.compress(msa);
    ASSERT_EQ(compressed_size_in_bytes, 5728); // fixed at what it was on first successfull compression
                                               // not hand-calculated as the values above

    // Decompress and check if decompression ° compression = identity 
    RangeList all_sites = { SiteRange(0, msa.num_sites()) };
    auto decompressed_msa = compressed_msa.decompress(all_sites);
    assert(msa.labels().size() == 38);
    for (string taxon: msa.labels()) {
        string no_question_marks_sequence = msa.at(taxon);
        replace(no_question_marks_sequence.begin(), no_question_marks_sequence.end(), '?', '-');
        ASSERT_EQ(decompressed_msa->at(taxon), no_question_marks_sequence);
    }
}

TEST(MSACompressionFromFile, CompressingRbcl) {
    PartitionedMSA part_msa;
    Tree tree;
    
    load_msa_and_build_parsimony_tree("test/data/rbcl.phy", part_msa, tree);
    
    const MSA& msa = part_msa.full_msa();

    // Compress
    MSATreeCompression compressed_msa(tree);
    size_t compressed_size_in_bytes = compressed_msa.compress(msa);
    ASSERT_EQ(compressed_size_in_bytes, 27558); // fixed at what it was on first successfull compression
                                                // not hand-calculated as the values above

    // Decompress and check if decompression ° compression = identity 
    RangeList all_sites = { SiteRange(0, msa.num_sites()) };
    auto decompressed_msa = compressed_msa.decompress(all_sites);
    assert(msa.labels().size() == 436);
    for (string taxon: msa.labels()) {
        ASSERT_EQ(decompressed_msa->at(taxon), msa.at(taxon));
    }
}