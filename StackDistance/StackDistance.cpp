// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"

template <class L>
void AvlTree<L>::destroy(AvlNode<L> * & tree)
{
	if (!tree)
		return;

	destroy(tree->left);
	destroy(tree->right);
	delete tree;
	tree = nullptr;
}

template <class L>
void AvlTree<L>::insert(AvlNode<L> * & tree, L & v) {
	if (!tree) {
		tree = new AvlNode<L>(v);
		return;
	}

	std::pair<L, L> & interval = tree->interval;
	assert(!(v >= interval.first && v <= interval.second));

	AvlNode<L> * & n1 = tree->left;
	AvlNode<L> * & n2 = tree->right;

	if (tree->holes < 0)
		return;

	if (v == interval.first - 1) {
		dist += interval.second - v + tree->rHoles;
		interval.first = v;
		if (n1) {
			std::pair<L, L> & temp = findMax(n1);
			if (v == temp.second + 1) {
				interval.first = temp.first;
				remove(n1, temp);
			}
		}
	}
	else if (v == interval.second + 1) {
		dist += tree->rHoles;
		interval.second = v;
		if (n2) {
			std::pair<L, L> & temp = findMin(n2);
			if (v == temp.first - 1) {
				interval.second = temp.second;
				remove(n2, temp);
			}
		}
	}
	else if (v < interval.first - 1) {
		dist += tree->rHoles + interval.second - interval.first + 1;
		insert(n1, v);
	}
	else if (v > tree->interval.second + 1)
		insert(n2, v);

	balance(tree);
}

template <class L>
void AvlTree<L>::insert(L & a) { dist = 0; insert(root, a); }

template <class L>
std::pair<L, L> & AvlTree<L>::findMin(AvlNode<L> * & tree)
{
	assert(tree);

	if (!tree->left)
		return tree->interval;
	else
		findMin(tree->left);
}

template <class L>
std::pair<L, L> & AvlTree<L>::findMax(AvlNode<L> * & tree)
{
	assert(tree);

	if (!tree->right)
		return tree->interval;
	else
		findMax(tree->right);
}

template <class L>
void AvlTree<L>::remove(AvlNode<L> * & tree, std::pair<L, L> & inter)
{
	if (!tree)
		return;

	if (inter.first > tree->interval.first)
		remove(tree->right, inter);
	else if (inter.first < tree->interval.first)
		remove(tree->left, inter);
	else if (tree->left && tree->right) {
		std::pair<L, L> & temp = findMin(tree->right);
		remove(tree, temp);
	}
	else {
		AvlNode<L> * old = tree;
		if (!tree->left && !tree->right)
			tree = nullptr;
		else
			tree = new AvlNode<L>(*(tree->left ? tree->left : tree->right));

		delete old;
	}

	balance(tree);
}

template <class L>
void AvlTree<L>::rotate(AvlNode<L> * & tree)
{
	if (!tree)
		return;

	if (tree->holes < 0)
		return;

	AvlNode<L> * n1 = tree->left;
	AvlNode<L> * n2 = tree->right;

	if (n2 && AvlNode<L>::getHeight(n1) < AvlNode<L>::getHeight(n2)) {
		tree->right = n2->left;
		tree->updateHoles();
		tree->updateHeight();

		n2->left = tree;
		n2->updateHoles();
		n2->updateHeight();

		tree = n2;
	}
	else if (n1 && AvlNode<L>::getHeight(n1) > AvlNode<L>::getHeight(n2)) {
		tree->left = n1->right;
		tree->updateHoles();
		tree->updateHeight();

		n1->right = tree;
		n1->updateHoles();
		n1->updateHeight();

		tree = n1;
	}
}

template <class L>
void AvlTree<L>::doubleRotate(AvlNode<L> * & tree)
{
	if (!tree)
		return;

	if (tree->holes < 0)
		return;

	AvlNode<L> * & n1 = tree->left;
	AvlNode<L> * & n2 = tree->right;

	if (AvlNode<L>::getHeight(n1) < AvlNode<L>::getHeight(n2))
		rotate(n2);
	else
		rotate(n1);

	rotate(tree);
}

template <class L>
void AvlTree<L>::balance(AvlNode<L> * & tree)
{
	if (!tree)
		return;

	AvlNode<L> * & n1 = tree->left;
	AvlNode<L> * & n2 = tree->right;
	int err = AvlNode<L>::getHeight(n1) - AvlNode<L>::getHeight(n2);
	if (err > 1) {
		if (AvlNode<L>::getHeight(n1->left) >= AvlNode<L>::getHeight(n1->right))
			rotate(tree);
		else
			doubleRotate(tree);
	}
	else if (err < -1) {
		if (AvlNode<L>::getHeight(n2->left) <= AvlNode<L>::getHeight(n2->right))
			rotate(tree);
		else
			doubleRotate(tree);
	}

	tree->updateHeight();
	tree->updateHoles();
}

Reader::Reader(std::string _path) {
	file.open(_path, std::ios::in);
	if (file.fail()) {
		return;
	}

	std::string temp;
	index = 0;
	bool root = true;

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		uint64_t addr;
		lineStream >> temp;
		lineStream >> std::hex >> addr;

		/* check a valid address */
		if (!addr)
			continue;

		long & value = addrMap[addr];

		++index;

		/* value is LONG_MAX under cold miss */
		if (!value) {
			value = index;
			continue;
		}

		/* update b of last reference */
		if (value < index) {
			if (root) {
				tree.setRoot(value);
				root = false;
			}
			else {
				tree.insert(value);
				std::cout << "distance is " << index - value - tree.dist << std::endl;
			}
		}

		value = index;
	}
}

int main()
{

	Reader reader("E:\\ShareShen\\gem5-origin\\statistics-results\\traceFile.txt");

	return 0;
}