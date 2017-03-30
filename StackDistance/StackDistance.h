#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

template <class T>
class AvlNode
{
	//int holes;
	//friend class AvlTree<T>;

public:
	int holes;
	int rHoles;
	std::pair<T, T> interval;
	int height;
	AvlNode<T> * left;
	AvlNode<T> * right;

	AvlNode(T & a) : holes(1), rHoles(0), height(0), left(nullptr), right(nullptr)
	{
		interval = std::make_pair(a, a);
	}

	AvlNode(AvlNode<T> & n) : holes(n.holes), rHoles(n.rHoles), height(n.height), left(n.left), right(n.right)
	{
		interval = n.interval;
	}

	~AvlNode()
	{
		delete left;
		delete right;
		left = nullptr;
		right = nullptr;
	}

	static int getHeight(AvlNode<T> * & node)
	{
		return node ? node->height : -1;
	}

	void updateHeight()
	{
		height = 1 + std::max(getHeight(left), getHeight(right));
	}

	void updateHoles()
	{
		rHoles = right ? right->holes : 0;

		holes =
			(left ? left->holes : 0) + rHoles +
			int(interval.second - interval.first) + 1;
	}
};

template <class L>
class AvlTree
{
	AvlNode<L> * root;
	friend class Reader;
public:
	int dist;

	AvlTree(L & v) : dist(0)
	{
		root = new AvlNode<L>(v);
	}

	AvlTree() : root(nullptr), dist(0) {};

	void setRoot(L & v)
	{
		root = new AvlNode<L>(v);
	}

	~AvlTree() { destroy(root); }

	void destroy(AvlNode<L> * & tree);

	void insert(AvlNode<L> * & tree, L & v);

	void insert(L & a);

	/*AvlNode<L> * & find(AvlNode<L> * & tree, int & v)
	{
	if (!tree)
	return nullptr;

	if (v < tree->holes)
	find(tree->left, v);
	else if (v > tree->holes)
	find(tree->right, v);
	else
	return tree;
	}*/

	std::pair<L, L> & findMin(AvlNode<L> * & tree);

	std::pair<L, L> & findMax(AvlNode<L> * & tree);

	void remove(AvlNode<L> * & tree, std::pair<L, L> & inter);

	void rotate(AvlNode<L> * & tree);

	void doubleRotate(AvlNode<L> * & tree);

	void balance(AvlNode<L> * & tree);
};

class Reader
{
private:
	std::ifstream file;
	std::string line;
	std::map <uint64_t, long> addrMap;
	long index;
	std::vector <bool> b;
	AvlTree<long> tree;

public:
	Reader(std::string _path);
};