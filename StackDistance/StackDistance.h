#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#include <unordered_map>

#define MISS_BAR 513

template <class U>
class Histogram
{
	std::map<U, int> bins;

public:
	Histogram() {};

	~Histogram() {};

	void sample(U x)
	{
		if (!bins[x]) {
			bins[x] = 1;
		}
		else
			++bins[x];
	}

	void print(std::ofstream * & f)
	{
		std::map<int, int>::iterator it;

		for (it = bins.begin(); it != bins.end(); ++it) {
			*f << "dist" << it->first << "\t\t";
			*f << it->second << std::endl;
		}
	}

	void changeAssoc()
	{

	}
};

template <class T>
class AvlNode
{
	//int holes;
	//friend class AvlTreeStack<T>;

public:
	/* this is the num of holes of this tree,
	including the holes of subtrees and the self interval. */
	int holes;
	/* this is the num of holes of entire right subtree. */
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

class AvlTreeStack
{
	AvlNode<long> * root;
	bool isRoot;
	Histogram<int> hist;
	std::map <uint64_t, long> addrMap;
	long index;
	int dist;

public:

	AvlTreeStack(long & v) : dist(0)
	{
		root = new AvlNode<long>(v);
		isRoot = true;
	}

	AvlTreeStack() : root(nullptr), dist(0), isRoot(false) {};

	void setRoot(long & v)
	{
		root = new AvlNode<long>(v);
	}

	~AvlTreeStack() { destroy(root); }

	void destroy(AvlNode<long> * & tree);

	void insert(AvlNode<long> * & tree, long & v);

	void insert(long & a);

	/*AvlNode<long> * & find(AvlNode<long> * & tree, int & v)
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

	std::pair<long, long> & findMin(AvlNode<long> * & tree);

	std::pair<long, long> & findMax(AvlNode<long> * & tree);

	void remove(AvlNode<long> * & tree, std::pair<long, long> & inter);

	void rotate(AvlNode<long> * & tree);

	void doubleRotate(AvlNode<long> * & tree);

	void balance(AvlNode<long> * & tree);

	void calStackDist(uint64_t addr)
	{
		long & value = addrMap[addr];

		++index;

		/* value is 0 under cold miss */
		if (!value) {
			hist.sample(MISS_BAR);
			value = index;
			return;
		}

		/* update b of last reference */
		if (value < index) {
			/* set a root of the tree */
			if (isRoot) {
				setRoot(value);
				isRoot = false;
			}
			else {
				/* insert a hole */
				insert(value);
				int stackDist = index - value - dist;
				stackDist = stackDist >= MISS_BAR ? MISS_BAR : stackDist;
				hist.sample(stackDist);
			}
		}

		value = index;
	}

	void print(std::string fname)
	{
		std::ofstream * file = new std::ofstream(fname);
		hist.print(file);
		file->close();
		delete file;
	}
};

class ListStack
{
	std::list<uint64_t> addrList;

	Histogram<int> hist;

public:
	ListStack() {};

	~ListStack() {};

	void calStackDist(uint64_t addr) 
	{
		std::list<uint64_t>::iterator it;
		/* start counting the distance from 1 */
		long distance = 1;
		bool found = false;

		for (it = addrList.begin(); it != addrList.end(); ++it)
		{
			if (*it == addr) {
				addrList.erase(it);
				hist.sample(distance);
				found = true;
				break;
			}
			else
				++distance;
		}

		/* cold miss */
		if (!found)
			hist.sample(MISS_BAR);

		if (addrList.size() >= MISS_BAR) {
			addrList.pop_back();
		}
		/* push new address to the head of the list */
		addrList.push_front(addr);
	};

	void print(std::string fname)
	{
		std::ofstream * file = new std::ofstream(fname);
		hist.print(file);
		file->close();
		delete file;
	}
};

class Reader
{
private:
	std::ifstream file;
	std::string line;
	AvlTreeStack avlTreeStack;
	ListStack listStack;

public:
	Reader(std::string _path);
};