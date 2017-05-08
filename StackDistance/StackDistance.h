#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <list>
#include <float.h>
#define MISS_BAR 1025

inline double power(const double base,const int index);

inline double biDistribution(const int m, const int n, const double p);

template <class U, class T>
class Histogram
{
	std::map<U, T> bins;
	std::vector<double> binsVec;
	std::vector<double> setDistr;
	std::vector<double> transBins;

public:
	Histogram() {};

	~Histogram() {};

	void sample(U x);

	bool intoVector();

	void changeAssoc(const int & cap, const int & blk, const int & assoc);

	void print(std::ofstream & file);
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
	std::map <uint64_t, long> addrMap;
	long index;
	int dist;

public:
	/* hist contains the stack distance distribution */
	Histogram<int, int> hist;

	AvlTreeStack(long & v) : dist(0)
	{
		root = new AvlNode<long>(v);
		isRoot = true;
	}

	AvlTreeStack() : root(nullptr), dist(0), isRoot(false) {};

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

	void calStackDist(uint64_t addr);

	bool transHist(const int cap, const int blk, const int assoc);
};

class ListStack
{
	std::list<uint64_t> addrList;

	Histogram<int, int> hist;

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
};

class Reader
{
private:
	AvlTreeStack avlTreeStack;
	ListStack listStack;

public:
	Reader(std::string _path, std::string _pathout);
};