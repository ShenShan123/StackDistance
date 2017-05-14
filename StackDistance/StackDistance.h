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


#define MISS_BAR 1024 * 4 + 1

inline double power(const double base,const int index, double coef = 1.0);

inline double biDistribution(const int m, const int n, const double p);

template <class B, class T>
class Histogram
{
	std::vector<B> binsVec;
	std::vector<T> transBins;
	std::vector<double> assocDist;
	/* number of sampling */
	B samples;
	/* total number of miss references */
	B misses;
	/* miss rate */
	double missRate;

public:
	Histogram() : misses(0), missRate(0.0), 
		          binsVec(MISS_BAR + 1, 0), 
		          transBins(MISS_BAR + 1, 0),
				  assocDist(MISS_BAR + 1, 0) {}

	~Histogram() {};

	void sample(B x);

	void changeAssoc(const int & cap, const int & blk, const int & assoc);

	void transToStackDist()
	{
		std::vector<double> frac;
		B temp = 0;
		/* calculate the fraction of reuse distance
		   is greater than i */
		for (int i = 0; i < binsVec.size(); ++i) {
			temp += binsVec[i];
			frac.push_back((double)(samples - temp) / samples);
		}

		for (int i = 1; i < binsVec.size(); ++i) {
			double sumFrac = 0.0;
			for (int j = 1; j <= i; ++j)
				sumFrac += frac[j];
			transBins[(int)std::round(sumFrac)] += binsVec[i];
		}
		transBins[0] = binsVec[0];
	}

	void print(std::ofstream & file);

	void calMissRate(const int & assoc)
	{
		for (int i = assoc; i < transBins.size(); ++i)
			misses += (B)std::round(transBins[i]);
		missRate = (double)misses / samples;
	}

	//void calMissRate(const int & assoc)
	//{
	//	for (int i = assoc; i < assocDist.size(); ++i)
	//		misses += assocDist[i];
	//	missRate = (double)misses / samples;
	//}

	void calMissRate(const int & cap, const int & blk)
	{
		int blkNum = cap / blk;
		for (int i = blkNum; i < transBins.size(); ++i)
			misses += transBins[i];
		missRate = (double)misses / samples;
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
	std::map <uint64_t, long> addrMap;
	long index;
	int dist;

public:
	/* hist contains the stack distance distribution */
	Histogram<int, long double> hist;

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

	void transHist(const int cap, const int blk, const int assoc);
};

class ListStack
{
	std::list<uint64_t> addrList;

	Histogram<int, double> hist;

public:
	ListStack() {};

	~ListStack() {};

	void calStackDist(uint64_t addr);
};

class ReuseDist
{
	std::map<uint64_t, long> addrMap;
	long index;

public:
	Histogram<int, int> hist;

	ReuseDist() {};
	~ReuseDist() {};

	void calReuseDist(uint64_t addr)
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
			int reuseDist = index - value - 1;
			reuseDist = reuseDist >= MISS_BAR ? MISS_BAR : reuseDist;
			hist.sample(reuseDist);
		}

		value = index;
	}

	void transHist()
	{
		hist.transToStackDist();
	}
};

class Reader
{
private:
	AvlTreeStack avlTreeStack;
	ListStack listStack;
	ReuseDist reuseDist;

public:
	Reader(std::string _path, std::string _pathout);
};