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
#include <random>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <time.h>

//#define REUSE
#define STACK
//#define SAMPLE

#define MISS_BAR 1024 * 2

/* for recording distribution into a Histogram, 
   Accur is the accuracy of transforming calculation */
template <class B = int, class Accur = double>
class Histogram
{
	/* Histogram implemented by std::map */
	std::map<long, B> binsMap;
	/* Histogram implemented by std::vector */
	std::vector<B> binsVec;
	/* Histogram implemented by std::vector, 
	   after fully-to-set-associative cache transformation */
	std::vector<Accur> binsTra;
	/* number of sampling */
	B samples;
	/* total number of miss references */
	B misses;
	/* cache miss rate */
	Accur missRate;

public:
	Histogram() : misses(0), missRate(0.0) {};

	~Histogram() {};

	void sample(B x);
	/* complete the histogram with std::vector, to fast calculation. */
	bool mapToVector();
	/* use boost binomial distribution to do a fully-to-set-associative cache transformation */
	void fullyToSetAssoc(const int & cap, const int & blk, const int & assoc);
	/* use boost poisson distribution to do a fully-to-set-associative cache transformation */
	void fullyToSetAssoc_Poisson(const int & cap, const int & blk, const int & assoc);
	/* transform reuse distance distribution to stack distance distribution */
	void reuseDistToStackDist();

	void print(std::ofstream & file);
	/* calculate the miss rate for LRU set associative cache */
	void calMissRate(const int & assoc);

	void calMissRate(const int & cap, const int & blk);
	/* calculate the miss rate for PLRU set associative cache */
	void calMissRate(const int & assoc, const bool plru);
};

template <class Accur>
class AvlNode
{
	//int holes;
	friend class AvlTreeStack;
	/* this is the num of holes of this tree,
	including the holes of subtrees and the self interval. */
	int holes;
	/* this is the num of holes of entire right subtree. */
	int rHoles;
	std::pair<Accur, Accur> interval;
	int height;
	AvlNode<Accur> * left;
	AvlNode<Accur> * right;

public:
	AvlNode(Accur & a) : holes(1), rHoles(0), height(0), left(nullptr), right(nullptr)
	{
		interval = std::make_pair(a, a);
	}

	AvlNode(AvlNode<Accur> & n) : holes(n.holes), rHoles(n.rHoles), height(n.height), left(n.left), right(n.right)
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

	static int getHeight(AvlNode<Accur> * & node)
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

/* for calculating stack distance distribution via AVL Tree, with no sampling*/
class AvlTreeStack
{
	std::map <uint64_t, long> addrMap;
	AvlNode<long> * root;
	/* the index of refs in memory trace */
	long index;
	/* holes between current ref and last ref with same address */
	int curHoles;

public:
	AvlTreeStack(long & v) : index(0), curHoles(0)
	{
		root = new AvlNode<long>(v);
	}

	AvlTreeStack() : root(nullptr), index(0), curHoles(0) {};

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

	/* find the minimal interval node */
	AvlNode<long> * & findMin(AvlNode<long> * & tree);
	
	/* find the maximal interval node */
	AvlNode<long> * & findMax(AvlNode<long> * & tree);

	void remove(AvlNode<long> * & tree, std::pair<long, long> & inter);

	void rotate(AvlNode<long> * & tree);

	void doubleRotate(AvlNode<long> * & tree);

	void balance(AvlNode<long> * & tree);

	void calReuseDist(uint64_t addr, Histogram<> & hist);
};

/* useless */
class ListStack
{
	std::list<uint64_t> addrList;

	Histogram<int, double> hist;

public:
	ListStack() {};

	~ListStack() {};

	void calStackDist(uint64_t addr);
};

/* do reuse distance statistics */
class ReuseDist
{
	std::map<uint64_t, long> addrMap;
	long index;

public:
	ReuseDist() {};
	~ReuseDist() {};

	void calReuseDist(uint64_t addr, Histogram<> & hist);
};

/* stack distance staticstics with sampling */
class SampleStack
{
	long index;

	typedef std::unordered_set<uint64_t> AddrSet;
	/* each watchpoint has a set to keep the unique mem refs */
	std::unordered_map<uint64_t, AddrSet> addrTable;
	int randNum;
	int sampleCounter;
	int expectSamples;
	int hibernInter;
	int sampleInter;
	long statusCounter;
	bool isSampling;

public:
	SampleStack() : sampleCounter(0), expectSamples(2000), hibernInter(3000000), 
		            sampleInter(1000000), isSampling(false), statusCounter(0) {};

	/* to generate a random number */
	int genRandom();

	void calStackDist(uint64_t addr, Histogram<> & hist);
};

class Reader
{
private:
#ifdef STACK
	AvlTreeStack avlTreeStack;
#endif

#ifdef REUSE
	ReuseDist reuseDist;
#endif

#ifdef SAMPLE
	SampleStack sampleStack;
#endif
	Histogram<> histogram;

public:
	Reader(std::string _path, std::string _pathout);
};