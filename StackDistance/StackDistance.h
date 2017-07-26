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
//#define LOG
static int64_t MISS_BAR = 0;

static int64_t Trunc = 0;

/* for recording distribution into a Histogram, 
   Accur is the accuracy of transforming calculation */
template <class B = int64_t, class Accur = double>
class Histogram
{
	/* Histogram implemented by std::map */
	std::map<long, B> binsMap;
	/* Histogram implemented by std::vector */
	//std::vector<B> binsVec;
	std::vector<Accur> binsVec;
	/* Histogram implemented by std::vector, 
	   after fully-to-set-associative cache transformation */
	std::vector<Accur> binsTra;
	/* number of sampling */
	B samples;
	/* total number of hit references */
	B hits;
	/* cache miss rate */
	Accur missRate;

public:
	Histogram() : samples(0), hits(0), missRate(0.0) {};

	~Histogram() {};

	void clear();

	void sample(B x);

	bool mapToVector();
	/* complete the histogram with std::vector, to fast calculation. */
	bool mapToVector(B * buffer, int bufSize);

	bool mapToVector(std::vector<B> & buffer);
	/* use boost binomial distribution to do a fully-to-set-associative cache transformation */
	Accur fullyToSetAssoc(const int & cap, const int & blk, const int & assoc);
	/* use boost poisson distribution to do a fully-to-set-associative cache transformation */
	Accur calLruMissRatePoisson(const int & cap, const int & blk, const int & assoc);
	/* transform reuse distance distribution to stack distance distribution */
	void reuseDistToStackDist();

	void print(std::ofstream & file);
	/* calculate the miss rate for LRU set associative cache via cdf of binomial distribution */
	Accur calLruMissRate(const int & cap, const int & blk, const int & assoc);
	
	Accur calPlruMissRate(const int & cap, const int & blk, const int & assoc);
	/* calculate the miss rate for PLRU set associative cache, 
	   first do set-association transforming, then calculate hit function */
	Accur calMissRate(const int & cap, const int & blk, const int & assoc, const bool plru);
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
	AvlNode(Accur & a);

	AvlNode(AvlNode<Accur> & n);

	~AvlNode();

	static int getHeight(AvlNode<Accur> * & node)
	{
		return node ? node->height : -1;
	}

	void updateHeight();

	void updateHoles();
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
	AvlTreeStack(long & v);

	AvlTreeStack();

	~AvlTreeStack() { destroy(root); }

	void destroy(AvlNode<long> * & tree);

	void clear();

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

	void calStackDist(uint64_t addr, Histogram<> & hist);
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
	/* hibernating interval size */
	int hibernInter;
	/* sampling interval size */
	int sampleInter;
	/* the counter starts a sampling or hibernating interval */
	uint64_t statusCounter;
	/* the counter chooses the random sample to record SD */
	uint64_t sampleCounter;
	enum State {sampling, hibernating};
	State state;
	double sampleRate;

public:
	SampleStack(double);

	SampleStack(int sSize, int hSize, double sRate);

	/* to generate a random number */
	int genRandom();

	void calStackDist(uint64_t addr, Histogram<> & hist);
};

class Reader
{
public:
	Reader(std::ifstream & fin);
};