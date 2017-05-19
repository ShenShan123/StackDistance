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

//#define REUSE
//#define STACK
#define SAMPLE

#define MISS_BAR 1024 * 2

inline double power(const double base,const int index, double coef = 1.0);

inline double biDistribution(const int m, const int n, const double p);

template <class B, class T>
class Histogram
{
	std::map<long, B> binsMap;
	std::vector<B> binsVec;
	std::vector<T> binsTra;
	/* number of sampling */
	B samples;
	/* total number of miss references */
	B misses;
	/* miss rate */
	double missRate;

public:
	Histogram() : misses(0), missRate(0.0)
		//binsVec(MISS_BAR + 1, 0), 
		//binsTra(MISS_BAR + 1, 0) 
	{};

	~Histogram() {};

	void sample(B x);

	bool intoVector();
	/* use self-defined binomial fucntion */
	void changeAssoc(const int & cap, const int & blk, const int & assoc);
	/* use boost binomial distribution */
	void changeAssoc2(const int & cap, const int & blk, const int & assoc);
	/* use boost poisson distribution */
	void changeAssoc3(const int & cap, const int & blk, const int & assoc);

	void transToStackDist();

	void print(std::ofstream & file);

	void calMissRate(const int & assoc);

	void calMissRate(const int & cap, const int & blk)
	{
		int blkNum = cap / blk;
		for (int i = blkNum; i < binsVec.size(); ++i)
			misses += binsVec[i];
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
	std::map <uint64_t, long> addrMap;
	long index;
	int dist;

public:
	AvlNode<long> * root;

	/* hist contains the stack distance distribution */
	Histogram<int, double> hist;

	AvlTreeStack(long & v) : dist(0)
	{
		root = new AvlNode<long>(v);
	}

	AvlTreeStack() : root(nullptr), dist(0) {};

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

	AvlNode<long> * & findMin(AvlNode<long> * & tree);

	AvlNode<long> * & findMax(AvlNode<long> * & tree);

	void remove(AvlNode<long> * & tree, std::pair<long, long> & inter);

	void rotate(AvlNode<long> * & tree);

	void doubleRotate(AvlNode<long> * & tree);

	void balance(AvlNode<long> * & tree);

	void calStackDist(uint64_t addr);

	void smapledStackDist(uint64_t addr);

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

class SampleStack
{
	long index;
	std::unordered_map<uint64_t, std::unordered_set<uint64_t>> addrTable;
	int randNum;
	int sampleCounter;
	int expectSamples;
	int hibernInter;
	int sampleInter;
	long statusCounter;
	bool isSampling;

public:
	Histogram<int, double> hist;

	SampleStack() : sampleCounter(0), expectSamples(2000), hibernInter(3000000), 
		            sampleInter(1000000), isSampling(false), statusCounter(0) {};

	int genRandom()
	{
		// construct a trivial random generator engine from a time-based seed:
		int seed = std::chrono::system_clock::now().time_since_epoch().count();
		static std::mt19937 engine(seed);
		static std::geometric_distribution<int> geom((double)expectSamples / sampleInter);
		int rand = 0;
		/* the random can not be 0 */
		while (!rand)
			rand = geom(engine);

		return rand;
	}

	void calStackDist(uint64_t addr)
	{
		++sampleCounter;
		++statusCounter;

		/* start sampling interval */
		if (!isSampling && statusCounter == hibernInter) {
			statusCounter = 0;
			sampleCounter = 0;
			isSampling = true;
			randNum = genRandom();
		}
		/* start hibernation interval */
		else if (isSampling && statusCounter == sampleInter) {
			statusCounter = 0;
			isSampling = false;
		}

		/* if we find a same address x in addrTable,
		   record its stack distance and
		   the sampling of x is finished */
		std::unordered_map<uint64_t, std::unordered_set<uint64_t>>::iterator pos;
		pos = addrTable.find(addr);
		if (pos != addrTable.end()) {
			int dist = pos->second.size() - 1;
			hist.sample(dist);
			addrTable.erase(addr);
		}

		std::unordered_map<uint64_t, std::unordered_set<uint64_t>>::iterator it = addrTable.begin();
		
		/* make sure the max size of addrTable */
		if (addrTable.size() > 500) {
			hist.sample(MISS_BAR);
			addrTable.erase(it->first);
		}

		/* record unique mem references between the sampled address x */
		for (it = addrTable.begin(); it != addrTable.end(); ) {
			it->second.insert(addr);
			/* if the set of sampled address x is too large,
			   erase it from  the table and record as MISS_BAR */
			if (it->second.size() > MISS_BAR) {
				std::unordered_map<uint64_t, std::unordered_set<uint64_t>>::iterator eraseIt = it;
				++it;
				hist.sample(MISS_BAR);
				addrTable.erase(eraseIt->first);
			}
			else
				++it;
		}

		/* if it is time to do sampling */
		if (isSampling && sampleCounter == randNum) {
			/* it is a new sampled address */
			assert(addrTable[addr].empty());
			addrTable[addr].insert(0);
			/* reset the sampleCounter and randNum to prepare next sample */
			sampleCounter = 0;
			randNum = genRandom();
			//randNum = 1;
		}
	}
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

	ListStack listStack;

public:
	Reader(std::string _path, std::string _pathout);
};