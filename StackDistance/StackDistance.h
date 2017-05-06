#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cmath>

#define MISS_BAR 513

class IntervalSetDistribution
{
public:
	/* setDistr are set distributions of all interval */
	std::vector<std::vector<double>> intSetDistr;
	IntervalSetDistribution() {};

	~IntervalSetDistribution() {};

	void recordSetDistr(const std::string _path)
	{
		std::ifstream file;
		file.open(_path, std::ios::in);

		if (file.fail()) {
			return;
		}

		std::string line;
		std::string temp;
		int samples;
		int numSamples = 0;
		int setIdx = -1;
		int i;
		int j = -1;

		while (std::getline(file, line)) {
			std::stringstream lineStream(line);
			lineStream >> temp;
			std::string idxStr = temp.substr(28);
			/* this is a set distribution of a new window */
			if (idxStr == "samples") {
				++j;
				i = -1;
				lineStream >> numSamples;
				std::vector<double> newSetDistr;
				intSetDistr.push_back(newSetDistr);
			}
			else {
				++i;
				int setIdx = std::stoi(idxStr);
				lineStream >> samples;
				if (setIdx == i)
					intSetDistr[j].push_back((double)samples / numSamples);
				else
					intSetDistr[j].push_back(0);
			}
		}

		file.close();
	}
};

template <class U, class T>
class Histogram
{
	std::map<U, T> origBins;
	std::vector<double> bins;

public:
	std::vector<double> transBins;

	Histogram() {};

	~Histogram() {};

	void clear()
	{
		origBins.clear();
		bins.clear();
		transBins.clear();
	}

	/* do sample, record every time if the distance is x */
	void sample(U x)
	{
		if (!origBins[x]) {
			origBins[x] = 1;
		}
		else
			++origBins[x];
	}

	bool intoVector()
	{
		std::map<U, T>::iterator it = origBins.begin();
		std::map<U, T>::iterator itLast = origBins.end();
		--itLast;
		assert(itLast->first == MISS_BAR);

		for (U i = 0; i <= itLast->first; ++i) {
			/* we keep 0 value bars */
			if (it->first != i)
				bins.push_back(0);
			else {
				bins.push_back((double)it->second);
				++it;
			}
			/* initialize transBins */
			transBins.push_back(0);
		}

		return bins.size() == MISS_BAR + 1;
	}

	void print(std::ofstream * & file)
	{
		for (int i = 1; i < bins.size(); ++i) {
			/* do not print zero value origBins */
			if (!transBins[i]) continue;

			*file << "dist" << i << "\t\t" << transBins[i] << std::endl;
		}
	}

	double combinations(const int m, const int n) {
		if (m < 0 || m > n) {
			return 0;
		}

		double iComb = 1;
		int i = 0;

		while (i < m) {
			++i;
			iComb *= n - i + 1;
			iComb /= i;
		}

		return iComb;
	}

	void changeAssoc(const std::vector<double> & setDistr)
	{
		assert(bins[MISS_BAR] && bins[1]);
		
		/* calculate the histogram for each set */
		for (int s = 0; s < setDistr.size(); ++s) {
			/* p is a probability of mapping to a specific set */
			const double & p = setDistr[s];

			/* transition of each bar */
			for (int i = 2; i < bins.size(); ++i) {
				if (!bins[i]) continue;
				/* each bar j before the current, need to be added,
				   it's a binominal distribution */
				for (int j = i - 1; j > 0; --j) {
					transBins[j] += p * bins[i] * combinations(j - 1, i - 1) *
						std::pow(p, j - 1) * std::pow(1 - p, i - j);
				}
				/* update current bar */
				transBins[i] += bins[i] * std::pow(p, i);
			}
		}
		/* for the original first bar, no need to change */
		transBins[1] += bins[1];
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
	std::map <uint64_t, long> addrMap;
	long index;
	int dist;

public:
	Histogram<int, int> hist;

	AvlTreeStack(long & v) : dist(0), index(0)
	{
		root = new AvlNode<long>(v);
	}

	AvlTreeStack() : root(nullptr), dist(0), index(0) {};

	void setRoot(long & v)
	{
		root = new AvlNode<long>(v);
	}

	~AvlTreeStack() { destroy(root); }

	void clear() 
	{
		destroy(root);
		addrMap.clear();
		hist.clear();
		index = 0;
		dist = 0;
	}

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
		++index;

		long & value = addrMap[addr];

		/* value is 0 under cold miss */
		if (!value) {
			hist.sample(MISS_BAR);
			value = index;
			return;
		}

		/* update b of last reference */
		if (value < index) {
			/* insert a hole */
			insert(value);
			int stackDist = index - value - dist;
			stackDist = stackDist >= MISS_BAR ? MISS_BAR : stackDist;
			hist.sample(stackDist);
		}

		value = index;
	}

	void transStack(const std::vector<double> & setDistr)
	{
		if (hist.intoVector())
			hist.changeAssoc(setDistr);
		else
			std::cout << "fail to intoVector()!!" << std::endl;
	}
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
	std::ofstream outFile;
	std::string line;
	std::vector<Histogram<int, int>> intHists;
	ListStack listStack;
	IntervalSetDistribution intSetDistr;

public:
	Reader(const std::string _path, const std::string _pathout, const std::string _pathset);

	void sumHists();
};