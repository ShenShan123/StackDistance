#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cmath>

#include <unordered_map>

#define MISS_BAR 513

template <class U, class T>
class Histogram
{
	std::map<U, T> bins;
	std::vector<double> binsVec;
	std::vector<double> setDistr;
	std::vector<double> transBins;

public:
	Histogram() {};

	void recordSetDistr(std::string & _path)
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
		int i = -1;

		while (std::getline(file, line)) {
			std::stringstream lineStream(line);
			lineStream >> temp;
			std::string idxStr = temp.substr(28);
			if (idxStr == "samples") {
				lineStream >> numSamples;
			}
			else {
				int setIdx = std::stoi(idxStr);
				lineStream >> samples;
				if (setIdx == i)
					setDistr.push_back((double)samples / numSamples);
				else
					setDistr.push_back(0);
			}
			++i;
		}

		file.close();
	}

	~Histogram() {};

	void sample(U x)
	{
		if (!bins[x]) {
			bins[x] = 1;
		}
		else
			++bins[x];
	}

	bool intoVector()
	{
		std::map<U, T>::iterator it = bins.begin();
		std::map<U, T>::iterator itLast = bins.end();
		--itLast;
		assert(itLast->first == MISS_BAR);

		for (U i = 0; i <= itLast->first; ++i) {
			/* we keep 0 value bars */
			if (it->first != i)
				binsVec.push_back(0);
			else {
				binsVec.push_back((double)it->second);
				++it;
			}
			/* initialize transBins */
			transBins.push_back(0);
		}

		return binsVec.size() == MISS_BAR + 1;
	}

	void print(std::ofstream * & file)
	{
		for (int i = 1; i < binsVec.size(); ++i) {
			/* do not print zero value bins */
			if (!transBins[i]) continue;

			*file << "dist" << i << "\t\t" << transBins[i] << std::endl;
		}
	}

	double combinations(int m, int n) {
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

	void changeAssoc()
	{
		assert(binsVec[MISS_BAR] && binsVec[1]);
		
		/* calculate the histogram for each set */
		for (int s = 0; s < setDistr.size(); ++s) {
			/* p is a probability of mapping to a specific set */
			double & p = setDistr[s];

			/* transition of each bar */
			for (int i = 2; i < binsVec.size(); ++i) {
				if (!binsVec[i]) continue;
				/* each bar j before the current, need to be added,
				   it's a binominal distribution */
				for (int j = i - 1; j > 0; --j) {
					transBins[j] += p * binsVec[i] * combinations(j - 1, i - 1) *
						std::pow(p, j - 1) * std::pow(1 - p, i - j);
				}
				/* update current bar */
				transBins[i] += binsVec[i] * std::pow(p, i);
			}
		}
		/* for the original first bar, no need to change */
		transBins[1] += binsVec[1];
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
	Histogram<int, int> hist;
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

	void setDistr(std::string _path)
	{
		hist.recordSetDistr(_path);
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

	void print(std::string & _path)
	{
		std::ofstream * file = new std::ofstream(_path);

		if (hist.intoVector())
			hist.changeAssoc();

		hist.print(file);
		file->close();
		delete file;
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
	std::string line;
	AvlTreeStack avlTreeStack;
	ListStack listStack;

public:
	Reader(std::string _path, std::string _pathout, std::string _pathset);
};