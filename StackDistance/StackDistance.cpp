// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"

inline double power(const double base, const int index)
{
	if (!index || base <= 0)
		return 0.0;

	double ans = 1.0;
	for (int i = 1; i <= index; ++i) {
		ans *= base;
		/* keep off underflow_error,
		actually we should throw an underflow exception */
		if (ans <= DBL_MIN)
			return 0.0;
	}

	return ans;
}

inline double biDistribution(const int m, const int n, const double p)
{
	if (m < 0 || m > n || p < 0 || p > 1)
		return 0.0;

	double ans = 1.0;
	int i = 1; /* we start from 1. */

	if (m >= n - m)
		while (i <= m) {
			double np = i <= (n - m) ? (1 - p) : 1.0;
			ans *= (n - i + 1) * p * np;
			ans /= i;
			++i;
			/* keep off underflow_error */
			if (ans <= DBL_MIN)
				return 0.0;
		}
	else
		while (i <= n - m) {
			if (i <= m) {
				ans *= (n - i + 1) * p * (1 - p);
				ans /= i;
			}
			else
				ans *= (1 - p);
			++i;
			/* keep off underflow_error */
			if (ans <= DBL_MIN)
				return 0.0;
		}
	return ans;
}

template <class B, class T>
void Histogram<B, T>::sample(B x)
{
	if (!binsVec[x]) {
		binsVec[x] = 1;
	}
	else
		++binsVec[x];
	/* calculate the total num of sampling */
	++samples;
}

//template <class B, class T>
//bool Histogram<B, T>::intoVector()
//{
//	std::map<B, T>::iterator it = bins.begin();
//	std::map<B, T>::iterator itLast = bins.end();
//	--itLast;
//	assert(itLast->first == MISS_BAR);
//
//	for (B i = 0; i <= MISS_BAR; ++i) {
//		/* we keep 0 value bars */
//		if (it->first != i)
//			binsVec.push_back(0);
//		else {
//			binsVec.push_back((double)it->second);
//			++it;
//		}
//		/* initialize transBins */
//		transBins.push_back(0);
//	}
//
//	return binsVec.size() == MISS_BAR + 1;
//}

template <class B, class T>
void Histogram<B, T>::changeAssoc(const int & cap, const int & blk, const int & assoc)
{
	assert(binsVec[MISS_BAR] && binsVec[1]);
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	double p = (double)1 / setNum;

	/* calculate the histogram for each set */
	for (int s = 0; s < setNum; ++s) {
		/* transition of each bar */
		for (int i = 1; i < binsVec.size(); ++i) {
			if (!binsVec[i]) continue;
			/* each bar j before the current, need to be added,
			it's a binominal distribution */
			for (int j = i - 1; j >= 0; --j) {
				transBins[j] += p * binsVec[i] * biDistribution(j, i, p);
			}
			/* update current bar */
			transBins[i] += p * binsVec[i] * power(p, i);
		}
	}
	/* for the original first bar, no need to change */
	transBins[0] += binsVec[0];
}

template <class B, class T>
void Histogram<B, T>::print(std::ofstream & file)
{
	file << "total_samples" << "\t" << samples << std::endl;
	file << "total_misses" << "\t" << misses << std::endl;
	file << "miss_rate" << "\t" << missRate << std::endl;

	for (int i = 0; i < transBins.size(); ++i) {
		/* do not print zero value bins */
		if (!transBins[i]) continue;
		file << "dist" << i << "\t" << transBins[i] << std::endl;
	}
}

void AvlTreeStack::destroy(AvlNode<long> * & tree)
{
	if (!tree)
		return;

	destroy(tree->left);
	destroy(tree->right);
	delete tree;
	tree = nullptr;
}


void AvlTreeStack::insert(AvlNode<long> * & tree, long & v) {
	if (!tree) {
		tree = new AvlNode<long>(v);
		return;
	}

	std::pair<long, long> & interval = tree->interval;
	assert(!(v >= interval.first && v <= interval.second));

	AvlNode<long> * & n1 = tree->left;
	AvlNode<long> * & n2 = tree->right;

	if (tree->holes < 0)
		return;

	if (v == interval.first - 1) {
		dist += interval.second - v + tree->rHoles;
		interval.first = v;
		if (n1) {
			std::pair<long, long> & temp = findMax(n1);
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
			std::pair<long, long> & temp = findMin(n2);
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


void AvlTreeStack::insert(long & a) { dist = 0; insert(root, a); }


std::pair<long, long> & AvlTreeStack::findMin(AvlNode<long> * & tree)
{
	assert(tree);

	if (!tree->left)
		return tree->interval;
	else
		findMin(tree->left);
}


std::pair<long, long> & AvlTreeStack::findMax(AvlNode<long> * & tree)
{
	assert(tree);

	if (!tree->right)
		return tree->interval;
	else
		findMax(tree->right);
}


void AvlTreeStack::remove(AvlNode<long> * & tree, std::pair<long, long> & inter)
{
	if (!tree)
		return;

	if (inter.first > tree->interval.first)
		remove(tree->right, inter);
	else if (inter.first < tree->interval.first)
		remove(tree->left, inter);
	else if (tree->left && tree->right) {
		std::pair<long, long> & temp = findMin(tree->right);
		remove(tree, temp);
	}
	else {
		AvlNode<long> * old = tree;
		if (!tree->left && !tree->right)
			tree = nullptr;
		else
			tree = new AvlNode<long>(*(tree->left ? tree->left : tree->right));

		delete old;
	}

	balance(tree);
}


void AvlTreeStack::rotate(AvlNode<long> * & tree)
{
	if (!tree)
		return;

	if (tree->holes < 0)
		return;

	AvlNode<long> * n1 = tree->left;
	AvlNode<long> * n2 = tree->right;

	if (n2 && AvlNode<long>::getHeight(n1) < AvlNode<long>::getHeight(n2)) {
		tree->right = n2->left;
		tree->updateHoles();
		tree->updateHeight();

		n2->left = tree;
		n2->updateHoles();
		n2->updateHeight();

		tree = n2;
	}
	else if (n1 && AvlNode<long>::getHeight(n1) > AvlNode<long>::getHeight(n2)) {
		tree->left = n1->right;
		tree->updateHoles();
		tree->updateHeight();

		n1->right = tree;
		n1->updateHoles();
		n1->updateHeight();

		tree = n1;
	}
}


void AvlTreeStack::doubleRotate(AvlNode<long> * & tree)
{
	if (!tree)
		return;

	if (tree->holes < 0)
		return;

	AvlNode<long> * & n1 = tree->left;
	AvlNode<long> * & n2 = tree->right;

	if (AvlNode<long>::getHeight(n1) < AvlNode<long>::getHeight(n2))
		rotate(n2);
	else
		rotate(n1);

	rotate(tree);
}


void AvlTreeStack::balance(AvlNode<long> * & tree)
{
	if (!tree)
		return;

	AvlNode<long> * & n1 = tree->left;
	AvlNode<long> * & n2 = tree->right;
	int err = AvlNode<long>::getHeight(n1) - AvlNode<long>::getHeight(n2);
	if (err > 1) {
		if (AvlNode<long>::getHeight(n1->left) >= AvlNode<long>::getHeight(n1->right))
			rotate(tree);
		else
			doubleRotate(tree);
	}
	else if (err < -1) {
		if (AvlNode<long>::getHeight(n2->left) <= AvlNode<long>::getHeight(n2->right))
			rotate(tree);
		else
			doubleRotate(tree);
	}

	tree->updateHeight();
	tree->updateHoles();
}

void AvlTreeStack::calStackDist(uint64_t addr)
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
		/* insert a hole */
		insert(value);
		int stackDist = index - value - dist - 1;
		/* if the stack distance is large than MISS_BAR, the reference is definitely missed. */
		stackDist = stackDist >= MISS_BAR ? MISS_BAR : stackDist;
		hist.sample(stackDist);
	}

	value = index;
}

void AvlTreeStack::transHist(const int cap, const int blk, const int assoc)
{
	hist.changeAssoc(cap, blk, assoc);
}

Reader::Reader(std::string _path, std::string _pathout) {
	std::ifstream file;
	std::ofstream fileOut;

	file.open(_path, std::ios::in);
	if (file.fail()) {
		std::cout << "File openning failed! " << std::endl;
		return;
	}

	fileOut.open(_pathout, std::ios::out);
	if (fileOut.fail()) {
		std::cout << "File openning failed! " << std::endl;
		return;
	}

	std::string temp;
	std::string line;

	std::cout << "Reading files and calculate stack distances..." << std::endl;
	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		uint64_t paddr;
		uint64_t vaddr;
		lineStream >> temp;
		lineStream >> std::hex >> paddr;
		lineStream >> temp;
		lineStream >> std::hex >> vaddr;

		uint64_t addr = vaddr;
		/* check a valid address */
		if (!addr)
			continue;
		
		/* call the methods to calculate stack distances */
		//listStack.calStackDist(addr);
		uint64_t mask = 63;
		addr = addr & (~mask);
		//avlTreeStack.calStackDist(addr);
		reuseDist.calReuseDist(addr);
	}
	//listStack.print("E:\\ShareShen\\gem5-origin\\m5out-se-x86\\perlbench.txt");
	std::cout << "transiting stack distances..." << std::endl;

	int cap = 32 * 1024;
	int blk = 64;
	int assoc = 4;
	//avlTreeStack.transHist(cap, blk, assoc);
	//avlTreeStack.hist.calMissRate(assoc);
	//avlTreeStack.hist.print(fileOut);
	reuseDist.transHist();
	reuseDist.hist.calMissRate(cap, blk);
	file.close();
	fileOut.close();
}

int main()
{
	//Histogram<int, double> hist;
	//hist.sample(1);
	//hist.sample(4);
	//hist.sample(4);
	//Histogram<int, double> hist("E:\\ShareShen\\gem5-origin\\m5out-se-x86\\perlbench-l1d32k256assoc-set-distribution.txt");
	Reader reader("E:\\ShareShen\\gem5-stable\\m5out-se-x86\\requtTraceFile-perlbench.txt", \
		"E:\\ShareShen\\gem5-stable\\m5out-se-x86\\perlbench-avl-4assoc.txt");

 	return 0;
}