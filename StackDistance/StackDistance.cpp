// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"

inline double power(const double base, const int pow, double coef)
{
	if (base <= 0)
		return 0.0;

	for (int i = 1; i <= pow; ++i) {
		coef *= base;
		/* keep off underflow_error,
		actually we should throw an underflow exception */
		if (coef <= DBL_MIN)
			return 0.0;
	}

	return coef;
}

inline double biDistribution(const int m, const int n, const double p)
{
	if (m < 0 || m > n || p < 0 || p > 1)
		return 0.0;
	else if (m == 0)
		return 1.0;

	double ans = 1.0;
	//int i = 1; /* we start from 1. */
	bool hasP = true;

	if (m >= n - m) {
		int i = n - m;
		while (i > 0) {
			ans *= (n - i + 1) * (1 - p) * p;
			ans /= i;
			--i;
			/* keep off underflow_error */
			if (ans <= DBL_MIN)
				return 0.0;
			else if (ans >= FLT_MAX && hasP) {
				ans = power(p, m - n + m, ans);
				hasP = false;
			}

			assert(ans < DBL_MAX);
		}
		if (hasP)
			ans = power(p, m - n + m, ans);
	}
	else {
		int i = m;
		while (i > 0) {
			ans *= (n - i + 1) * p * (1 - p);
			ans /= i;
			--i;
			/* keep off underflow_error */
			if (ans <= DBL_MIN)
				return 0.0;
			else if (ans >= FLT_MAX && hasP) {
				ans = power(1 - p, n - m - m, ans);
				hasP = false;
			}
			assert(ans < DBL_MAX);
		}
		if (hasP)
			ans = power(1 - p, n - m - m, ans);
	}
	return ans;
}

template <class B, class T>
void Histogram<B, T>::sample(B x)
{
	if (!binsMap[x]){ 
		binsMap[x] = 1;
	}
	else
		++binsMap[x];
	/* calculate the total num of sampling */
	++samples;
}

template <class B, class T>
bool Histogram<B, T>::intoVector()
{
	std::map<long, B>::iterator last = binsMap.end();
	/* point to the last element */
	long maxSize = (--last)->first;

	binsVec.reserve(maxSize);

	long temp = 0;
	std::map<long, B>::iterator it = binsMap.begin();

	for (int i = 0; i <= maxSize; ++i)
		if (temp == it->first) {
			binsVec.push_back(it->second);
			++it;
			++temp;
		}
		else {
			binsVec.push_back(0);
			++temp;
		}

		return binsVec.size() == last->first + 1;
}

template <class B, class T>
void Histogram<B, T>::transToStackDist()
{
	std::vector<double> frac;
	B temp = 0;
	/* calculate the fraction of reuse distance
	is greater than i */
	for (int i = 0; i < binsVec.size(); ++i) {
		temp += binsVec[i];
		frac.push_back((double)(samples - temp) / samples);
	}

	std::vector<T> transTemp(binsTra);

	for (int i = 1; i < binsVec.size(); ++i) {
		double sumFrac = 0.0;
		for (int j = 1; j <= i; ++j)
			sumFrac += frac[j];
		transTemp[(int)std::round(sumFrac)] += binsVec[i];
	}
	transTemp[0] = binsVec[0];

	binsVec.clear();
	binsVec = transTemp;
}

template <class B, class T>
void Histogram<B, T>::changeAssoc(const int & cap, const int & blk, const int & assoc)
{
	binsTra.resize(binsVec.size(), 0);
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	double p = (double)1 / setNum;

	/* calculate the histogram for each set */
	for (int s = 0; s < setNum; ++s) {
		/* transition of each bar */
		for (int i = 1; i < binsVec.size(); ++i) {
			if (!binsVec[i]) continue;

			/* each bar j before the current, need to be added,
			it's a binomial distribution */
			boost::math::binomial binom(i, p);
			
			for (int j = i - 1; j >= 0; --j)
				//binsTra[j] += p * binsVec[i] * biDistribution(j, i, p);
				binsTra[j] += p * binsVec[i] * boost::math::pdf(binom, j);

			/* update current bar */
			binsTra[i] += p * binsVec[i] * power(p, i);
		}
	}
	/* for the original first bar, no need to change */
	binsTra[0] += binsVec[0];
}

template <class B, class T>
void Histogram<B, T>::changeAssoc2(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	double p = (double)1 / setNum;

	/* if a mem ref's stack distance <= assoc - 1, it is a cache hit */
	for (int i = 0; i < assoc; ++i)
		binsTra.push_back(binsVec[i]);
	
	binsTra.resize(binsVec.size());

	/* transition of each bar */
	for (int i = assoc; i < binsVec.size(); ++i) {
		if (!binsVec[i]) continue;

		/* each bar j before the current, need to be added,
		it's a binomial distribution */
		boost::math::binomial binom(i, p);

		/* start from j=i, because i is the current bar */
		for (int j = i; j >= 0; --j)
			binsTra[j] += binsVec[i] * boost::math::pdf(binom, j);
	}
}

template <class B, class T>
void Histogram<B, T>::changeAssoc3(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	double p = (double)1 / setNum;

	/* if a mem ref's stack distance <= assoc - 1, it is a cache hit */
	for (int i = 0; i < assoc; ++i)
		binsTra.push_back(binsVec[i]);

	binsTra.resize(binsVec.size());

	/* transition of each bar */
	for (int i = assoc; i < binsVec.size(); ++i) {
		if (!binsVec[i]) continue;

		/* each bar j before the current, need to be added,
		it's a binomial distribution */
		boost::math::poisson_distribution<double> pois(i * p);

		/* start from j=i, because i is the current bar */
		for (int j = i; j >= 0; --j)
			binsTra[j] += binsVec[i] * boost::math::pdf(pois, j);
	}
}

template <class B, class T>
void Histogram<B, T>::calMissRate(const int & assoc)
{
	for (int i = assoc; i < binsTra.size(); ++i)
		misses += (B)std::round(binsTra[i]);
	missRate = (double)misses / samples;
}

template <class B, class T>
void Histogram<B, T>::print(std::ofstream & file)
{
	file << "total_samples" << "\t" << samples << std::endl;
	file << "total_misses" << "\t" << misses << std::endl;
	file << "miss_rate" << "\t" << missRate << std::endl;

	for (int i = 0; i < binsTra.size(); ++i) {
		/* do not print zero value bins */
		if (!binsTra[i]) continue;
		file << "dist" << i << "\t" << binsTra[i] << std::endl;
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
			std::pair<long, long> & temp = findMax(n1)->interval;
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
			std::pair<long, long> & temp = findMin(n2)->interval;
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

AvlNode<long> * & AvlTreeStack::findMin(AvlNode<long> * & tree)
{
	assert(tree);

	if (!tree->left)
		return tree;
	else
		findMin(tree->left);
}

AvlNode<long> * & AvlTreeStack::findMax(AvlNode<long> * & tree)
{
	assert(tree);

	if (!tree->right)
		return tree;
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

	/* the tree has two children , 
	   replace its content with the min subnode and delete the min node */
	else if (tree->left && tree->right) {
		AvlNode<long> * & minNode = findMin(tree->right);
		tree->interval = minNode->interval;
		remove(tree->right, tree->interval);
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

void AvlTreeStack::smapledStackDist(uint64_t addr)
{
	++index;
}

void AvlTreeStack::transHist(const int cap, const int blk, const int assoc)
{
	//hist.changeAssoc(cap, blk, assoc);
	hist.changeAssoc2(cap, blk, assoc);
	//hist.changeAssoc3(cap, blk, assoc);
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
	int interval = 1;

	std::cout << "Reading files and calculate stack distances..." << std::endl;

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		uint64_t paddr;
		uint64_t vaddr;
		lineStream >> temp;
		
		if (temp == "interval") {
			if (interval)
				--interval;
			else
				break;
		}

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

#ifdef STACK
		avlTreeStack.calStackDist(addr);
#endif

#ifdef REUSE
		reuseDist.calReuseDist(addr);
#endif

#ifdef SAMPLE
		sampleStack.calStackDist(addr);
#endif
	}
	//listStack.print("E:\\ShareShen\\gem5-origin\\m5out-se-x86\\perlbench.txt");
	std::cout << "transiting stack distances..." << std::endl;

	int cap = 64 * 1024;
	int blk = 64;
	int assoc = 2;
#ifdef STACK
	bool succ = avlTreeStack.hist.intoVector();
	avlTreeStack.transHist(cap, blk, assoc);
	avlTreeStack.hist.calMissRate(assoc);
	//avlTreeStack.hist.print(fileOut);
#endif

#ifdef REUSE
	reuseDist.transHist();
	reuseDist.hist.changeAssoc(cap, blk, assoc);
	reuseDist.hist.calMissRate(assoc);
	//reuseDist.hist.calMissRate(cap, blk);
#endif

#ifdef SAMPLE
	bool succ = sampleStack.hist.intoVector();
	sampleStack.hist.changeAssoc2(cap, blk, assoc);
	sampleStack.hist.calMissRate(assoc);
#endif // SAMPLE

	file.close();
	fileOut.close();
}

int main()
{
	/*SampleStack sampleStack;
	sampleStack.calStackDist(8);
	sampleStack.calStackDist(20);
	sampleStack.calStackDist(21);
	sampleStack.calStackDist(21);
	sampleStack.calStackDist(23);
	sampleStack.calStackDist(24);
	sampleStack.calStackDist(20);
	sampleStack.calStackDist(8);*/

	AvlTreeStack avlTreeStack;
	long temp;
	temp = 4;
	avlTreeStack.insert(temp);
	temp = 1;
	avlTreeStack.insert(temp);
	temp = 7;
	avlTreeStack.insert(temp);
	temp = 4;
	std::pair<long, long> d = std::make_pair(4, 4);
	avlTreeStack.remove(avlTreeStack.root, d);

	Reader reader("E:\\ShareShen\\gem5-stable\\m5out-se-x86\\cactusADM\\cactusADM-trace-part.txt", \
		"E:\\ShareShen\\gem5-stable\\m5out-se-x86\\cactusADM\\cactusADM-avl-2assoc.txt");

 	return 0;
}