// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"

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
bool Histogram<B, T>::mapToVector()
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
void Histogram<B, T>::reuseDistToStackDisth()
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
void Histogram<B, T>::fullyToNWay(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	double p = (double)1 / setNum;

	binsTra.clear();

	/* if a mem ref's stack distance <= assoc - 1, it is a cache hit */
	for (int i = 0; i < assoc; ++i)
		binsTra.push_back(binsVec[i]);

	binsTra.resize(binsVec.size());

	/* transition of each bar */
	for (int i = assoc; i < binsVec.size(); ++i) {
		if (!binsVec[i]) continue;

		/* start from j=i, because i is the current bar */
		for (int j = i; j >= 0; --j)
			binsTra[j] += binsVec[i] * biDistribution(j, i, p);
	}
}

template <class B, class T>
void Histogram<B, T>::fullyToNWay2(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	double p = (double)1 / setNum;

	binsTra.clear();

	/* if a mem ref's stack distance <= assoc - 1, it is a cache hit */
	for (int i = 0; i < assoc; ++i)
		binsTra.push_back(binsVec[i]);
	
	binsTra.resize(binsVec.size());

	try {
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
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}
}

template <class B, class T>
void Histogram<B, T>::fullyToNWay3(const int & cap, const int & blk, const int & assoc)
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
	misses = 0;

	for (int i = assoc; i < binsTra.size(); ++i)
		misses += (B)std::round(binsTra[i]);
	missRate = (double)misses / samples;
}

template <class B, class T>
void Histogram<B, T>::calMissRate(const int & cap, const int & blk)
{
	int blkNum = cap / blk;
	for (int i = blkNum; i < binsVec.size(); ++i)
		misses += binsVec[i];
	missRate = (double)misses / samples;
}

/* to calculate C(a, k) / C(b, k) */
inline double combinationRatio(int b, int a, int k)
{
	if (b < a || !a || !b || !k)
		throw std::exception("Error: combinationRatio: wrong arguments");
	
	/* This is for the condition that wouldn't occur. */
	if (k > a) return 0.0;


	int loops = b - a;
	double result = 1.0;
	while (loops--)
		result *= (double)(b - k) / (b--);

	if (result <= DBL_MIN)
		throw std::exception("Error: combinationRatio: underflow");

	return result;
}

inline double serialAccesPerm()
{
}

template <class B, class T>
void Histogram<B, T>::calMissRate(const int & assoc, const bool plru)
{
	if (!plru)
		return;

	int tempAssoc = assoc - 1;
	/* secureDist = log2(assoc),
	   the least not evcited distance, called secure distance */
	int secureDist = 0;
	while (tempAssoc) {
		tempAssoc >>= 1;
		++secureDist;
	}

	/* Now we calculate the approximate miss rate in an assumed LRU cache
	   that the number of blocks equals to association,
	   and references with stack distance <= log2(assoc) must be hit. */
	double tempMisses = 0.0;
	double tempHits = 0.0;
	for (int i = 0; i <= secureDist; ++i)
		tempHits += binsTra[i];
	for (int i = secureDist + 1; i < assoc; ++i)
		tempMisses += binsTra[i];
	double lastAccesMiss = tempMisses / (tempMisses + tempHits);

	double log = secureDist;
	/* probability for accesses being in-order
	   no need to use BOOST */
	double inOrderAccesP = 1.0;
	while (log)
		inOrderAccesP *= log--;
	inOrderAccesP = 1 / inOrderAccesP;

	double missesD = 0.0;

	try {
		for (int i = secureDist + 1; i < assoc; ++i) {
			/* probability for eviction */
			double p = 1.0;
			/* probability for accessing in-order */
			double pInOrder = 1.0;
			/* accesses before current cousins */
			double front = 0.0;
			for (int j = 0; j < secureDist; ++j) {
				double c = combinationRatio(assoc - 1, assoc - 1 - (1 << j), i - 1);
				p *= 1 - c;

				//pInOrder *= 1 - std::pow((double)j / (j + 1), (double)(1 << j) * (i - 1) / (assoc - 1));
				
				/* number of accesses in the serialization in same subtree */
				double cousins = std::max(std::floor((1 << j) * (i - 1) / (assoc - 1)), 1.0);
				/* tgammma_ratio(a, b) for (a-1)! / (b-1)! */
				pInOrder *= cousins * boost::math::tgamma_ratio(front + cousins, front + 1);

				front += cousins;
				if (pInOrder < DBL_MIN)
					throw std::exception("Error: calMissRate: underflow");
			}

			pInOrder /= boost::math::tgamma(i);
			if (pInOrder >= 1)
				throw std::exception("Error: calMissRate: wrong calculation");
			missesD += binsTra[i] * p * pInOrder * lastAccesMiss;
		}
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
		return;
	}

	misses = std::ceil(missesD);
	/* references with SD > assoc always miss. */
	for (int i = assoc; i < binsVec.size(); ++i)
		misses += binsTra[i];

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
		interval.first = v;
		if (n1) {
			std::pair<long, long> & temp = findMax(n1)->interval;
			if (v == temp.second + 1) {
				interval.first = temp.first;
				remove(n1, temp);
			}
		}
		/* update the holes between current address */
		curHoles += interval.second - v + tree->rHoles;
	}

	else if (v == interval.second + 1) {
		interval.second = v;
		if (n2) {
			std::pair<long, long> & temp = findMin(n2)->interval;
			if (v == temp.first - 1) {
				interval.second = temp.second;
				remove(n2, temp);
			}
		}
		curHoles += tree->rHoles;
	}

	else if (v < interval.first - 1) {
		insert(n1, v);
		curHoles += tree->rHoles + interval.second - interval.first + 1;
	}

	else if (v > tree->interval.second + 1)
		insert(n2, v);

	balance(tree);
}

void AvlTreeStack::insert(long & a) 
{	
	/* everytime insert a address, get its total holes */
	curHoles = 0; 
	insert(root, a);
}

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

void AvlTreeStack::insert(uint64_t addr)
{
	long & value = addrMap[addr];

	++index;

	/* value is 0 under cold miss */
	if (!value) {
#ifdef AVL_HIST
		hist.sample(MISS_BAR);
#endif
		value = index;
		return;
	}

	/* update b of last reference */
	if (value < index) {
		/* insert a hole */
		insert(value);
		int stackDist = index - value - curHoles - 1;
		/* if the stack distance is large than MISS_BAR, the reference is definitely missed. */
		stackDist = stackDist >= MISS_BAR ? MISS_BAR : stackDist;
#ifdef AVL_HIST
		hist.sample(stackDist);
#endif
	}

	value = index;
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
	int interval = 0;
	bool start = false;

	std::cout << "Reading files and calculate stack distances..." << std::endl;

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		uint64_t paddr;
		uint64_t vaddr;
		lineStream >> temp;
		
		if (temp == "interval") {
			++interval;
			if (interval == 2) {
				start = true;
				continue;
			}
			else if (interval > 2)
				break;
			else
				continue;
		}

		if (!start)
			continue;

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
		avlTreeStack.insert(addr);
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

	int cap = 32 * 1024;
	int blk = 64;
	int assoc = 16;
#ifdef STACK
	bool succ = avlTreeStack.hist.mapToVector();
	avlTreeStack.hist.fullyToNWay2(cap, blk, assoc);
	avlTreeStack.hist.calMissRate(assoc);
	avlTreeStack.hist.calMissRate(assoc, true);
#endif

#ifdef REUSE
	reuseDist.transHist();
	reuseDist.hist.fullyToNWay2(cap, blk, assoc);
	reuseDist.hist.calMissRate(assoc);
	//reuseDist.hist.calMissRate(cap, blk);
#endif

#ifdef SAMPLE
	bool succ = sampleStack.hist.mapToVector();
	clock_t starts = clock();
	sampleStack.hist.fullyToNWay(cap, blk, assoc);
	double t = clock() - starts;
	sampleStack.hist.calMissRate(assoc);

	start = clock();
	sampleStack.hist.fullyToNWay2(cap, blk, assoc);
	double t2 = clock() - starts;
	sampleStack.hist.calMissRate(assoc);
#endif // SAMPLE

	file.close();
	fileOut.close();
}

int main()
{
	//double c = boost::math::tgamma((double)50);
	double c = boost::math::tgamma_ratio(1, 1);
	Reader reader("E:\\ShareShen\\gem5-stable\\m5out-se-x86\\cactusADM\\cactusADM-trace-part.txt", \
		"E:\\ShareShen\\gem5-stable\\m5out-se-x86\\cactusADM\\cactusADM-avl-2assoc.txt");

 	return 0;
}