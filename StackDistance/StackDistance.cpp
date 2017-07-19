// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"
#include "NtlTest.h"

template <class B, class Accur>
void Histogram<B, Accur>::clear()
{
	binsMap.clear();
	binsVec.clear();
	binsTra.clear();
	samples = 0;
	misses = 0;
	hits = 0;
	missRate = 0;
}

template <class B, class Accur>
void Histogram<B, Accur>::sample(B x)
{
	if (!binsMap[x])
		binsMap[x] = 1;
	else
		++binsMap[x];
	/* calculate the total num of sampling */
	++samples;
}

//template <class B, class Accur>
//bool Histogram<B, Accur>::mapToVector()
//{
//	auto last = binsMap.end();
//	/* point to the last element */
//	long maxSize = (--last)->first;
//	/* reserve some memory */
//	binsVec.reserve(maxSize);
//
//	auto it = binsMap.begin();
//
//	for (long i = 0; i <= maxSize; ++i)
//		if (i == it->first) {
//			binsVec.push_back(it->second);
//			++it;
//		}
//		else
//			binsVec.push_back(0);
//
//		return binsVec.size() == last->first + 1;
//}


template <class B, class Accur>
bool Histogram<B, Accur>::mapToVector(B * buffer, int bufSize)
{
	assert(bufSize);
	/* reserve some memory */
	binsVec.reserve(MISS_BAR);

	/* first bin is SD=0 */
	binsVec.push_back(buffer[0]);
	/* calculate the sum of all samples */
	samples = buffer[0];

	/* i is the index of buffer and starts from 1 */
	for (long i = 1; i < bufSize; ++i)
		/* fill up with 0 for the blank bins */
		if (!buffer[i])
			binsVec.resize(binsVec.size() + std::exp2(i - 1));
		else {
			for (int j = std::exp2(i - 1); j < std::exp2(i) && j <= MISS_BAR; ++j)
				binsVec.push_back((Accur)buffer[i] / std::exp2(i - 1));
			samples += buffer[i];
		}

	return binsVec.size() == MISS_BAR;
}

template <class B, class Accur>
bool Histogram<B, Accur>::mapToVector(std::vector<B> & buffer)
{
	assert(buffer.size());
	/* reserve some memory */
	binsVec.reserve(MISS_BAR + 1);

	/* first bin is SD=0 */
	binsVec.push_back(buffer[0]);
	/* calculate the sum of all samples */
	samples = buffer[0];

	/* i is the index of buffer and starts from 1 */
	for (long i = 1; i < buffer.size(); ++i)
		/* fill up with 0 for the blank bins */
		if (!buffer[i])
			binsVec.resize(binsVec.size() + std::exp2(i - 1));
		else {
			for (int j = std::exp2(i - 1); j < std::exp2(i) && j <= MISS_BAR; ++j)
				binsVec.push_back((Accur)buffer[i] / std::exp2(i - 1));
			samples += buffer[i];
		}

		return binsVec.size() == MISS_BAR + 1;
}

template <class B, class Accur>
void Histogram<B, Accur>::reuseDistToStackDist()
{
	std::vector<Accur> frac;
	B temp = 0;
	/* calculate the fraction of reuse distance
	is greater than i */
	for (int i = 0; i < binsVec.size(); ++i) {
		temp += binsVec[i];
		frac.push_back((Accur)(samples - temp) / samples);
	}

	std::vector<B> transTemp(binsVec);

	for (int i = 1; i < binsVec.size(); ++i) {
		Accur sumFrac = 0.0;
		for (int j = 1; j <= i; ++j)
			sumFrac += frac[j];
		/* now we use the probabilities in RD to transform into SD, F(r)=s */
		transTemp[(int) std::round(sumFrac)] += binsVec[i];
	}
	transTemp[0] = binsVec[0];

	binsVec.clear();
	binsVec = transTemp;
}

/* calculate log2(s) + 1 */
template<class T>
inline T log2p1(T s)
{
	T result = 0;
	while (s) {
		s >>= 1;
		++result;
	}

	return result;
}

template <class B, class Accur>
Accur Histogram<B, Accur>::fullyToSetAssoc(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	Accur p = (Accur)1 / setNum;
	binsTra.clear();

	/* secureDist = log2(assoc),
	the least not evcited SD, called secure distance */
	int secureDist = log2p1(assoc) - 1;

	/* if a mem ref's SD <= secureDist, it is a cache hit, 
	   we don't need to transform these bar */
	for (int i = 0; i <= secureDist; ++i)
		binsTra.push_back(binsVec[i]);
	binsTra.resize(binsVec.size());

	try {
		/* transition of each bar */
		for (int i = secureDist + 1; i < binsVec.size(); ++i) {
			if (!binsVec[i]) continue;

			/* each bar j before the current, need to be added,
			it's a binomial distribution */
			boost::math::binomial binom(i, p);

			/* start from j=0, ref with SD >= assoc + secureDist is definitely miss in PLRU cache,
			   so we don't need to calculate these bars */
			for (int j = 0; j <= i && j < assoc + secureDist; ++j)
				binsTra[j] += binsVec[i] * boost::math::pdf(binom, (Accur)j);

			//for (int j = 0; j <= k; ++j)
			//	binsTra[log2p1(j)] += (Accur)binsVec[i] / std::exp2(i - 1) * boost::math::pdf(binom, (Accur)j);// shen
		}

		hits = 0;
		for (int i = 0; i < assoc; ++i) // shen
			hits += (B)std::round(binsTra[i]);
		return 1 - (Accur)hits / samples;
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}
}

template <class B, class Accur>
void Histogram<B, Accur>::fullyToSetAssoc_Poisson(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	Accur p = (Accur)1 / setNum;

	/* if a mem ref's stack distance <= assoc - 1, it is a cache hit */
	for (int i = 0; i < assoc; ++i)
		binsTra.push_back(binsVec[i]);

	binsTra.resize(binsVec.size());

	/* transition of each bar */
	for (int i = assoc; i < binsVec.size(); ++i) {
		if (!binsVec[i]) continue;

		/* each bar j before the current, need to be added,
		it's a binomial distribution */
		boost::math::poisson_distribution<Accur> pois(i * p);

		/* start from j=i, because i is the current bar */
		for (int j = i; j >= 0; --j)
			binsTra[j] += binsVec[i] * boost::math::pdf(pois, j);
	}
}

template <class B, class Accur>
Accur Histogram<B, Accur>::calLruMissRate(const int & cap, const int & blk, const int & assoc)
{
	hits = 0;

	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	Accur p = (Accur)1 / setNum;

	/* we use cdf to calculate the probability that 
	   other SDs will be less than assoc after transforming */
	for (int i = assoc; i < binsVec.size(); ++i) {
		boost::math::binomial binom(i, p);
		hits += (B)std::round(binsVec[i] * boost::math::cdf(binom, (Accur)assoc - 1));
	}

	/* add up the SD < assoc fraction */
	for (int i = 0; i < assoc; ++i)
		hits += (B)std::round(binsVec[i]);

	missRate = 1 - (Accur)hits / samples;
	return missRate;
}

/* PLRU hit function is a piecewise linear function. */
template <class B, class Accur>
Accur Histogram<B, Accur>::calPlruMissRate(const int & cap, const int & blk, const int & assoc)
{
	fullyToSetAssoc(cap, blk, assoc);

	/* secureDist = log2(assoc),
	the least not evcited distance, called secure distance */
	int secureDist = log2p1(assoc) - 1;

	/* we assume that refs with SD <= secureDist are hits, SD > secureDist is a miss,
	   so the miss rate is SD(secureDist + 1) / sum(SD, secureDist + 1, 0). */
	Accur tempMisses = 0.0;
	Accur tempHits = 0.0;
	for (int i = 0; i <= secureDist; ++i)
		tempHits += binsTra[i];
	tempMisses = binsTra[secureDist + 1];
	Accur lastAccesMiss = tempMisses / (tempMisses + tempHits);

	//Accur tempMisses = 0.0;
	//Accur tempHits = 0.0;
	//for (int i = 0; i <= secureDist; ++i)
	//	tempHits += binsTra[i];
	//for (int i = secureDist + 1; i < assoc; ++i)
	//	tempMisses += binsTra[i];
	//Accur lastAccesMiss = tempMisses / (tempMisses + tempHits);

	int tempAssoc = assoc - 1;
	/* p1, a point with SD=secureDist+1 */
	Accur p1 = 1.0;
	/* 2^0 * 2^1 * ... * 2^(secureDist-1) / A(assoc-1, secureDist) */
	for (int i = 0; i < secureDist; ++i) {
		p1 *= (1 << i);
		p1 /= tempAssoc--;
	}
	/* the last ref is a miss */
	p1 *= lastAccesMiss;
	/* hit probability */
	p1 = 1 - p1;

	/* p2, a point with SD=assoc */
	Accur p2 = 1.0;
	/* accesses before current cousins */
	Accur front = 0.0;
	try {
		/* calculate the permutations */
		for (int j = 0; j < secureDist; ++j) {
			/* number of accesses in the serialization in same subtree */
			Accur cousins = (1 << j);
			/* tgammma_ratio(a, b) for (a-1)! / (b-1)! */
			p2 *= cousins * boost::math::tgamma_ratio(front + cousins, front + 1);
			front += cousins;
		}
		/* divide all possible permutation */
		p2 /= boost::math::tgamma((Accur)assoc);
		/* hit probability */
		p2 = 1 - p2;
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}

	hits = 0;
	/* the gradient of first piece */
	Accur k1 = (p2 - p1) / (assoc - secureDist - 1);
	for (int i = secureDist + 1; i <= assoc; ++i) {
		Accur y = k1 * (i - secureDist - 1) + p1;
		hits += y * binsTra[i];
		//std::cout << y << std::endl;
	}
	/* p3, with SD=assoc + assoc/2 - 2, hit probability */
	Accur p3 = 1 / (Accur)(assoc - 1);
	/* the gradient of second piece */
	Accur k2 = (p3 - p2) / (assoc / 2 - 2);
	for (int i = assoc + 1; i <= assoc * 3 / 2 - 2; ++i) {
		Accur y = k2 * (i - assoc) + p2;
		hits += y * binsTra[i];
		//std::cout << y << std::endl;
	}
	
	for (int i = 0; i <= secureDist; ++i)
		hits += binsTra[i];

	missRate = 1 - (Accur)hits / samples;
	return missRate;
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

template <class B, class Accur>
Accur Histogram<B, Accur>::calMissRate(const int & cap, const int & blk, const int & assoc, const bool plru)
{
	if (!plru)
		return calLruMissRate(cap, blk, assoc);
	else
		return calPlruMissRate(cap, blk, assoc);

	////fullyToSetAssoc(cap, blk, assoc);

	//int tempAssoc = assoc - 1;
	///* secureDist = log2(assoc),
	//   the least not evcited distance, called secure distance */
	//int secureDist = 0;
	//while (tempAssoc) {
	//	tempAssoc >>= 1;
	//	++secureDist;
	//}

	///* Now we calculate the approximate miss rate in an assumed LRU cache
	//   that the number of blocks equals to association,
	//   and references with stack distance <= log2(assoc) must be hit. */
	////Accur tempMisses = 0.0;
	////Accur tempHits = 0.0;
	////for (int i = 0; i <= secureDist; ++i)
	////	tempHits += binsTra[i];
	////for (int i = secureDist + 1; i < assoc; ++i)
	////	tempMisses += binsTra[i];
	////Accur lastAccesMiss = tempMisses / (tempMisses + tempHits);

	//Accur lastAccesMiss = 0.15;
	//Accur missesD = 0.0;
	//try {
	//	for (int i = secureDist + 1; i < assoc; ++i) {
	//		/* probability for eviction */
	//		Accur p = 1.0;
	//		/* probability for accessing in-order */
	//		Accur pInOrder = 1.0;
	//		/* accesses before current cousins */
	//		Accur front = 0.0;
	//		for (int j = 0; j < secureDist; ++j) {
	//			Accur c = combinationRatio(assoc - 1, assoc - 1 - (1 << j), i - 1);
	//			p *= 1 - c;
	//			/* number of accesses in the serialization in same subtree */
	//			Accur cousins = std::max(std::floor((1 << j) * (i - 1) / (assoc - 1)), 1.0);
	//			/* tgammma_ratio(a, b) for (a-1)! / (b-1)! */
	//			pInOrder *= cousins * boost::math::tgamma_ratio(front + cousins, front + 1);

	//			front += cousins;
	//		}
	//		/* divide all possible permutation */
	//		pInOrder /= boost::math::tgamma(i);
	//		if (pInOrder >= 1)
	//			throw std::exception("Error: calMissRate: wrong calculation");
	//		missesD += p * pInOrder * lastAccesMiss; //* binsTra[i];
	//		std::cout << 1 - p * pInOrder * lastAccesMiss << std::endl;
	//	}
	//}
	//catch (std::exception e) {
	//	std::cout << e.what() << std::endl;
	//	return -INFINITY;
	//}

	//misses = std::floorl(missesD);
	///* references with SD > assoc always miss. */
	//for (int i = assoc; i < binsVec.size(); ++i)
	//	misses += binsTra[i];

	//missRate = (Accur)misses / samples;

	//return missRate;
}

template <class B, class Accur>
void Histogram<B, Accur>::print(std::ofstream & file)
{
	auto it = binsMap.begin();
	auto last = --binsMap.end();
	for (int i = 0; i < last->first; ++i)
		/* print zero value bins */
		if (it->first != i)
			file << "0 ";
		else {
			file << it->second << " ";
			++it;
		}
		file << "\n";

}

template <class Accur>
AvlNode<Accur>::AvlNode(Accur & a) : holes(1), rHoles(0), height(0), left(nullptr), right(nullptr)
{
	interval = std::make_pair(a, a);
}

template <class Accur>
AvlNode<Accur>::AvlNode(AvlNode<Accur> & n) : holes(n.holes), rHoles(n.rHoles), height(n.height), left(n.left), right(n.right)
{
	interval = n.interval;
}

template <class Accur>
AvlNode<Accur>::~AvlNode()
{
	delete left;
	delete right;
	left = nullptr;
	right = nullptr;
}

template <class Accur>
void AvlNode<Accur>::updateHeight()
{
	height = 1 + std::max(getHeight(left), getHeight(right));
}

template <class Accur>
void AvlNode<Accur>::updateHoles()
{
	rHoles = right ? right->holes : 0;

	holes =
		(left ? left->holes : 0) + rHoles +
		int(interval.second - interval.first) + 1;
}

AvlTreeStack::AvlTreeStack(long & v) : index(0), curHoles(0)
{
	root = new AvlNode<long>(v);
}

AvlTreeStack::AvlTreeStack() : root(nullptr), index(0), curHoles(0) {};

void AvlTreeStack::destroy(AvlNode<long> * & tree)
{
	if (!tree)
		return;

	destroy(tree->left);
	destroy(tree->right);
	delete tree;
	tree = nullptr;
}

void AvlTreeStack::clear()
{
	addrMap.clear();
	destroy(root);
	index = -1;
	curHoles = -1;
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

void AvlTreeStack::calStackDist(uint64_t addr, Histogram<> & hist)
{
	long & value = addrMap[addr];

	++index;

	/* value is 0 under cold miss */
	if (!value) {
		value = index;
		hist.sample(log2p1(MISS_BAR));
		return;
	}

	/* update b of last reference */
	if (value < index) {
		/* insert a hole */
		insert(value);
		int stackDist = index - value - curHoles - 1;
		/* if the stack distance is large than MISS_BAR, the reference is definitely missed. */
		stackDist = stackDist >= MISS_BAR ? MISS_BAR : stackDist;
		hist.sample(log2p1(stackDist));
	}

	value = index;
}

void ReuseDist::calReuseDist(uint64_t addr, Histogram<> & hist)
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

SampleStack::SampleStack(int sSize, int hSize, double sRate) : sampleCounter(0), statusCounter(0), hibernInter(hSize),
sampleInter(sSize), state(State::hibernating), sampleRate(sRate) {};

SampleStack::SampleStack(double sRate) : sampleCounter(0), statusCounter(0), hibernInter(0),
sampleInter(0), state(State::sampling), sampleRate(sRate) {};

int SampleStack::genRandom()
{
	// construct a trivial random generator engine from a time-based seed:
	std::chrono::system_clock::rep seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 engine(seed);
	static std::geometric_distribution<int> geom(sampleRate);
	int rand = 0;
	/* the random can not be 0 */
	while (!rand)
		rand = geom(engine);

	return rand;
}

void SampleStack::calStackDist(uint64_t addr, Histogram<> & hist)
{
	/* start a new sampling interval */
	if (state == State::hibernating && !statusCounter) {
		statusCounter = sampleInter;
		sampleCounter = genRandom();
		state = State::sampling;
	}
	/* start hibernation interval */
	else if (state == State::sampling && !statusCounter) {
		statusCounter = hibernInter;
		state = State::hibernating;
	}

	/* if we find a same address x in addrTable,
	record its stack distance and the sampling of x is finished */
	auto pos = addrTable.find(addr);
	if (pos != addrTable.end()) {
		hist.sample(log2p1(pos->second.size() - 1));
		addrTable.erase(addr);
	}

	auto it = addrTable.begin();

	/* make sure the max size of addrTable */
	if (addrTable.size() > sampleRate * sampleInter * 2) {
		hist.sample(log2p1(MISS_BAR));
		addrTable.erase(it->first);
	}

	/* record unique mem references between the sampled address x */
	for (it = addrTable.begin(); it != addrTable.end(); ) {
		it->second.insert(addr);
		/* if the set of sampled address x is too large,
		erase it from  the table and record as MISS_BAR */
		if (it->second.size() > MISS_BAR) {
			auto eraseIt = it;
			++it;
			hist.sample(log2p1(MISS_BAR));
			addrTable.erase(eraseIt->first);
		}
		else
			++it;
	}

	/* if it is time to do sampling */
	if (state == State::sampling && !sampleCounter) {
		/* it is a new sampled address */
		assert(!addrTable[addr].size());
		addrTable[addr].insert(0);
		/* reset the sampleCounter and randNum to prepare next sample */
		sampleCounter = genRandom();
	}

	--statusCounter;
	if (state == State::sampling)
		--sampleCounter;
}

Reader::Reader(std::ifstream & fin) {
	std::ofstream fout("missrate-gcc-166-100M.txt");
	Histogram<int64_t> histogram;
	int binSize = log2p1(MISS_BAR) + 1;
	int64_t * temp = new int64_t[binSize];
	std::vector<int64_t> buffer(binSize);

	while (!fin.eof()) {
		fin.read((char *)temp, binSize * sizeof(int64_t));
		if (temp[binSize - 1] < 0) break;
		//for (int i = 0; i < binSize; ++i)
		//	std::cout << temp[i] << " ";
		//std::cout << "\n";
		for (int i = 0; i < binSize; ++i) {
			buffer[i] += temp[i];
		}


		///* save temp bins into histogram */
		//bool succ = histogram.mapToVector(temp, binSize);
		//
		//int cap = 32 * 1024;
		//int blk = 64;
		//int assoc = 16;

		//double lruMissRate = histogram.calMissRate(cap, blk, assoc, false);
		////double lruTraMissRate = histogram.fullyToSetAssoc(cap, blk, assoc);
		//double plruMissRate = histogram.calMissRate(cap, blk, assoc, true);

		//fout << std::setprecision(6) << lruMissRate << " " << plruMissRate << std::endl;
	}

	delete[] temp;

	/* save temp bins into histogram */
	bool succ = histogram.mapToVector(buffer);

	int cap = 32 * 1024;
	int blk = 64;
	int assoc = 16;

	double lruMissRate = histogram.calMissRate(cap, blk, assoc, false);
	//double lruTraMissRate = histogram.fullyToSetAssoc(cap, blk, assoc);
	double plruMissRate = histogram.calMissRate(cap, blk, assoc, true);

	fout << std::setprecision(6) << lruMissRate << " " << plruMissRate << std::endl;
	fout.close();
}

int main()
{
	//vector a = generateA();
	//double z = 7.0;
	//double c = lanzocsConstWithE(z);
	//double d = boost::math::lanczos::lanczos13m53::lanczos_sum_expG_scaled(z);
	//boost::math::binomial binom(100, 0.4);
	//double g = boost::math::pdf(binom, 20);
	//Histogram<> hist;
	//hist.sample(log2p1(0));
	//hist.sample(log2p1(8));
	//hist.sample(log2p1(13));
	//bool c = hist.mapToVector();
	//missProb(8);
	//std::ofstream out("E:\\ShareShen\\dump.txt");
	//Histogram<> hist;
	//hist.sample(9);
	//hist.sample(10);
	//hist.print(out);
	//out.close();
	//int a[10] = { 1,2,3 };
	//void * q = a;
	//uint64_t c = reinterpret_cast<uint64_t> (q);
	std::ifstream fin;
	fin.open("E:\\ShareShen\\pin-3.2-81205-gcc-linux\\source\\tools\\memTraceSimple\\gcc-trace\\gcc-166-0.0002-100M.txt", std::ios::binary | std::ios::in);
	if (fin.fail()) {
		std::cout << "SDD file openning failed! " << std::endl;
		return 0;
	}
	Reader reader(fin);

 	return 0;
}