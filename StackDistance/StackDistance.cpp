// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"
#include "NtlTest.h"

//template <class T, class Lanczos, class Policy>
//T ibeta_power_terms(T a,
//	T b,
//	T x,
//	T y,)
//{
//	T result;
//
//	T c = a + b;
//
//	// combine power terms with Lanczos approximation:
//	T agh = static_cast<T>(a + G - 0.5);
//	T bgh = static_cast<T>(b + G - 0.5);
//	T cgh = static_cast<T>(c + G - 0.5);
//	result = lanzocsConstWithE(c) / (lanzocsConstWithE(a) * lanzocsConstWithE(b));
//	//result *= prefix;
//	// combine with the leftover terms from the Lanczos approximation:
//	result *= sqrt(bgh / boost::math::constants::e<T>());
//	result *= sqrt(agh / cgh);
//
//	// l1 and l2 are the base of the exponents minus one:
//	T l1 = (x * b - y * agh) / agh;
//	T l2 = (y * a - x * bgh) / bgh;
//	if (((std::min)(fabs(l1), fabs(l2)) < 0.2))
//	{
//		// when the base of the exponent is very near 1 we get really
//		// gross errors unless extra care is taken:
//		if ((l1 * l2 > 0) || ((std::min)(a, b) < 1))
//		{
//			//
//			// This first branch handles the simple cases where either: 
//			//
//			// * The two power terms both go in the same direction 
//			// (towards zero or towards infinity).  In this case if either 
//			// term overflows or underflows, then the product of the two must 
//			// do so also.  
//			// *Alternatively if one exponent is less than one, then we 
//			// can't productively use it to eliminate overflow or underflow 
//			// from the other term.  Problems with spurious overflow/underflow 
//			// can't be ruled out in this case, but it is *very* unlikely 
//			// since one of the power terms will evaluate to a number close to 1.
//			//
//			if (fabs(l1) < 0.1)
//			{
//				result *= exp(a * boost::math::log1p(l1, pol));
//				BOOST_MATH_INSTRUMENT_VARIABLE(result);
//			}
//			else
//			{
//				result *= pow((x * cgh) / agh, a);
//				BOOST_MATH_INSTRUMENT_VARIABLE(result);
//			}
//			if (fabs(l2) < 0.1)
//			{
//				result *= exp(b * boost::math::log1p(l2, pol));
//				BOOST_MATH_INSTRUMENT_VARIABLE(result);
//			}
//			else
//			{
//				result *= pow((y * cgh) / bgh, b);
//				BOOST_MATH_INSTRUMENT_VARIABLE(result);
//			}
//		}
//		else if ((std::max)(fabs(l1), fabs(l2)) < 0.5)
//		{
//			//
//			// Both exponents are near one and both the exponents are 
//			// greater than one and further these two 
//			// power terms tend in opposite directions (one towards zero, 
//			// the other towards infinity), so we have to combine the terms 
//			// to avoid any risk of overflow or underflow.
//			//
//			// We do this by moving one power term inside the other, we have:
//			//
//			//    (1 + l1)^a * (1 + l2)^b
//			//  = ((1 + l1)*(1 + l2)^(b/a))^a
//			//  = (1 + l1 + l3 + l1*l3)^a   ;  l3 = (1 + l2)^(b/a) - 1
//			//                                    = exp((b/a) * log(1 + l2)) - 1
//			//
//			// The tricky bit is deciding which term to move inside :-)
//			// By preference we move the larger term inside, so that the
//			// size of the largest exponent is reduced.  However, that can
//			// only be done as long as l3 (see above) is also small.
//			//
//			bool small_a = a < b;
//			T ratio = b / a;
//			if ((small_a && (ratio * l2 < 0.1)) || (!small_a && (l1 / ratio > 0.1)))
//			{
//				T l3 = boost::math::expm1(ratio * boost::math::log1p(l2, pol), pol);
//				l3 = l1 + l3 + l3 * l1;
//				l3 = a * boost::math::log1p(l3, pol);
//				result *= exp(l3);
//				BOOST_MATH_INSTRUMENT_VARIABLE(result);
//			}
//			else
//			{
//				T l3 = boost::math::expm1(boost::math::log1p(l1, pol) / ratio, pol);
//				l3 = l2 + l3 + l3 * l2;
//				l3 = b * boost::math::log1p(l3, pol);
//				result *= exp(l3);
//				BOOST_MATH_INSTRUMENT_VARIABLE(result);
//			}
//		}
//		else if (fabs(l1) < fabs(l2))
//		{
//			// First base near 1 only:
//			T l = a * boost::math::log1p(l1, pol)
//				+ b * log((y * cgh) / bgh);
//			if ((l <= tools::log_min_value<T>()) || (l >= tools::log_max_value<T>()))
//			{
//				l += log(result);
//				if (l >= tools::log_max_value<T>())
//					return policies::raise_overflow_error<T>(function, 0, pol);
//				result = exp(l);
//			}
//			else
//				result *= exp(l);
//			BOOST_MATH_INSTRUMENT_VARIABLE(result);
//		}
//		else
//		{
//			// Second base near 1 only:
//			T l = b * boost::math::log1p(l2, pol)
//				+ a * log((x * cgh) / agh);
//			if ((l <= tools::log_min_value<T>()) || (l >= tools::log_max_value<T>()))
//			{
//				l += log(result);
//				if (l >= tools::log_max_value<T>())
//					return policies::raise_overflow_error<T>(function, 0, pol);
//				result = exp(l);
//			}
//			else
//				result *= exp(l);
//			BOOST_MATH_INSTRUMENT_VARIABLE(result);
//		}
//	}
//	else
//	{
//		// general case:
//		T b1 = (x * cgh) / agh;
//		T b2 = (y * cgh) / bgh;
//		l1 = a * log(b1);
//		l2 = b * log(b2);
//		BOOST_MATH_INSTRUMENT_VARIABLE(b1);
//		BOOST_MATH_INSTRUMENT_VARIABLE(b2);
//		BOOST_MATH_INSTRUMENT_VARIABLE(l1);
//		BOOST_MATH_INSTRUMENT_VARIABLE(l2);
//		if ((l1 >= tools::log_max_value<T>())
//			|| (l1 <= tools::log_min_value<T>())
//			|| (l2 >= tools::log_max_value<T>())
//			|| (l2 <= tools::log_min_value<T>())
//			)
//		{
//			// Oops, under/overflow, sidestep if we can:
//			if (a < b)
//			{
//				T p1 = pow(b2, b / a);
//				T l3 = a * (log(b1) + log(p1));
//				if ((l3 < tools::log_max_value<T>())
//					&& (l3 > tools::log_min_value<T>()))
//				{
//					result *= pow(p1 * b1, a);
//				}
//				else
//				{
//					l2 += l1 + log(result);
//					if (l2 >= tools::log_max_value<T>())
//						return policies::raise_overflow_error<T>(function, 0, pol);
//					result = exp(l2);
//				}
//			}
//			else
//			{
//				T p1 = pow(b1, a / b);
//				T l3 = (log(p1) + log(b2)) * b;
//				if ((l3 < tools::log_max_value<T>())
//					&& (l3 > tools::log_min_value<T>()))
//				{
//					result *= pow(p1 * b2, b);
//				}
//				else
//				{
//					l2 += l1 + log(result);
//					if (l2 >= tools::log_max_value<T>())
//						return policies::raise_overflow_error<T>(function, 0, pol);
//					result = exp(l2);
//				}
//			}
//			BOOST_MATH_INSTRUMENT_VARIABLE(result);
//		}
//		else
//		{
//			// finally the normal case:
//			result *= pow(b1, a) * pow(b2, b);
//			BOOST_MATH_INSTRUMENT_VARIABLE(result);
//		}
//	}
//
//	BOOST_MATH_INSTRUMENT_VARIABLE(result);
//
//	return result;
//}

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

template <class B, class Accur>
bool Histogram<B, Accur>::mapToVector()
{
	std::map<long, B>::iterator last = binsMap.end();
	/* point to the last element */
	long maxSize = (--last)->first;
	/* reserve some memory */
	binsVec.reserve(maxSize);

	std::map<long, B>::iterator it = binsMap.begin();

	for (long i = 0; i <= maxSize; ++i)
		if (i == it->first) {
			binsVec.push_back(it->second);
			++it;
		}
		else
			binsVec.push_back(0);

		return binsVec.size() == last->first + 1;
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

template <class B, class Accur>
Accur Histogram<B, Accur>::fullyToSetAssoc(const int & cap, const int & blk, const int & assoc)
{
	int setNum = cap / (blk * assoc);
	/* p is a probability of mapping to a specific set */
	Accur p = (Accur)1 / setNum;
	binsTra.clear();

	int tempAssoc = assoc - 1;
	/* secureDist = log2(assoc),
	the least not evcited SD, called secure distance */
	int secureDist = 0;
	while (tempAssoc) {
		tempAssoc >>= 1;
		++secureDist;
	}

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
		}

		hits = 0;
		for (int i = 0; i < assoc; ++i)
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
		hits += binsVec[i];

	missRate = 1 - (Accur)hits / samples;
	return missRate;
}

/* PLRU hit function is a piecewise linear function. */
template <class B, class Accur>
Accur Histogram<B, Accur>::calPlruMissRate(const int & cap, const int & blk, const int & assoc)
{
	fullyToSetAssoc(cap, blk, assoc);

	int tempAssoc = assoc - 1;
	/* secureDist = log2(assoc),
	the least not evcited distance, called secure distance */
	int secureDist = 0;
	while (tempAssoc) {
		tempAssoc >>= 1;
		++secureDist;
	}

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

	tempAssoc = assoc - 1;
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
		std::cout << y << std::endl;
	}
	/* p3, with SD=assoc + assoc/2 - 2, hit probability */
	Accur p3 = 1 / (Accur)(assoc - 1);
	/* the gradient of second piece */
	Accur k2 = (p3 - p2) / (assoc / 2 - 2);
	for (int i = assoc + 1; i <= assoc * 3 / 2 - 2; ++i) {
		Accur y = k2 * (i - assoc) + p2;
		hits += y * binsTra[i];
		std::cout << y << std::endl;
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

void AvlTreeStack::calReuseDist(uint64_t addr, Histogram<> & hist)
{
	long & value = addrMap[addr];

	++index;

	/* value is 0 under cold miss */
	if (!value) {
		value = index;
		hist.sample(MISS_BAR);
		return;
	}

	/* update b of last reference */
	if (value < index) {
		/* insert a hole */
		insert(value);
		int stackDist = index - value - curHoles - 1;
		/* if the stack distance is large than MISS_BAR, the reference is definitely missed. */
		stackDist = stackDist >= MISS_BAR ? MISS_BAR : stackDist;
		hist.sample(stackDist);
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

int SampleStack::genRandom()
{
	// construct a trivial random generator engine from a time-based seed:
	std::chrono::system_clock::rep seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 engine(seed);
	static std::geometric_distribution<int> geom((double)expectSamples / sampleInter);
	int rand = 0;
	/* the random can not be 0 */
	while (!rand)
		rand = geom(engine);

	return rand;
}

void SampleStack::calStackDist(uint64_t addr, Histogram<> & hist)
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
	record its stack distance and the sampling of x is finished */
	std::unordered_map<uint64_t, AddrSet>::iterator pos;
	pos = addrTable.find(addr);
	if (pos != addrTable.end()) {
		hist.sample(pos->second.size() - 1);
		addrTable.erase(addr);
	}

	std::unordered_map<uint64_t, AddrSet>::iterator it = addrTable.begin();

	/* make sure the max size of addrTable */
	if (addrTable.size() > expectSamples) {
		hist.sample(MISS_BAR);
		addrTable.erase(it->first);
	}

	/* record unique mem references between the sampled address x */
	for (it = addrTable.begin(); it != addrTable.end(); ) {
		it->second.insert(addr);
		/* if the set of sampled address x is too large,
		erase it from  the table and record as MISS_BAR */
		if (it->second.size() > MISS_BAR) {
			std::unordered_map<uint64_t, AddrSet>::iterator eraseIt = it;
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
		assert(!addrTable[addr].size());
		addrTable[addr].insert(0);
		/* reset the sampleCounter and randNum to prepare next sample */
		sampleCounter = 0;
		randNum = genRandom();
		//randNum = 1;
	}
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
	/* the index of interval */
	int i = 2;

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		uint64_t paddr;
		uint64_t vaddr;
		lineStream >> temp;
		
		if (temp == "interval") {
			++interval;
			if (interval == i) {
				start = true;
				continue;
			}
			else if (interval > i)
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
		avlTreeStack.calReuseDist(addr, histogram);
#endif

#ifdef REUSE
		reuseDist.calReuseDist(addr, histogram);
#endif

#ifdef SAMPLE
		sampleStack.calStackDist(addr, histogram);
#endif
	}
	//listStack.print("E:\\ShareShen\\gem5-origin\\m5out-se-x86\\perlbench.txt");
	std::cout << "transiting stack distances..." << std::endl;

	int cap = 32 * 1024;
	int blk = 64;
	int assoc = 16;
	/* transforming and calculating miss rate */
	bool succ = histogram.mapToVector();

#ifdef REUSE
	histogram.reuseDistToStackDist();
	histogram.calMissRate(512);
#endif // REUSE

	double lruMissRate = histogram.calMissRate(cap, blk, assoc, false);
	//double lruTraMissRate = histogram.fullyToSetAssoc(cap, blk, assoc);
	double plruMissRate2 = histogram.calMissRate(cap, blk, assoc, true);
	double plruMissRate = histogram.calMissRate(cap, blk, assoc, true);

	file.close();
	fileOut.close();
}

void missProb(int assoc)
{
	int tempAssoc = assoc - 1;
	/* secureDist = log2(assoc),
	the least not evcited distance, called secure distance */
	int secureDist = 0;
	while (tempAssoc) {
		tempAssoc >>= 1;
		++secureDist;
	}

	double p = 1.0;
	tempAssoc = assoc - 1;

	for (int i = 0; i < secureDist; ++i) {
		p *= (1 << i);
		p /= tempAssoc--;
	}

	double missRate = 0.2;
	p *= missRate;
	std::cout << 1 - p << std::endl;

	for (int i = 0; i < assoc * 3 / 2; ++i) {
		p = (1 - p);
		double t = std::pow(0.5, secureDist);
		p *= t;
		missRate -= 0.02;
		p *= missRate;
		std::cout << 1 - p << std::endl;
	}
}

int main()
{
	//vector a = generateA();
	//double z = 7.0;
	//double c = lanzocsConstWithE(z);
	//double d = boost::math::lanczos::lanczos13m53::lanczos_sum_expG_scaled(z);
	//boost::math::binomial binom(100, 0.4);
	//double g = boost::math::pdf(binom, 20);
	/*Histogram<> hist;
	hist.calMissRate(0, 0, 8, true);*/
	
	Reader reader("E:\\ShareShen\\gem5-stable\\m5out-se-x86\\cactusADM\\cactusADM-trace-part.txt", \
		"E:\\ShareShen\\gem5-stable\\m5out-se-x86\\cactusADM\\cactusADM-avl-2assoc.txt");

 	return 0;
}