// StackDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StackDistance.h"


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
	if (tree == nullptr) {
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

Reader::Reader(const std::string _path, const std::string _pathout, const std::string _pathset) {
	/* read set distributions of all intervals */
	intSetDistr.recordSetDistr(_pathset);

	file.open(_path, std::ios::in);
	if (file.fail()) {
		return;
	}

	std::string temp;
	bool root = true;
	AvlTreeStack * avlTreeStack = nullptr;
	int it = 0;

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		uint64_t paddr;
		uint64_t vaddr;
		lineStream >> temp;
		/* a new interval starts */
		if (temp == "interval") {
			/* if there is an old interval, 
			   we transit the histogram and push it into intHists. */
			if (avlTreeStack != nullptr) {
				avlTreeStack->transStack(intSetDistr.intSetDistr[it]);
				intHists.push_back(avlTreeStack->hist);
				++it;
				delete avlTreeStack;
			}

			avlTreeStack = new AvlTreeStack;
			std::cout << "interval " << it << std::endl;
			continue;
		}

		assert(avlTreeStack != nullptr);

		lineStream >> std::hex >> paddr;
		lineStream >> temp;
		lineStream >> std::hex >> vaddr;

		uint64_t addr = paddr;
		/* check a valid address */
		if (!addr)
			continue;
		
		/* call the methods to calculate stack distances */
		//listStack.calStackDist(addr);
		uint64_t mask = 63;
		addr = addr & (~mask);
		avlTreeStack->calStackDist(addr);
	}
	/* the last interval */
	avlTreeStack->transStack(intSetDistr.intSetDistr[it]);
	intHists.push_back(avlTreeStack->hist);
	delete avlTreeStack;
	
	outFile.open(_pathout, std::ios::out);
	if (outFile.fail())
		return;
	
	sumHists();
	//listStack.print("E:\\ShareShen\\gem5-origin\\m5out-se-x86\\perlbench.txt");
	file.close();
	outFile.close();
}

void Reader::sumHists()
{
	std::vector<double> sumBins;
	for (int b = 0; b < intHists[0].transBins.size(); ++b) {
		double temp = 0.0;

		for (int i = 0; i < intHists.size(); ++i)
			temp += intHists[i].transBins[b];

		sumBins.push_back(temp);
		if (!temp) continue;
		outFile << "dist::" << b << "\t" << temp << std::endl;
	}
}

int main()
{
	//IntervalSetDistribution intSetDistr;
	//intSetDistr.recordSetDistr("E:\\ShareShen\\gem5-stable\\m5out-se-x86\\perlbench-32k4assoc-setDistr.txt");

	Reader reader("E:\\ShareShen\\gem5-stable\\m5out-se-x86\\requtTraceFile-perlbench.txt", \
		"E:\\ShareShen\\gem5-stable\\m5out-se-x86\\perlbench-avl-4assoc.txt", \
		"E:\\ShareShen\\gem5-stable\\m5out-se-x86\\perlbench-32k4assoc-setDistr.txt");

	return 0;
}