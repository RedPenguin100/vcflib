#include "BedReader.h"

#include <iostream>

#include "split.h"


basic_string<char> strip(string const& str, char const* separators = " \t")
{
	string::size_type const first = str.find_first_not_of(separators);
	return (first == string::npos) ? string()
		: str.substr(first, str.find_last_not_of(separators) - first + 1);
}

void parseRegion(
	const string& region,
	string& startSeq,
	int& startPos,
	int& stopPos) {

	const size_t foundFirstColon = region.find(':');

	// we only have a single string, use the whole sequence as the target
	if (foundFirstColon == string::npos) {
		startSeq = region;
		startPos = 0;
		stopPos = -1;
	}
	else {
		startSeq = region.substr(0, foundFirstColon);
		string sep = "..";
		size_t foundRangeSep = region.find(sep, foundFirstColon);
		if (foundRangeSep == string::npos) {
			sep = "-";
			foundRangeSep = region.find('-', foundFirstColon);
		}
		if (foundRangeSep == string::npos) {
			startPos = atoi(region.substr(foundFirstColon + 1).c_str());
			// differ from bamtools in this regard, in that we process only
			// the specified position if a range isn't given
			stopPos = startPos + 1;
		}
		else {
			startPos = atoi(region.substr(foundFirstColon + 1, foundRangeSep - foundFirstColon).c_str());
			// if we have range sep specified, but no second number, read to the end of sequence
			if (foundRangeSep + sep.size() != region.size()) {
				stopPos = atoi(region.substr(foundRangeSep + sep.size()).c_str()); // end-exclusive, bed-format
			}
			else {
				//stopPos = reference.sequenceLength(startSeq);
				stopPos = -1;
			}
		}
	}
}

vector<BedTarget> BedReader::entries()
{
	vector<BedTarget> entries;

	if (!isOpen()) {
		cerr << "bed targets file is not open" << endl;
		exit(1);
	}

	string line;
	while (std::getline(file, line)) {
		vector<string> fields = split(line, " \t");
		BedTarget entry(strip(fields[0]),
			atoi(strip(fields[1]).c_str()),
			atoi(strip(fields[2]).c_str()),
			(fields.size() >= 4) ? strip(fields[3]) : "");
		entries.push_back(entry);
	}

	return entries;

}

BedTarget::BedTarget(const string& s)
{
	parseRegion(s, seq, left, right);
}

vector<BedTarget*> BedReader::targetsContained(const BedTarget& target)
{
	vector<Interval<size_t, BedTarget*> > results = intervals[target.seq].findContained(target.left, target.right);
	vector<BedTarget*> contained;
	contained.reserve(results.size());
	for (const auto& res : results) {
		contained.push_back(res.value);
	}
	return contained;
}

vector<BedTarget*> BedReader::targetsOverlapping(const BedTarget& target)
{
	vector<Interval<size_t, BedTarget*> > results = intervals[target.seq].findOverlapping(target.left, target.right);
	vector<BedTarget*> overlapping;
	overlapping.reserve(results.size());
	for (const auto& r : results) {
		overlapping.push_back(r.value);
	}
	return overlapping;
}

void BedReader::addTargets(vector<BedTarget>& targets)
{
	map<string, vector<Interval<size_t, BedTarget*> > > intervalsBySeq;
	for (vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t) {
		intervalsBySeq[t->seq].push_back(Interval<size_t, BedTarget*>(1 + t->left, t->right, &*t));
	}
	for (map<string, vector<Interval<size_t, BedTarget*> > >::iterator s = intervalsBySeq.begin(); s != intervalsBySeq.end(); ++s) {
		intervals[s->first] = IntervalTree<size_t, BedTarget*>((vector<Interval<size_t, BedTarget*> >&&)s->second);
	}
}

void BedReader::open(const string& fname)
{
	file.open(fname.c_str());
	_isOpen = true;
	targets = entries();
	map<string, vector<Interval<size_t, BedTarget*> > > intervalsBySeq;
	for (vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t) {
		intervalsBySeq[t->seq].push_back(Interval<size_t, BedTarget*>(1 + t->left, t->right, &*t));
	}
	for (map<string, vector<Interval<size_t, BedTarget*> > >::iterator s = intervalsBySeq.begin(); s != intervalsBySeq.end(); ++s) {
		intervals[s->first] = IntervalTree<size_t, BedTarget*>((vector<Interval<size_t, BedTarget*> >&&)s->second);
	}
}
