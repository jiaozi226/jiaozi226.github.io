// Copyright 2014 Sina Inc. All rights reserved.
// Author: yanbing3@staff.sina.com.cn (Yan-Bing Bai)
// Description: Simulator

#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>

using std::set;
using std::vector;
using std::string;
using std::ifstream;
using std::cin;
using std::cout;
using std::endl;


// Individual structure
class Ind {
	public:
		Ind(int _id, int _pos, int _stra):id_(_id), pos_(_pos), stra_(_stra), payoff_(0.0f), payoff_sum_(0.0f), payoff_cnt_(0) {}
  public:
	  int id_;
		int pos_;
		int stra_;
		double payoff_;
		double payoff_sum_;
		int payoff_cnt_;
};

double rand_double();
int rand_int(int max);
int rand_step(int range);

/*
 * Positions: 0 ~ num_pos
 * Strategies: 0 ~ num_stra
 * v: probability to remain the same position as parent
 * u: probability to remain the same strategy as parent
*/
void help() {
  std::cout << "You should add command line params as the following order:" << endl;
	std::cout << "num_individuals num_positions num_strategies prob_remain_strategy prob_remain_position cycle_count strategy_matrix position_range";
	std::cout << "Example: 100 30 3 0.5 0.5 10000 mat.dat";
	std::cout << endl;
}
int main(int argc, char** argv) {
  // Parse args
	int num_ind = 0, num_pos = 0, num_stra = 0, cycle_count = 0, range = 0;
	double u = 0.5f, v = 0.5f;
	string mat_fname;
	if (argc < 9) {
	  help();
		return 0;
	}
	num_ind = atoi(argv[1]);
	num_pos = atoi(argv[2]);
	num_stra = atoi(argv[3]);
	assert(num_stra == 2);
	u = atof(argv[4]);
	v = atof(argv[5]);
	cycle_count = atoi(argv[6]);
	mat_fname = argv[7];
	range = atoi(argv[8]);

	// Build data structures
	// Individuals and pos
	vector<Ind> ind_arr;
	ind_arr.reserve(num_ind);
	vector<set<int> > pos_arr;
	pos_arr.resize(num_pos);
	for (int i = 0; i < num_ind; i++) {
		Ind tmp(i, rand_int(num_pos), rand_int(num_stra));
		ind_arr.push_back(tmp);
		pos_arr[tmp.pos_].insert(i);
	}
	// Load in strategy matrix
	vector<vector<double> > stra_matrix;
	ifstream mat_in(mat_fname.c_str());
	for (int i = 0; i < num_stra; i++) {
		vector<double> tmp;
	  for (int j = 0; j < num_stra; j++) {
		  double f;
			mat_in >> f;
			tmp.push_back(f);
		}
		stra_matrix.push_back(tmp);
	}

	// payoff results for each round
	vector<double> strategy_ratio(num_stra, 0.0f);

  // icc icd statistics
	double icc_sum = 0.0;
	double icd_sum = 0.0;
	// Begin cycling
	for (int k = 1; k <= cycle_count; k++) {
	  // clear all payoffs and do statistics
		int num_ind_strategy_c = 0;
		int num_ind_strategy_d = 0;
		for (int i = 0; i < num_ind; i++) {
			// clear payoffs
		  ind_arr[i].payoff_ = 0.0f;
		  ind_arr[i].payoff_sum_ = 0.0f;
			ind_arr[i].payoff_cnt_ = 0;
			// Do the statistics
		  if (ind_arr[i].stra_ == 0) {
				num_ind_strategy_c += 1;
			} else {
			  num_ind_strategy_d += 1;
			}
		}
		// check
		assert(num_ind_strategy_c + num_ind_strategy_d == num_ind);

		// icc and icd in current round
		int icc = 0;
		int icd = 0;
		// Compute payoff for each individual
		int check_ind_num = 0;
	  for (int pos = 0; pos < num_pos; pos++) {
			vector<int> id_list;
			id_list.reserve(pos_arr[pos].size());
			for (set<int>::iterator it = pos_arr[pos].begin(); it != pos_arr[pos].end(); ++it) {
			  id_list.push_back(*it);
			}
			check_ind_num += id_list.size();
			// For each pair of individuals
			for (int i = 0; i < id_list.size(); i++) {
				int idx1 = id_list[i];
			  for (int j = i+1; j < id_list.size(); j++) {
				  int idx2 = id_list[j];
					// ind_arr[idx1].payoff_sum_ += stra_matrix[ind_arr[idx1].stra_][ind_arr[idx2].stra_];
					// ind_arr[idx2].payoff_sum_ += stra_matrix[ind_arr[idx2].stra_][ind_arr[idx1].stra_];
					ind_arr[idx1].payoff_+= 1;
					ind_arr[idx2].payoff_+= 1;
					ind_arr[idx1].payoff_cnt_++;
					ind_arr[idx2].payoff_cnt_++;

					if (ind_arr[idx1].stra_ == 0 && ind_arr[idx2].stra_ == 0) {
					  icc += 2;
					}
					if (ind_arr[idx1].stra_ !=  ind_arr[idx2].stra_) {
					  icd += 1;
					}
				}
				if (ind_arr[idx1].stra_ == 0) {
				  icc += 1;
				}
			}
		}
		// Just for check
		assert(check_ind_num == num_ind);
		// std::cout << icc << " " << icd << std::endl;
		icc_sum += (icc * num_ind_strategy_d);
		icd_sum += (icd * num_ind_strategy_d);

		//
		vector<double> payoff_sum(num_stra, 0.0f);
		vector<double> strategy_cnt(num_stra, 0);
		// Average and stat the sum
		for (int i = 0; i < num_ind; i++) {
			Ind& ind = ind_arr[i];
			payoff_sum[ind.stra_] += ind.payoff_;
			// Count the number of individuals for each strategy
			strategy_cnt[ind.stra_] += 1;
		}
		// Save result
		double total_payoff = 0.0f;
		for (int i = 0; i < num_stra; i++) {
			total_payoff += payoff_sum[i];

			// Compute the ratio for each stategy and add to total ratio
			double ratio = (double)strategy_cnt[i] / num_ind;
			strategy_ratio[i] += ratio;
		}
		// Generate a new ind
		vector<double> prob_list(num_ind, 0.0f);
		for (int i = 0; i < num_ind; i++) {
		  prob_list[i] = ind_arr[i].payoff_ / total_payoff;
		}
		double toss = rand_double();
		int idx = 0;
		double prob_sum = 0.0f;
		for (int i =0 ; i < num_ind; i++) {
		  prob_sum += prob_list[i];
			if (prob_sum >= toss) {
			  idx = i;
				break;
			}
		}
		Ind new_ind = ind_arr[idx];
		if (rand_double() < v) {
		  int step = rand_step(range);
			int old_pos = new_ind.pos_;
			int new_pos = (num_pos + old_pos + step) % num_pos;
			new_ind.pos_ = new_pos;
		}
		if (rand_double() < u) {
		  new_ind.stra_ = rand_int(num_stra);
		}
		// Delete an old one
		int del_idx = rand_int(num_ind);
		if (del_idx < 0 || del_idx >= num_ind) {
		  std::cout << "del idx:" << del_idx;
		}
		assert(del_idx >= 0 && del_idx < num_ind);
		set<int>& id_set = pos_arr[ind_arr[del_idx].pos_];
		set<int>::iterator del_it = id_set.find(del_idx);
		assert(del_it != id_set.end());
		id_set.erase(del_it);
		new_ind.id_ = del_idx;
		ind_arr[del_idx] = new_ind;
		pos_arr[new_ind.pos_].insert(del_idx);

		for (int i = 0; i < num_ind; i++) {
		  assert(i == ind_arr[i].id_);
		}

		// output info if neccessary
	  if (k % 1000000 == 0) {
		  std::cout << "result up to round " << k << std::endl;
			// Output the final result
			double sum = 0.0;
			std::cout << "avg ratio for each strategy:";
			for (int i = 0; i < num_stra; i++) {
			  std::cout << strategy_ratio[i] / k<< " ";
				sum += (strategy_ratio[i]/k);
			}
			std::cout << " sum=" << sum << endl;

			// output icc sum and icd sum
  		std::cout << "avg (icc * num_d):" << icc_sum / k << endl;
			std::cout << "avg (icd * num_d):" << icd_sum / k << endl;
			std::cout << "(avg icc*num_d) / (avg icd*num_d) = " << icc_sum / icd_sum << endl;
			std::cout << std::endl;
			std::cout << std::endl;
		}
	}
	// Output the final result
	double sum = 0.0;
	std::cout << "avg ratio for each strategy:";
	for (int i = 0; i < num_stra; i++) {
	  strategy_ratio[i] /= cycle_count;
	   std::cout << strategy_ratio[i] << " ";
		sum += strategy_ratio[i];
	}
	std::cout << endl;

	// output icc sum and icd sum
  std::cout << "avg (icc * num_d):" << icc_sum / cycle_count << endl;
	std::cout << "avg (icd * num_d):" << icd_sum / cycle_count << endl;
	std::cout << "(avg icc*num_d) / (avg icd*num_d) = " << icc_sum / icd_sum << endl;
}
// Return a double in range [0,1)
double rand_double() {
  return (double)(rand()) / ((double)RAND_MAX +1 );
}
// Return a int in range [0,max)
int rand_int(int max) {
  return (int)(rand_double() * max) % max;
}
int rand_step(int range) {
  int step = rand_int(range+1);
	if (rand() % 2 == 0) {
	  return step;
	} else {
	  return 0 - step;
	}
}
