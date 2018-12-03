#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <string>

#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

// # cells
int num_cells, num_cells_actual;
// area of each cell
vector<int> cells_area;

// cell array, storing nets that are on this cell
vector<unordered_set<int>> cells;
// # nets
int num_nets;
// net array, storing cells that are on this net
vector<unordered_set<int>> nets;
// area of A or B < constraint
int area_constraint, area_diff_constraint;
// maximum possible gain
int p_max;

int total_area, A_area, B_area;
// two-way partition
unordered_set<int> A, B;
// block each cell belongs to
vector<int> cells_block;
// nets distribution, A=0, B=1
vector<vector<int>> nets_distribution;
// gain of each cell
vector<int> cells_gain;

// k-way partition
vector<unordered_set<int>> partitions;
// # partitions
int num_partitions;
// groups cells belong to
vector<int> cells_group;

// SA
int min_cost;
vector<unordered_set<int>> min_cost_partitions;

void ReadCellsFile(string filename)
{
    ifstream file;
    file.open(filename);

    file >> num_cells;
    cells_area = vector<int>(num_cells);

    total_area = 0;
    int id, area;
    unordered_set<int> init_partitions;
    for (int i = 0; i < num_cells; i++) {
        file >> id >> area;
        cells_area[id] = area;
        total_area += area;
        if (area == 0)
            num_cells_actual--;

        init_partitions.insert(i);
    }

    partitions.emplace_back(init_partitions);
    area_diff_constraint = total_area / 10;

    file.close();
}

void ReadNetsFile(string filename)
{
    ifstream file;
    file.open(filename);

    file >> num_nets;

    cells.insert(cells.begin(), num_cells, unordered_set<int>());
    nets.insert(nets.begin(), num_nets, unordered_set<int>());
    vector<int> temp(2, 0);

    for (int i = 0; i < num_nets; i++) {
        int num_pins, pin;
        file >> num_pins;
        for (int j = 0; j < num_pins; j++) {
            file >> pin;
            cells[pin].insert(i);
            nets[i].insert(pin);
        }

        nets_distribution.emplace_back(temp);
    }

    p_max = 0;
    for (int i = 0; i < num_cells; i++) {
        int p = cells[i].size();
        if (p > p_max)
            p_max = p;
    }

    file.close();
}

void ReadAuxFile(string aux_file)
{
    ifstream file;
    file.open(aux_file);

    string cells_file, nets_file;
    file >> cells_file >> nets_file >> area_constraint;

    file.close();

    string path = aux_file.substr(0, aux_file.find_last_of("/") + 1);
    ReadCellsFile(path + cells_file);
    ReadNetsFile(path + nets_file);
}

void InitialCut(const unordered_set<int> &partition)
{
    for (int i = 0; i < num_nets; i++) {
        nets_distribution[i][0] = 0;
        nets_distribution[i][1] = 0;
    }

    cells_block = vector<int>(num_cells, -1);
    A = unordered_set<int>();
    A_area = 0;

    int half_total_area = total_area / 2;
    bool done = false;
    for (const unordered_set<int> &net : nets) {
        for (int i : net) {
            if (partition.find(i) == partition.end() || cells_block[i] != -1)
                continue;

            if (A_area + cells_area[i] < half_total_area) {
                cells_block[i] = 0;
                A_area += cells_area[i];
                A.insert(i);
                for (int j : cells[i])
                    nets_distribution[j][0]++;
            }
            else {
                done = true;
                break;
            }
            if (done)
                break;
        }
    }

    B_area = total_area - A_area;

    B = unordered_set<int>();
    for (int i : partition) {
        if (cells_block[i] == -1) {
            B.insert(i);
            cells_block[i] = 1;
            for (int j : cells[i])
                nets_distribution[j][1]++;
        }
    }
}

void InitialSACut()
{
    for (int i = 0; i < num_nets; i++) {
        nets_distribution[i][0] = 0;
        nets_distribution[i][1] = 0;
    }

    A_area = 0;
    for (int i : A) {
        A_area += cells_area[i];
        for (int j : cells[i])
            nets_distribution[j][0]++;
    }

    cells_block = vector<int>(num_cells, 0);

    B_area = total_area - A_area;
    for (int i : B) {
        B_area += cells_area[i];
        for (int j : cells[i])
            nets_distribution[j][1]++;
    }
}

int InitialGain(const unordered_set<int> &partition, vector<unordered_set<int>> &gain_bucket)
{
    for (int i = 0; i < 2 * p_max + 1; i++)
        gain_bucket[i] = unordered_set<int>();

    cells_gain = vector<int>(num_cells);
    int max_gain = -1;
    for (int i : partition) {
        int g = 0;
        int block = A.find(i) != A.end() ? 0 : 1;

        for (int j : cells[i]) {
            // if F(n) = 1, g = g + 1
            if (nets_distribution[j][block] == 1)
                g++;
            // if T(n) = 0, g = g - 1
            if (nets_distribution[j][1 - block] == 0)
                g--;
        }

        // shift gain distribution to 0 ~ 2p_max + 1, no negative number
        g += p_max;
        if (g > max_gain)
            max_gain = g;
        gain_bucket[g].insert(i);
        cells_gain[i] = g;
    }

    return max_gain;
}

void UpdateGain(int base_cell, vector<int> &cells_lock_, 
                unordered_set<int> &A_, unordered_set<int> &B_, vector<vector<int>> &nets_distribution_,
                vector<unordered_set<int>> &gain_bucket_, vector<int> &cells_gain_)
{
    // Front block
    int F = A.find(base_cell) != A.end() ? 0 : 1;
    unordered_set<int> &F_block = F == 0 ? A_ : B_;
    // To block
    int T = 1 - F;
    unordered_set<int> &T_block = T == 1 ? B_ : A_;

    // lock the base cell and complement its block
    gain_bucket_[cells_gain_[base_cell]].erase(base_cell);
    cells_lock_[base_cell] = 1;
    F_block.erase(base_cell);
    T_block.insert(base_cell);

    for (int n : cells[base_cell]) {
        if (nets_distribution_[n][T] == 0) {
            for (int c : nets[n]) {
                if (cells_lock_[c] != 1) {
                    gain_bucket_[cells_gain_[c]++].erase(c);
                    gain_bucket_[cells_gain_[c]].insert(c);
                }
            }
        }
        else if (nets_distribution_[n][T] == 1) {
            for (int c : nets[n]) {
                if (T_block.find(c) != T_block.end() && cells_lock_[c] != 1) {
                    gain_bucket_[cells_gain_[c]--].erase(c);
                    gain_bucket_[cells_gain_[c]].insert(c);
                }
            }
        }

        nets_distribution_[n][F]--;
        nets_distribution_[n][T]++;

        if (nets_distribution_[n][F] == 0) {
            for (int c : nets[n]) {
                if (cells_lock_[c] != 1) {
                    gain_bucket_[cells_gain_[c]--].erase(c);
                    gain_bucket_[cells_gain_[c]].insert(c);
                }
            }
        }
        else if (nets_distribution_[n][F] == 1) {
            for (int c : nets[n]) {
                if (F_block.find(c) != F_block.end() && cells_lock_[c] != 1) {
                    gain_bucket_[cells_gain_[c]++].erase(c);
                    gain_bucket_[cells_gain_[c]].insert(c);
                }
            }
        }
    }
}

void Iteration(const unordered_set<int> &partition, int max_gain, vector<unordered_set<int>> &gain_bucket)
{
    // maximum partial sum of gains
    int G = 0;

    do {
        unordered_set<int> A_temp(A);
        unordered_set<int> B_temp(B);
        vector<vector<int>> nets_distribution_temp = nets_distribution;
        vector<unordered_set<int>> gain_bucket_temp = gain_bucket;
        vector<int> cells_gain_temp(cells_gain);
        // cell lock
        vector<int> cells_lock(num_cells, 1);
        // cells not belong to this partition are lock
        for (int i : partition)
            cells_lock[i] = 0;
        vector<int> cells_lock_temp = cells_lock;
        vector<int> cells_block_temp = cells_block;

        int A_area_temp = A_area, B_area_temp = B_area;

        // gain of each move
        vector<int> gains;
        // area difference after each move
        vector<int> area_diffs;
        // record move cells
        vector<int> moved_cells;
        // record max gains
        vector<int> max_gains;

        // total # cells moves
        int partition_size = partition.size();
        for (int i = 0; i < partition_size; i++) {
            int base_cell = -1;
            int area_diff;
            while (base_cell == -1) {
			    for (int j : gain_bucket_temp[max_gain]) {
                    if (cells_block_temp[j] == 0)
                        area_diff = abs((A_area_temp - cells_area[j]) - (B_area_temp + cells_area[j]));
                    else
                        area_diff = abs((A_area_temp + cells_area[j]) - (B_area_temp - cells_area[j]));

                    if (area_diff < area_diff_constraint) {
                        base_cell = j;
                        break;
                    }
                }

                if (max_gain == 0 || base_cell != -1)
                    break;

                max_gain--;
            }

            if (base_cell == -1)
                break;

            // shift back to original distribution (-pmax ~ pmax)
            gains.emplace_back(max_gain - p_max);
            area_diffs.emplace_back(area_diff);
            moved_cells.emplace_back(base_cell);

            UpdateGain(base_cell, cells_lock_temp, A_temp, B_temp, nets_distribution_temp, gain_bucket_temp, cells_gain_temp);

            if (cells_block_temp[base_cell] == 0) {
                A_area_temp -= cells_area[base_cell];
                B_area_temp += cells_area[base_cell];
            }
            else {
                A_area_temp += cells_area[base_cell];
                B_area_temp -= cells_area[base_cell];
            }
            cells_block_temp[base_cell] = 1 - cells_block_temp[base_cell];

            for (int j = p_max * 2; j >= 0; j--) {
                if (gain_bucket_temp[j].size() != 0) {
                    max_gain = j;
                    break;
                }
            }
        }

        G = -p_max - 1;
        int num_moves = -1, partial_sum = 0;
        int size = gains.size();
        for (int i = 0; i < size; i++) {
            partial_sum += gains[i];
            if (partial_sum > G) {
                G = partial_sum;
                num_moves = i;
            }
            else if (partial_sum == G) {
                if (area_diffs[i] < area_diffs[num_moves])
                    num_moves = i;
            }
        }
        num_moves++;

        if (G <= 0)
            break;

        for (int i = 0; i < num_moves; i++) {
            int base_cell = moved_cells[i];

            UpdateGain(base_cell, cells_lock, A, B, nets_distribution, gain_bucket, cells_gain);

            // complement base cell block and update A, B area
            if (cells_block[base_cell] == 0) {
                A_area -= cells_area[base_cell];
                B_area += cells_area[base_cell];
            }
            else {
                A_area += cells_area[base_cell];
                B_area -= cells_area[base_cell];
            }
            cells_block[base_cell] = 1 - cells_block[base_cell];
        }

        max_gain = InitialGain(partition, gain_bucket);
    } while (G > 0);
}

void RecursivePartition()
{
    bool done = true;
    do {
        done = true;
        vector<unordered_set<int>> partitions_temp;
        for (const unordered_set<int> &partition : partitions) {
            total_area = 0;
            for (const int cell : partition)
                total_area += cells_area[cell];

            if (total_area <= area_constraint) {
                partitions_temp.emplace_back(partition);
                continue;
            }

            done = false;
            area_diff_constraint = total_area / 10;

            InitialCut(partition);
            vector<unordered_set<int>> gain_bucket(p_max * 2 + 1);
            int max_gain = InitialGain(partition, gain_bucket);
            Iteration(partition, max_gain, gain_bucket);

            partitions_temp.emplace_back(A);
            partitions_temp.emplace_back(B);
        }
        partitions = partitions_temp;
    } while(!done);

    num_partitions = partitions.size();
}

int CalculateCost()
{
    cells_group = vector<int>(num_cells, 0);

    for (int i = 0; i < num_cells; i++) {
        int group = 0;
        for (int j = 0; j < num_partitions; j++) {
            if (partitions[j].find(i) != partitions[j].end()) {
                group = j;
                break;
            }
        }
        cells_group[i] = group;
    }

    int cost = 0;
    for (int i = 0; i < num_nets; i++) {
        unordered_set<int> span_group;
        for (int cell : nets[i])
            span_group.insert(cells_group[cell]);
        int span = span_group.size();
        cost += (span - 1) * (span - 1);
    }

    return cost;
}

void Output(string filename)
{
    ofstream file;
    file.open(filename);

    int cost = CalculateCost();
    file << cost << '\n';
    file << num_partitions << '\n';

    for (int i = 0; i < num_cells; i++)
        file << cells_group[i] << '\n';

    file.close();
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        cout << "Usage: \n";
        cout << "  ./k-way_min-cut_partition <path/to/aux_file>\n";
        exit(1);
    }
    string aux_file = argv[1];
    ReadAuxFile(aux_file);

    RecursivePartition();

    Output(aux_file + ".out");

    return 0;
}
