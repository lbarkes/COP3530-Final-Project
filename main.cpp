//
//  Project2
//

#include <iostream>
#include <vector>
#include <forward_list>
#include <unordered_map>
#include <queue>    // for std::priority_queue
#include <utility>  // for std::pair
#include <climits>
#include <cmath>    // for std::ceil
using namespace std;

class Realm {
private:
    string charm;
    int max_incantations;
    vector<int> magi_LIS_partial_sums; // LIS = Longest Increasing Subsquence
    void longest_inc_subseq(const vector<int> & X) {
        int N = X.size();
        vector<int> P(N,0);
        vector<int> M(N + 1,0);
        for (int i = 0; i <= N; i++) {
            M[i] = i;
        }
        int L = 0;
        for (int i = 0; i < N; i++) {
            int lo = 1, hi = L;
            while (lo <= hi) {
                int mid = ceil((lo + hi)/2.0);
                if (X[M[mid]] < X[i]) {
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }
            int new_l = lo;
            P[i] = M[new_l-1];
            M[new_l] = i;
            if (new_l > L) {
                L = new_l;
            }
        }
        vector<int> S(L,0);
        int k = M[L];
        for (int i = L - 1; i >= 0; i--) {
            S[i] = X[k];
            k = P[k];
        }
        vector<int> partial_sums(L,0);
        partial_sums[0] = S[0];
        for (int i = 1; i < S.size(); i++) {
            partial_sums[i] = partial_sums[i - 1] + S[i];
        }
        magi_LIS_partial_sums = partial_sums;
        max_incantations = L;
    }
public:
    Realm(string charm, vector<int> available_magi) {
        this->charm = charm;
        longest_inc_subseq(available_magi);
    }
    string get_charm() const {
        return charm;
    }
    int get_max_incantations() const {
        return max_incantations;
    }
    int get_gems(int magi_amnt) const {
        return magi_LIS_partial_sums[magi_amnt - 1];
    }
    int get_edit_distance(const Realm & other) const {
        string word1 = this->get_charm(), word2 = other.get_charm();
        int subcost;
        int m = word1.length(), n = word2.length();
        vector<vector<int> > table(m + 1, vector<int> (n + 1));
        for (int i = 1; i <= m; i++)
            table[i][0] = i; // setting the columns
        for (int j = 1; j <= n; j++)
            table[0][j] = j; //setting the rows
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (word1[i - 1] == word2[j - 1])
                    subcost=0;
                else
                    subcost=1;
                table[i][j] = min(table[i - 1][j - 1] + subcost,
                                  min(table[i][j - 1] + 1, table[i - 1][j] + 1));
                //this is finding the minimum of each action which is replace,insertion or delete
            }
        }
        return table[m][n];
    }
};

struct RealmEdge {
    string neighbor; // the charm of the adjacent realm
    int edit_dist;
    int gems;
    RealmEdge(string edge_connected, int edit_dist, int gems) {
        this->neighbor = edge_connected;
        this->edit_dist = edit_dist;
        this->gems = gems;
    }
};
class RealmGraph {
private:
    unordered_map<string, forward_list<RealmEdge>> buckets;
public:
    RealmGraph(const vector<Realm> & realms){
        for (int i = 0; i < realms.size(); i++) {
            forward_list<RealmEdge> neighbors;
            for (int j = 0; j < realms.size(); j++) {
                if (i == j) { // don't create self loops
                    continue;
                }
                int edit_dist = realms[i].get_edit_distance(realms[j]);
                int avail_incant = realms[i].get_max_incantations();
                if (avail_incant >= edit_dist) { // add as a neighbor only if reachable
                    int gem_amt = realms[i].get_gems(edit_dist);
                    neighbors.push_front(RealmEdge(realms[j].get_charm(), edit_dist, gem_amt));
                }
            }
            buckets[realms[i].get_charm()] = neighbors;
        }
    }
    const forward_list<RealmEdge> & get_neighbors(string realm) {
        return buckets[realm];
    }
    vector<string> get_realms() {
        vector<string> keys;
        keys.reserve(buckets.size());
        for (auto key_value : buckets) {
            keys.push_back(key_value.first);
        }
        return keys;
    }
};
bool comp(pair<string,int> a, pair<string,int> b) {
    return a.second > b.second;
}
void print_shortest_distance(RealmGraph & graph, string source, string destination) {
    // store the total distances and gems used while traversing the graph
    unordered_map<string, pair<int,int>> totals;
    // distance of infinity means edge hasn't been explored yet
    for (string realm : graph.get_realms()) {
        totals[realm] = make_pair(INT_MAX, INT_MAX);
    }
    totals[source] = make_pair(0, 0); // total distance & gems number from source node to itself should be 0
    priority_queue<pair<string,int>,vector<pair<string,int>>,decltype(&comp)> pq(&comp);
    pq.push(make_pair(source, 0));
    while(!pq.empty()) {
        pair<string,int> clst_realm = pq.top();
        pq.pop();
        string clst_realm_name = clst_realm.first;
        int clst_realm_dist = clst_realm.second;
        if (clst_realm_name == destination) { // if destination node popped from pq, we are done.
            break;
        }
        for (RealmEdge edge : graph.get_neighbors(clst_realm_name)) {
            string crnt_nbr = edge.neighbor;
            int crnt_nbr_dist = edge.edit_dist;
            if (totals[crnt_nbr].first >= clst_realm_dist + crnt_nbr_dist) {
                // add neighbor to pq only if its total distance is smaller than previously considered
                if (totals[crnt_nbr].first > clst_realm_dist + crnt_nbr_dist) {
                    pq.push(make_pair(crnt_nbr, clst_realm_dist + crnt_nbr_dist));
                }
                totals[crnt_nbr].first = clst_realm_dist + crnt_nbr_dist; // update best distance
                int crnt_gem_amnt = totals[crnt_nbr].second;
                if (crnt_gem_amnt > totals[clst_realm_name].second + edge.gems) { //check if gem usage improved
                    totals[crnt_nbr].second = totals[clst_realm_name].second + edge.gems; // store new min gem
                }
            }
        }
    }
    if (totals[destination].first != INT_MAX) {
        cout << totals[destination].first << " " << totals[destination].second << endl;
    } else {
        cout << "IMPOSSIBLE" << endl;
    }
}

int main(int argc, const char * argv[]) {
    int realm_num;
    cin >> realm_num;
    vector<Realm> realms; // stores all the realms
    realms.reserve(realm_num);
    for (int i = 0; i < realm_num; i++) {
        string charm;
        cin >> charm;
        int magi_num;
        cin >> magi_num;
        vector<int> total_magi(magi_num);
        for (int j = 0; j < magi_num; j++) {
            cin >> total_magi[j];
        }
        Realm r(charm, total_magi);
        realms.push_back(r);
    }
    string source, destination;
    cin >> source >> destination;
    RealmGraph graph(realms);
    print_shortest_distance(graph, source, destination);
    print_shortest_distance(graph, destination, source);
    return 0;
}