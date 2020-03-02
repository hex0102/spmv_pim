#include <cstdlib>
#include <bits/stdc++.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;

#define PRINT_FREQUENCY 1000

struct snap_entry{
    int row_id;
    int col_id;
};

struct candidate_entry{
    int conflicted;
    int size;
};

/*
std::vector<int> intersection(std::vector<int> &v1,
                                      std::vector<int> &v2){
    std::vector<int> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}
*/

//        for (auto i: current_col_ind)
//            std::cout << i << ' ';

int column_packing(vector<int>& ptr_array, vector<int>& ind_array, int n_unique, \
    vector< vector<int> >& packing_group, vector< vector<int> >&  packing_group_entries){

    cout<<">>>Start Column Packing..."<<endl;
    candidate_entry init_entry = {1, 0};

    for(int i=0; i != n_unique; i++){
        if(i%PRINT_FREQUENCY==0){
            cout<<">>Processing column #"<<i<<endl;
            cout<<">>Current total # of column groups is "<< packing_group.size()<<endl;
        }
        vector<int>   current_col_ind(&ind_array[ptr_array[i]], &ind_array[ptr_array[i+1]]);
        int flag = 0;
        vector<int> current_v(1, i);
        if(packing_group.size()==0){
            packing_group.push_back(current_v);
            packing_group_entries.push_back(current_col_ind);
            continue;
        }
        vector<int> candidate_groups(packing_group.size(), 0);
        int sum = 0;
        #pragma omp parallel for shared(candidate_groups, packing_group, packing_group_entries) reduction(+: sum)
        for(int j=0; j < packing_group.size(); j++){
            vector<int> col_group = packing_group[j];
            vector<int> col_group_rowidx = packing_group_entries[j];
            vector<int> intersect_v;
            std::set_intersection(current_col_ind.begin(),current_col_ind.end(), \
                col_group_rowidx.begin(),col_group_rowidx.end(), back_inserter(intersect_v));
            if(intersect_v.size() > 0){
                sum += 0;
                candidate_groups[j] = 100000000;
            }else{
                sum += 1;
                candidate_groups[j] = packing_group[j].size();
            }
        }




        // cout<<sum<<endl;
        if(sum == 0){
            packing_group.push_back(current_v);
            packing_group_entries.push_back(current_col_ind);
        }else{
            int min_idx = std::distance(candidate_groups.begin(), std::min_element(candidate_groups.begin(), candidate_groups.end()));
            packing_group[min_idx].push_back(i);
            std::vector<int> dst;
            std::merge(packing_group_entries[min_idx].begin(), packing_group_entries[min_idx].end(), current_col_ind.begin(), current_col_ind.end(), std::back_inserter(dst));
            packing_group_entries[min_idx] = dst;
        }

    }


    ofstream output_packing_file;

    output_packing_file.open("packing_result.log");
    if (!output_packing_file.is_open()) {
        cout << "Uable to open output file" << endl;
    }

    for (int i = 0; i < packing_group.size(); i++) {
        for (int j = 0; j < packing_group[i].size(); j++){
            cout << packing_group[i][j] << " ";
        }
        output_packing_file<<packing_group[i].size()<< " ";
        cout << endl;
    }
    output_packing_file<<endl;
    for (int i = 0; i < packing_group_entries.size(); i++) {
        cout << packing_group_entries[i].size() << " ";
        output_packing_file << packing_group_entries[i].size() << " ";
    }



}


// sort the coo into row first order
bool sort_rowfirst(const snap_entry &a, const snap_entry &b)
{
    if (a.row_id != b.row_id)
        return a.row_id < b.row_id;
    else{
        return a.col_id < b.col_id;
    }
}

// sort the coo into col first order
bool sort_colfirst(const snap_entry &a, const snap_entry &b)
{
    if (a.col_id != b.col_id)
        return a.col_id < b.col_id;
    else{
        return a.row_id < b.row_id;
    }
}



// split a string with a predefined delimiter
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}



int main(int argc, const char * argv[]){

    int n_edges;
    string filename;

    if (argc < 3){
        cout<<"Usage: ./main snap_txt n_edges"<<endl;
        exit(1);
    }else{
	    filename = argv[1];
        n_edges = atoi(argv[2]);
    }

    cout<<"./main "<< filename  <<"\t"<< n_edges <<endl;

    size_t pos = filename.find(".txt");

    string outfile = filename.substr(0, pos) + ".hpp2";
    cout << "Output to file: " << outfile << endl;
    snap_entry init_entry = {0, 0};
    // Transposed snap_coo.   e.g., (dest_id, src_id)
    vector<snap_entry> snap_coo_t(n_edges, init_entry);

    vector<int> ind_t(n_edges*2, 0);

    string line;
    ifstream input_snap_file;

    input_snap_file.open(filename);
    if (!input_snap_file.is_open()) {
        cout << "Uable to open input file" << endl;
    }
    else {
        cout<<">>>start reading actual graph data..."<<endl;
        int i = 0;
        while (! input_snap_file.eof()){
            getline(input_snap_file, line);
            if (line.find("#") != string::npos) {
                continue; //skip the header of the graph files
            }else{
		        if(line.length() > 1){

                    vector<string> v;
                    if (line.find("\t") != string::npos){
                        v = split(line, string("\t"));
                    }else if(line.find(" ") != string::npos){
                        v = split(line, string(" "));
                    }else{
                        cout<<"unknow seperating delimiter..."<<endl;
                        exit(1);
                    }
                    // snap_coo_t[i] = {stoi(v[1]), stoi(v[0])};
                    ind_t[2*i] = stoi(v[1]);
                    ind_t[2*i + 1] = stoi(v[0]);
                    // cout<<i<<"\t"<<snap_coo_t[i].row_id<<" "<<snap_coo_t[i].col_id<<"\t"<<endl;
			        i+=1;
		        }
            }
        }
	    cout << ">>> Read " << i << " total NNZs..." << endl;
    }
    input_snap_file.close();

    int n_unique = std::set<int>( ind_t.begin(), ind_t.end() ).size();
    std::cout << "Number of unique elements is "<< n_unique << std::endl;
    std::vector<int> reference_v(n_unique) ; // vector with n_unique ints.
    std::iota(std::begin(reference_v), std::end(reference_v), 0); // Fill with 0, 1, ..., n_unique-1.
    //cout << ">>> start " << reference_v[0]<<endl;

    vector<int> keys = ind_t;
    sort( keys.begin(), keys.end() );
    keys.erase( unique( keys.begin(), keys.end() ), keys.end() );

    assert( reference_v.size() == keys.size() );


    // create a map from two vectors
    std::map< int, int > mergedMap;
    std::transform( keys.begin(), keys.end(), reference_v.begin(), \
        std::inserter(mergedMap, mergedMap.end() ), std::make_pair<int const&,int const&> );

    for (int i = 0; i < ind_t.size(); ++i) {
        ind_t[i]=mergedMap[ind_t[i]];
    }

    for (int i = 0; i < ind_t.size()/2; ++i){
        snap_coo_t[i].row_id = ind_t[2*i];
        snap_coo_t[i].col_id = ind_t[2*i+1];
    }


    cout<<"sorting..."<<endl;
    std::sort(snap_coo_t.begin(), snap_coo_t.end(), sort_colfirst);


    // should also eliminate duplicated entries...
    auto comp = [] ( const snap_entry& lhs, const snap_entry& rhs ) {return (lhs.row_id == rhs.row_id)&&(lhs.col_id == rhs.col_id);};
    auto last = std::unique(snap_coo_t.begin(), snap_coo_t.end(),comp);
    snap_coo_t.erase(last, snap_coo_t.end());

    // cout<<snap_coo_t[0].row_id<<" "<<snap_coo_t[0].col_id<<endl;
    ofstream output_snap_file;

    vector<int> val_array;
    vector<int> ind_array;
    vector<int> ptr_array;


    output_snap_file.open(outfile);
    // write in CSC format
    if (!output_snap_file.is_open()) {
        cout << "Uable to open output file" << endl;
    }
    else {
        for (auto i = snap_coo_t.begin(); i != snap_coo_t.end(); ++i) {
            output_snap_file << "1 ";
            val_array.push_back(1);
        }
        output_snap_file << "\n";
        for (auto i = snap_coo_t.begin(); i != snap_coo_t.end(); ++i) {
            output_snap_file << i->row_id << " ";
            ind_array.push_back(i->row_id);
        }
        output_snap_file << "\n";
        output_snap_file << "0 ";
        ptr_array.push_back(0);

        int cur_ptr = 0;
        int cc = 0;
        for (int cur_col = 0; cur_col != n_unique; ++cur_col) {
            while( snap_coo_t[cc].col_id == cur_col){
                cur_ptr++;
                cc += 1;

            }
            ptr_array.push_back(cur_ptr);
            output_snap_file << cur_ptr << " ";
        }
        output_snap_file << "\n";
    }
    output_snap_file.close();


    vector< vector<int> > packing_group;
    vector< vector<int> > packing_group_entries;
    column_packing(ptr_array, ind_array, n_unique, packing_group, packing_group_entries);

    return 0;
}
