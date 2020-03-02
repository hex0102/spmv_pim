import numpy as np
import sys
import scipy.sparse
from scipy.sparse import *
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import concurrent.futures


def collision_check(group_id):
    j = group_id
    col_group = packing_group[j]
    col_group_rowidx = packing_group_entries[j]
    confliced = 1
    if len(np.intersect1d(current_col_ind, col_group_rowidx)) > 0:
        return 1, 0
    else:
        return 0, len(col_group)


# inputs: a csc matrix:
# outputs: column groups:
def column_packing(m_input):
    if not isinstance(m_input, csc_matrix):
        m_input = csc_matrix(m_input) # m.indices  m.data m.indptr


    for i in range(len(m_input.indptr)-1):

        if i%1000 == 0:
            print(i)
            print(len(packing_group))
        global current_col_ind
        current_col_ind = m_input.indices[m_input.indptr[i]: m_input.indptr[i+1]]
        flag = 0
        if len(packing_group) == 0:
            packing_group.append([i])
            packing_group_entries.append(current_col_ind)
            continue

        candidate_groups = [] # first column (group id) | second column (n_cols in the group)

        n_group = len(packing_group)
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for out1, out2 in executor.map(collision_check, range(n_group)):
                candidate_groups.append([out1, out2])
                pass

        candidate_groups = np.asarray(candidate_groups)

        if np.min(candidate_groups, axis=0)[0] == 1:
            packing_group.append([i])
            packing_group_entries.append(current_col_ind)
        else:
            possible_group_idx = np.where(candidate_groups[:, 0] == 0)
            len_per_groups = candidate_groups[possible_group_idx, 1]
            selected_idx = possible_group_idx[0][np.argmin(len_per_groups)]
            if selected_idx.size > 1:
                selected_idx = np.random.choice(selected_idx.size, 1, replace=False)
            packing_group[selected_idx].append(i)
            packing_group_entries[selected_idx] = np.concatenate((packing_group_entries[selected_idx], current_col_ind))


        '''
        for j in range(len(packing_group)):
            col_group = packing_group[j]
            col_group_rowidx = packing_group_entries[j]
            # m_input.indices[m_input.indptr[col_group[0]]: m_input.indptr[col_group[0]+1]]
            # for m in range(1, len(col_group)):
            #   col_group_rowidx = np.concatenate((col_group_rowidx, m_input.indices[m_input.indptr[col_group[m]]: m_input.indptr[col_group[m]+1]]))
            if len(np.intersect1d(current_col_ind, col_group_rowidx)) > 0:
                continue
            else:
                flag = 1
                col_group.append(i)
                packing_group_entries[j] = np.concatenate((col_group_rowidx, current_col_ind))
                break
        if flag == 0:
            packing_group.append([i])
            packing_group_entries.append(current_col_ind)
        '''




def fix_hist_step_vertical_line_at_end(ax):
    axpolygons = [poly for poly in ax.get_children() if isinstance(poly, mpl.patches.Polygon)]
    for poly in axpolygons:
        poly.set_xy(poly.get_xy()[:-1])


packing_group = []
packing_group_entries = []
if __name__ == "__main__":

    #fig, ax = plt.subplots(figsize=(8, 4))

    if len(sys.argv) == 3:
        # snap file
        filename = sys.argv[1]
        n_edge = int(sys.argv[2])
    else:
        print("Incorrect number of argv...")
        exit(0)

    output_filename = filename[0:-4] + '.stat'
    coo_list_t = np.zeros((n_edge, 2))
    with open(filename) as f:
        curr_line = f.readline()
        cc = 0
        while(curr_line):
            if "#" in curr_line:
                curr_line = f.readline()
                continue
            list_num = curr_line.split()
            coo_list_t[cc, 0] = int(list_num[1])
            coo_list_t[cc, 1] = int(list_num[0])
            cc += 1
            curr_line = f.readline()

    coo_list_t = np.unique(coo_list_t, axis=0)
    n_node_total = np.max(coo_list_t)
    data = np.ones_like((coo_list_t[:, 0])).astype(int)
    unique_nodes = np.sort(np.unique(coo_list_t))
    n_node_used = unique_nodes.size

    dictionary = dict(zip(unique_nodes, np.arange(unique_nodes.size)))
    coo_list_t=np.vectorize(dictionary.get)(coo_list_t)

    spmv_coo = scipy.sparse.coo_matrix((data, (coo_list_t[:,0], coo_list_t[:,1])))
    spmv_csr = spmv_coo.tocsr()
    indptr = spmv_csr.indptr
    entry_n_row = indptr[1:] - indptr[:-1]

    spmv_csc = spmv_coo.tocsc()
    indptr = spmv_csc.indptr
    entry_n_col = indptr[1:] - indptr[:-1]

    for q in [50, 90, 95, 100]:
        print("{}%% percentile: {}".format(q, np.percentile(entry_n_row, q)))
    for q in [50, 90, 95, 100]:
        print("{}%% percentile: {}".format(q, np.percentile(entry_n_col, q)))

    fig, ax = plt.subplots(1, 2, sharey=True, tight_layout=True)

    df1 = pd.DataFrame(entry_n_row)
    df2 = pd.DataFrame(entry_n_col)
    ax[0] = df1.plot.hist(ax=ax[0], bins=1000, cumulative=True, density=True, histtype='step', legend=False)
    ax[1] = df2.plot.hist(ax=ax[1], bins=1000, cumulative=True, density=True, histtype='step', legend=False)
    fix_hist_step_vertical_line_at_end(ax[0])
    fix_hist_step_vertical_line_at_end(ax[1])
    plt.xlabel('row/col size')
    plt.show()


    column_packing(spmv_csc)

    with open(output_filename, 'w+') as f:
        for i in range(len(packing_group)):
            current_group_size = len(packing_group[i])
            f.write(str(current_group_size)+"\n")




    pass
