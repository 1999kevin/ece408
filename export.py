#!/usr/bin/python
# -*- coding: utf-8 -*-

import csv
import numpy as np
from scipy.sparse import csr_matrix


#csv格式：z-元素size，m-num_row，n-num_col， x-data
#z,m, n
#x
def export_full_r(filename, x):
    with open(filename, 'w', newline='') as csvfile:
        m, n = x.shape
        # z = np.array([np.dtype(np.double).itemsize], dtype='int32')
        writer = csv.writer(csvfile)
        writer.writerow([m,n])
        for i in range(m):
            writer.writerow(x[i])

#csv格式：
#N，m,n
#v
#j
#r
def export_rcs(filename, A):
    with open(filename,"w", newline='') as csvfile:
        A = csr_matrix(A, dtype=np.double)

        m, n = A.shape
        N = A.nnz
        v = A.data
        j = A.indices
        r = A.indptr

        writer = csv.writer(csvfile)

        #先写入columns_name
        writer.writerow([N,m,n])
        #写入多行用writerows
        writer.writerows([v,j,r])

def export_coo(filename, A):
    with open(filename,"w", newline='') as csvfile:
        A = csr_matrix(A, dtype=np.double)

        m, n = A.shape
        N = A.nnz
        v = A.data
        j = A.indices
        r = A.indptr

        r_ = np.zeros(len(j))
        number_row = 0
        r_index = 0
        i=0
        while(i<len(j)):
            if(i<=r[r_index+1]-1):
                r_[i] = number_row
            else:
                number_row+=1
                r_index+=1
                r_[i] = number_row
            i+=1

        writer = csv.writer(csvfile)

        #先写入columns_name
        writer.writerow([N,m,n])
        #写入多行用writerows
        writer.writerows([v,j,r_])

#z-total number of element,
#rank
#n_phy
#n-
#row -非零data
#csv格式：
#itemsize,rank,n_phy,n
#data
def export_sb_toe_r(filename, row):
    rank = len(blocks.shape)

    # only support 1-D filter for now --- more work for general case
    assert rank == 1

    n = row.shape
    n_phy = list(map(len, blocks.nonzero()))
    with open(filename,"w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([np.dtype(np.double).itemsize,rank,n_phy[0],n[0]])

        z = np.array(row, dtype=np.double)
        writer.writerow(z)

if __name__ == '__main__':
    N = 100
    L = 10
    M = 20

    #res = N

    #P_HT_fname = 'P_HT.csv'
    e_fname = 'data/e.csv'
    H_fname = 'data/H.csv'
    C_fname = 'data/C.csv'

    e = np.random.randn(N, L)
    export_full_r(e_fname, e)

    #H = np.random.randn(M, N)
    import numpy as np
    import scipy.sparse as ss
    import random
# 随机产生行、列坐标和值
    select_list1 = range(0,20)
    a = []
    for i in range(200):
        a.append(random.choice(select_list1))
    select_list2 = range(0,100)
    b = []
    for i in range(200):
        b.append(random.choice(select_list2))
    c = np.random.randn(200)
#a = random.sample(range(0, 19), 20)
#b = random.sample(range(0, 99), 200)
#c = random.sample(range(0, 1999), 200)
# 将list数据转为array数组
    rows = np.array(a)
    cols = np.array(b)
    v = np.array(c)
# coo_matrix函数生成稀疏矩阵
    sparseM = ss.coo_matrix((v,(rows,cols)))
# todense将稀疏矩阵转为完全阵
    H = sparseM.todense()
    export_coo(H_fname, H)

#row0 = np.ones((N, ))
    row0 = np.hstack((np.arange(50, 0, -1), np.zeros((50,))))
    blocks = row0
    export_sb_toe_r(C_fname, blocks)

#C = toeplitz(row0)
