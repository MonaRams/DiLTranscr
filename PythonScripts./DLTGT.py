print(__doc__)

from time import time

import matplotlib as m

#m.use('Agg')  # to not to show plots
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import re
import os
import warnings
import random
import copy
from numpy import array
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn import metrics
#DLT
from sklearn.decomposition import MiniBatchDictionaryLearning
from sklearn.decomposition import DictionaryLearning
from sklearn.metrics import confusion_matrix
# ICA and other methods
from sklearn.datasets import load_digits
from sklearn.decomposition import FastICA
from sklearn.decomposition import NMF
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

#other
import time
from datetime import datetime
import operator

#############################################################################
#fix directories
os.chdir("..")
mainDr = os.getcwd()
os.chdir("Data")
dataDir = os.getcwd()

#############################################################################
# get file names of all prepped files
fls = os.listdir()
dataFls = []
for fl in fls:
    if "DatPREPPED" in fl:#) and not("_geneNames" in fl) and not("_typeCol" in fl) and not("normalised" in fl) and ".txt" in fl:# and "847" in fl
        dataFls.append(fl)

#############################################################################
# prepare empty tables (...)
allPath_allKmeansResults = {}  # list()
allPath_allKmeansARI = pd.DataFrame(
    columns=['RS', 'Natoms', 'Nsparse',  'ARI', 'MI', 'RecError'])  # list()
# set global parameters
mthds = ['DLT', 'PCA', 'ICA', 'NMF', 'tSNE', 'UMAP']
posDict = False
posCode = False
alphaParam = 1
plotSummary = True
predefinedNrAt = []
allRs = np.array(range(0, 10))
######################################################################################
######################################################################################
Xs = {}
#time.sleep(3600)
for path in dataFls:
    ################################################################
    # data upload
    ################################################################
    os.chdir(dataDir)
    if not '.txt' in path or '.csv' in path:
        drs = os.listdir()
        ptrn = path + '.txt'
        path = [i for i in drs if ptrn in i][0]

    print('Loading data ' + path + ' ...')

    name = re.sub('.txt|.csv|PREPPED|Dat', '', path)

    X = pd.read_csv(path, sep="\t", header=0)
    info = pd.read_csv(name + 'SampleInfoPREPPED.txt', sep="\t", header=0)
    infoCols = pd.read_csv(name + 'TypeColsPREPPED.txt', sep="\t", header=0)
    print('done')

    infoColA = infoCols.values.flatten()[0]
    true_labels = pd.Categorical(info[infoColA], categories=info[infoColA].unique()).codes
    nGroups = len(np.unique(info[infoColA])) # OR read in

    print("Automated dictCol")

    dictCols = np.array(range(1, 2 * nGroups + 1))

    if max(dictCols)<4:
        dictCols = np.array(range(1, 5))

    if X.shape[0] > X.shape[1]:
        X = X.T
        print('transposed input')
    Xs[name] = X

    if not ('alldictCols' in locals()):
        alldictCols = {}
    alldictCols[name] = dictCols
    if not ('allTrueLabels' in locals()):
        allTrueLabels = {}
    allTrueLabels[name] = true_labels

    if not ('allNGroups' in locals()):
        allNGroups = {}
    allNGroups[name] = nGroups

    for mthd in mthds:  # ['ICA', 'NMF', 'PCA']:
        print('#################################################################################')
        print('#################################################################################')
        print('#################################################################################')
        print('Current method is ' + mthd)
        print('#################################################################################')
        print('#################################################################################')
        print('#################################################################################')

        X = Xs[name]
        ##############################
        # 2. normalisation ###
        ##############################
        print('Normalising...')
        if mthd=='NMF':
            X = X.values
            mn = X.min(axis=1, keepdims=True)
            mx = X.max(axis=1, keepdims=True)
            X = (X - mn) / (mx - mn)
        else:
            # standardise sample-wise (sum1)
            sms = np.sum(X, axis=1)
            X = X.T
            X /= sms
            X = X.T
            # standardise gene-wise (c&s)
            sd = np.std(X, axis=0)
            X -= np.mean(X, axis=0)
            X /= sd
            if mthd == list(set(mthds)-set(['NMF']))[0]:
                print('Writing to file...')
                X.to_csv('normalisedData' + '_' + name + '.txt', sep='\t', header=True)
            X = X.values
        print('done')


        os.chdir(mainDr)
        if not mthd + 'Results' in os.listdir():
            os.mkdir(mthd + 'Results')
        os.chdir(mthd + 'Results')
        resDrMthd = os.getcwd()

        allKmeansResults = {}  # list()
        allKmeansARI = pd.DataFrame(
            columns=['RS', 'Natoms', 'Nsparse',  'ARI', 'MI', 'RecError'])

        #make directory name and define transform algorithm
        if mthd == 'DLT':
            if posCode:
                transfAlg = 'lars'
            else:
                # transfAlg = 'omp'
                # transfAlg = 'lasso_cd'
                transfAlg = 'omp'  # subDir #lars

            if posDict:
                drNam = name + '_' + str(alphaParam) + '_Pdict' + '_' + transfAlg
            else:
                if posCode:
                    drNam = name + '_' + str(alphaParam) + '_Pcode' + '_' + transfAlg
                else:
                    drNam = name + '_' + str(alphaParam) + '_' + transfAlg
        else:
            if mthd == 'ICA':
                drNam = name + '_ica'
            else:
                if mthd == 'NMF':
                    drNam = name + '_nmf'
                else:
                    if mthd == 'PCA':
                        drNam = name + '_pca'
                    else:
                        if mthd == 'tSNE':
                            drNam = name + '_tsne'
                        else:
                            if mthd == 'UMAP':
                                drNam = name + '_umap'

        if not os.path.exists(drNam):
            os.makedirs(drNam)
        os.chdir(drNam)
        ####################################################################################################
        ####################################################################################################
        #computations for all different rarndom seeds
        for rs in allRs:
            print('#####################################################')
            print('#####################################################')
            print('Random state is ' + str(rs) + ' OF ' + str(allRs))
            print('#####################################################')
            print('#####################################################')
            maxARI = 0
            maxARIShort = 0
            MBD = False #will be overwritten as True in case a minibatch dictionary is learned
            dictCols = alldictCols[name]
            dictCols = np.array(dictCols)
            dictCols = dictCols[dictCols < X.shape[0]]
            # if mthd == 'tSNE':
            #    dictCols = dictCols[dictCols < 4]
            if mthd == 'UMAP':
                dictCols = dictCols[dictCols > 1]

            for col in dictCols:
            #for mthd in mthds:
                print('####################')
                print('Col = ' + str(col))
                print('####################')
                ##################################################
                ##################################################
                # learning the dictionary
                ##################################################
                ##################################################
                t0 = time.time()

                height, width = X.shape
                if mthd == 'DLT':
                    if height * width > 50000000 or col > 100:  # 500000000
                        print('#####################################################')
                        print('Learning the mini batch dictionary with ' + str(col) + ' atoms...')
                        if posDict:
                            dico = MiniBatchDictionaryLearning(n_components=col,
                                                               alpha=alphaParam,
                                                               n_iter=500,
                                                               random_state=rs,  # )
                                                               positive_dict=True,
                                                               transform_algorithm=transfAlg)
                            # , transform_algorithm = 'lars')positive_code=True,
                        else:
                            if posCode:
                                dico = MiniBatchDictionaryLearning(n_components=col, alpha=alphaParam,
                                                                   n_iter=500,
                                                                   random_state=rs,  # )
                                                                   positive_code=True,
                                                                   transform_algorithm=transfAlg)
                            else:
                                print('MINBATCH POS NG DICT')
                                dico = MiniBatchDictionaryLearning(n_components=col,
                                                                   alpha=alphaParam,
                                                                   n_iter=500,
                                                                   random_state=rs,
                                                                   transform_algorithm=transfAlg)

                        thisName = name + '_MBD' + str(col)
                        MBD = True
                    else:
                        print('#####################################################')
                        print('Learning the dictionary with ' + str(col) + ' atoms on entire data...')
                        t0 = time.time()
                        if posDict:
                            dico = DictionaryLearning(n_components=col, alpha=alphaParam,
                                                      positive_dict=True)  # , max_iter=50)  #
                            # #, transform_algorithm = 'lars')  # max_iter=1000positive_code=True,

                            thisName = name + '_D' + str(col)
                        else:
                            if posCode:
                                dico = DictionaryLearning(n_components=col,
                                                          alpha=alphaParam,
                                                          positive_code=True,
                                                          transform_algorithm='lars')
                                # #, transform_algorithm = 'lars')  # max_iter=1000positive_code=True,
                                # #, transform_algorithm = 'lars')  # max_iter=1000positive_code=True,
                            else:
                                dico = DictionaryLearning(n_components=col,alpha=alphaParam)
                                thisName = name + '_D' + str(col)

                        '''
                        Solves the optimization problem:

                        (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
                                    (U,V)
                                    with || V_k ||_2 = 1 for all  0 <= k < n_components

                        '''
                else:
                    if mthd == 'ICA':
                        ica = FastICA(n_components=col, random_state=rs)
                    else:
                        if mthd == 'NMF':
                            nmf = NMF(n_components=col, random_state=rs)
                        else:
                            if mthd == 'PCA':
                                pca = PCA(n_components=col)
                            else:
                                if mthd == 'tSNE':
                                    tsne = TSNE(n_components=2, random_state=rs, perplexity=10 * col)
                                else:
                                    if mthd == 'UMAP':
                                        ump = umap.UMAP(n_components=2, random_state=rs, n_neighbors=col)

                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                print("Time is " + current_time)


                ####################################################
                # ACTUAL DICTIONARY LEARNING
                ####################################################

                print('Computing the decomposition')
                ok = True
                if mthd == 'DLT':
                    with warnings.catch_warnings():
                        warnings.simplefilter('error')
                        try:
                            V = dico.fit(X).components_  # time consuming step
                        except Warning as e:
                            print('Warning during dictLearn:')
                            print(e)
                            ok = False  # True  #
                    if ok:
                        print("DLT worked")
                    else:
                        print("Error during dictionary learning")
                else:
                    if mthd == 'ICA':
                        try:
                            code = ica.fit_transform(X)
                        except:
                            print("Error ", sys.exc_info()[0], " in ICA occurred.")
                            ok = False
                        if ok:
                            V = ica.mixing_
                            V = V.T
                            sp = pd.DataFrame(code)
                    else:
                        if mthd == 'NMF':
                            code = nmf.fit_transform(X)
                            V = nmf.components_
                            sp = pd.DataFrame(code)
                        else:
                            if mthd == 'PCA':
                                code = pca.fit_transform(X)
                                V = np.matmul(X.T, code)
                                V = V.T
                                sp = pd.DataFrame(code)
                            else:
                                if mthd == 'tSNE':
                                    code = tsne.fit_transform(X)
                                    sp = pd.DataFrame(code)
                                else:
                                    if mthd == 'UMAP':
                                        code = ump.fit_transform(X)
                                        sp = pd.DataFrame(code)
                    sp.to_csv('sparse' + str(col) + '_' + name + '_RS' + str(rs) + '_' + mthd + '.txt', sep=' ')

                if ok:
                    dt = time.time() - t0
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    print('done in %.2f s ' % dt + '(' + str(round(dt / 60, 2)) + ' min)')
                    print("Time is " + current_time)


                    ##################################################
                    ##################################################
                    # Reconstruction using dictionary
                    ##################################################
                    ##################################################
                    #set number of atoms
                    if mthd != 'tSNE' and mthd != 'UMAP':
                        nrAts = np.array(range(1, V.shape[0] + 1))
                        # nrAts = np.array([14])
                        if len(predefinedNrAt) > 0:
                            nrAts = predefinedNrAt
                        nrAts = np.unique(nrAts[nrAts <= V.shape[0]])
                        nrAts = np.unique(nrAts[nrAts > 0])

                    sparseCode = {}

                    transform_algorithms = list()
                    if mthd == 'DLT':
                        for nrAt in nrAts:
                            toApp = tuple([str('Sparse Coding Algorithm\n' + str(nrAt) + ' atom'), transfAlg,
                                           # lasso_lars lasso_cd lars omp threshold
                                           {'transform_n_nonzero_coefs': nrAt}])
                            transform_algorithms.append(toApp)
                    else:
                        toApp = tuple([1, 1, 1])
                        transform_algorithms.append(toApp)

                    for title, transform_algorithm, kwargs in transform_algorithms:
                        if mthd == 'DLT':
                            print('\n' + title + '...')
                            # reconstructions[title] = X.copy()
                            t0 = time.time()
                            dico.set_params(transform_algorithm=transform_algorithm, **kwargs)
                            code = dico.transform(X)

                            sp = pd.DataFrame(code)
                            sp.to_csv('sparse' + str(col) + '_' + name + '_RS' + str(rs) + '_DiL.txt', sep=' ')

                            sparseCode[title] = code
                            rec = np.dot(code, V)
                            figTitle = name + '_' + transform_algorithm + str(
                                list(kwargs.values())[0]) + '_OrRecDiff' + '.png'
                            nSparse = kwargs['transform_n_nonzero_coefs']
                            dictKey = name + "_" + str(nSparse) + '_RS_' + str(rs)
                        else:
                            t0 = time.time()
                            sparseCode[name] = code
                            if mthd != 'tSNE' and mthd != 'UMAP':
                                rec = np.dot(code, V)
                            else:
                                rec = 0

                            figTitle = name + '_' + mthd + '.png'
                            dictKey = name + '_' + mthd + '_' + str(col) + '_RS_' + str(rs)
                            nSparse = col

                        if mthd != 'tSNE' or mthd != 'UMAP':
                            recError = np.linalg.norm(X - rec)
                            #totalErr = np.linalg.norm(X - rec) + alphaParam * np.linalg.norm(code, 1)
                            dt = time.time() - t0

                        computeCluster = True
                        if computeCluster:
                            print('Clustering..')
                            t0 = time.time()

                            kmeansRes = KMeans(n_clusters=len(np.unique(true_labels)), random_state=rs).fit(
                                code)

                            kmCl = {}
                            kmCl[' '] = KMeans(n_clusters=len(np.unique(true_labels)), random_state=rs).fit(
                                code).labels_
                            if not len(np.unique(true_labels)) == X.shape[0]:
                                kmCl['plus1'] = KMeans(n_clusters=len(np.unique(true_labels)) + 1,
                                                       random_state=rs).fit(
                                    code).labels_
                                kmCl['plus2'] = KMeans(n_clusters=len(np.unique(true_labels)) + 2,
                                                       random_state=rs).fit(
                                    code).labels_

                            # computing Adusted Rand INdex (ARI)
                            resARIs = []
                            for i in kmCl.keys():
                                thisARI = metrics.adjusted_rand_score(true_labels, kmCl[i])
                                thisMI = metrics.adjusted_mutual_info_score(true_labels, kmCl[i])
                                resARIs.append(thisARI)
                            if mthd != 'tSNE' and mthd != 'UMAP':
                                print('ARIS: ' + str(np.round(resARIs, 2)) + ' (' + str(V.shape[0]),
                                      'atoms, ' + mthd + ', ' + 'rs = ' + str(rs) + ')')
                            else:
                                print('ARIS: ' + str(np.round(resARIs, 2)) + ' (' + str(col),
                                      'atoms, ' + mthd + ', ' + 'rs = ' + str(rs) + ')')

                            thisARI = metrics.adjusted_rand_score(true_labels, kmCl[' '])
                            thisMI = metrics.adjusted_mutual_info_score(true_labels, kmCl[' '])

                            if thisARI >= maxARI:
                                maxARI = thisARI
                                print('\n--------------------- MAX ARI HERE (ARI) ---------------------\n')

                            # save Dict to file
                            if thisARI >= maxARI:
                                if len(allRs) <= 10 and mthd != 'tSNE' and mthd != 'UMAP':
                                    np.savetxt(
                                        name + '_' + mthd + '_' + str(col) + 'Comp_RS' + str(rs) + '.txt', V)

                            allKmeansResults[dictKey] = kmCl[' ']
                            newCol = np.round([[rs, col, nSparse,, thisARI,
                                                thisMI, recError]], 2)
                            allKmeansARI = allKmeansARI.append(
                                pd.DataFrame(newCol, columns=allKmeansARI.columns, index=[dictKey]))

                            # save to file
                            if len(allRs) == 1:
                                sp = pd.DataFrame(code)
                                sp.to_csv('sparse' + str(nSparse) + '_' + name + '_RS' + str(rs) + '_' + mthd + '.txt',
                                          sep=' ')
                                # "_ARIshort" + str(round(thisARIshort, 2)) + "_ARI" + str(round(thisARI, 2)) +

                            if len(Xs.keys()) < 20 and rs == 0 and plotSummary:
                                ##################################################
                                # plot eucl dist of sparse codes
                                ##################################################
                                from sklearn.metrics.pairwise import euclidean_distances

                                if mthd != 'DLT':
                                    thisName = name + '_' + mthd + '_' + str(int(col)) + '_' + str(int(nSparse))
                                if MBD:
                                    thisName = name + '_MBD_' + str(int(col)) + '_' + str(int(nSparse))
                                else:
                                    thisName = name + '_D_' + str(int(col)) + '_' + str(int(nSparse))

                                x = euclidean_distances(np.round(code, 4), np.round(code, 4))
                                plt.figure(figsize=(4 * 5, 4 * 3.3))
                                plt.imshow(x, vmin=x.min(), vmax=x.max(), cmap=plt.cm.gray,
                                           interpolation='nearest')
                                plt.savefig('eucl_' + name + '_' + str(nSparse) + '.png')
                                plt.close()
                        else:
                            if len(allRs) == 1:
                                sp = pd.DataFrame(code)
                                sp.to_csv('sparse' + str(nSparse) + '_' + name + '.txt', sep=' ')

            allKmeansARI['Natoms'] = allKmeansARI['Natoms'].astype(int)
            allKmeansARI['Nsparse'] = allKmeansARI['Nsparse'].astype(int)
            allKmeansARI['RS'] = allKmeansARI['RS'].astype(int)
            saveAllKmeansARI = allKmeansARI
            saveAllKmeansARI = round(saveAllKmeansARI, 3)
            saveAllKmeansARI.to_csv('ARIs_' + name + 'RS_' + str(min(saveAllKmeansARI['RS'])) + '.txt', sep='\t')