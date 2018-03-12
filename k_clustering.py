import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import copy
from data_arrange import ListToFormattedString
import  time


def k_clustering(filename):
    colmap = {1: 'r', 2: 'g', 3: 'b', 4:'c', 5:'m', 6:'y'}

    uq_pdvty = []
    allele=  []

    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if "V" in row['TYPE']:
                if float(row['UNIQUE_PRODUCTIVITY']) != 0.:
                    uq_pdvty.append(float(row['UNIQUE_PRODUCTIVITY']))
                    allele.append("".join([row['ALLELE'],row['FUNCTIONALITY']]))

        
    df = pd.DataFrame({'name': allele,
                       'x': uq_pdvty,
                       'y': [1]*len(uq_pdvty)
                       })
    np.random.seed(200)


    fig = plt.figure(figsize=(300, 4))
    ax = plt.axes()
    plt.scatter(df['x'], df['y'], linewidth=0.5, edgecolor='k',s=7,  color='k')
    plt.xlim(min(uq_pdvty), max(uq_pdvty))
    plt.show()

    error = []
    for k in range(1,7):
        # centroids[i] = [x, y]
        centroids = {
            i+1: [(max(uq_pdvty)-min(uq_pdvty))*np.random.random_sample()+min(uq_pdvty), 1]
            for i in range(k)
        }


        ## Assignment Stage
        def assignment(df, centroids):
            for i in centroids.keys():
                # sqrt((x1 - x2)^2 - (y1 - y2)^2)
                df['distance_from_{}'.format(i)] = (
                    np.sqrt(
                        (df['x'] - centroids[i][0]) ** 2
                    )
                )
            centroid_distance_cols = ['distance_from_{}'.format(i) for i in centroids.keys()]
            df['closest'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
            df['error']= df.loc[:, centroid_distance_cols].min(axis =1)
            df['closest'] = df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
            df['color'] = df['closest'].map(lambda x: colmap[x])
            return df
        df = assignment(df, centroids)

        ## Update Stage
        old_centroids = copy.deepcopy(centroids)
        def update(k):
            for i in centroids.keys():
                centroids[i][0] = np.mean(df[df['closest'] == i]['x'])
            return k
        centroids = update(centroids)
#        print  centroids

        ## Repeat Assigment Stage
        df = assignment(df, centroids)

        # Continue until all assigned categories don't change any more
        while True:
            closest_centroids = df['closest'].copy(deep=True)
            centroids = update(centroids)
            df = assignment(df, centroids)
            if closest_centroids.equals(df['closest']):
                break

        rank = {}
        sortedrank = sorted(centroids.values()[i][0] for i in range(len(centroids)))
        for nthplace in range(len(sortedrank)):
            prize = "level%d"%(nthplace+1)
            for n in centroids.keys():
                if sortedrank[nthplace]==centroids[n][0]:
                    rank[prize] = colmap[n]
                    
        #calculate error for k
        error.append(sum(df['error'].values.tolist()))

    plt.plot(range(1,7), error)
    plt.show()

                    
    k = input("elbow number: ")
    df = pd.DataFrame({'name': allele,
                       'x': uq_pdvty,
                       'y': [1]*len(uq_pdvty)
                       })
    np.random.seed(200)

    # centroids[i] = [x, y]
    centroids = {
        i+1: [(max(uq_pdvty)-min(uq_pdvty))*np.random.random_sample()+min(uq_pdvty), 1]
        for i in range(k)
    }


    ## Assignment Stage
    def assignment(df, centroids):
        for i in centroids.keys():
            # sqrt((x1 - x2)^2 - (y1 - y2)^2)
            df['distance_from_{}'.format(i)] = (
                np.sqrt(
                    (df['x'] - centroids[i][0]) ** 2
                )
            )
        centroid_distance_cols = ['distance_from_{}'.format(i) for i in centroids.keys()]
        df['closest'] = df.loc[:, centroid_distance_cols].idxmin(axis=1)
        df['error']= df.loc[:, centroid_distance_cols].min(axis =1)
        df['closest'] = df['closest'].map(lambda x: int(x.lstrip('distance_from_')))
        df['color'] = df['closest'].map(lambda x: colmap[x])
        return df

    df = assignment(df, centroids)
    #    print(df.head())

    ## Update Stage


    old_centroids = copy.deepcopy(centroids)

    def update(k):
        for i in centroids.keys():
            centroids[i][0] = np.mean(df[df['closest'] == i]['x'])
        return k

    centroids = update(centroids)

    ## Repeat Assigment Stage

    df = assignment(df, centroids)

    # Continue until all assigned categories don't change any more
    while True:
        closest_centroids = df['closest'].copy(deep=True)
        centroids = update(centroids)
        df = assignment(df, centroids)
        if closest_centroids.equals(df['closest']):
            break
    fig = plt.figure(figsize=(300, 4))
    ax = plt.axes()
    plt.scatter(df['x'], df['y'], color=df['color'], alpha=0.5, linewidth=0.5, edgecolor='k',s=7)
    for i in centroids.keys():
        plt.scatter(*centroids[i], color=colmap[i], s=7)
    plt.xlim(min(uq_pdvty), max(uq_pdvty))
    plt.show()

    rank = {}
    sortedrank = sorted(centroids.values()[i][0] for i in range(len(centroids)))
    for nthplace in range(len(sortedrank)):
        prize = "level%d"%(nthplace+1)
        for n in centroids.keys():
            if sortedrank[nthplace]==centroids[n][0]:
                rank[prize] = colmap[n]
                
    df.to_csv(''.join([filename[0:10],ListToFormattedString(rank.values()),'_out','.csv']))
