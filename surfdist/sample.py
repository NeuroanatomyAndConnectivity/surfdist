
def find_idx_match(simple_vertices, complex_vertices):
    '''
    Thanks to juhuntenburg. 
    Functions taken from https://github.com/juhuntenburg/brainsurfacescripts

    Finds those points on the complex mesh that correspoind best to the
    simple mesh while forcing a one-to-one mapping
    '''
    # make array for writing in final voronoi seed indices
    voronoi_seed_idx = np.zeros((simple_vertices.shape[0],), dtype='int64')-1
    missing = np.where(voronoi_seed_idx == -1)[0].shape[0]
    mapping_single = np.zeros_like(voronoi_seed_idx)

    neighbours = 0
    col = 0

    while missing != 0:

        neighbours += 100
        # find nearest neighbours
        inaccuracy, mapping = spatial.KDTree(
            complex_vertices).query(simple_vertices, k=neighbours)
        # go through columns of nearest neighbours until unique mapping is
        # achieved, if not before end of neighbours, extend number of
        # neighbours
        while col < neighbours:
            # find all missing voronoi seed indices
            missing_idx = np.where(voronoi_seed_idx == -1)[0]
            missing = missing_idx.shape[0]
            if missing == 0:
                break
            else:
                # for missing entries fill in next neighbour
                mapping_single[missing_idx] = np.copy(
                    mapping[missing_idx, col])
                # find unique values in mapping_single
                unique, double_idx = np.unique(
                    mapping_single, return_inverse=True)
                # empty voronoi seed index
                voronoi_seed_idx = np.zeros(
                    (simple_vertices.shape[0],), dtype='int64')-1
                # fill voronoi seed idx with unique values
                for u in range(unique.shape[0]):
                    # find the indices of this value in mapping
                    entries = np.where(double_idx == u)[0]
                    # set the first entry to the value
                    voronoi_seed_idx[entries[0]] = unique[u]
                # go to next column
                col += 1

    return voronoi_seed_idx, inaccuracy


def competetive_fast_marching(vertices, graph, seeds):
    '''
    Label all vertices on highres mesh to the closest seed vertex
    using a balanced binary search tree
    '''
    # make a labelling container to be filled with the search tree
    # first column are the vertex indices of the complex mesh
    # second column are the labels from the simple mesh
    # (-1 for all but the corresponding points for now)
    labels = np.zeros((vertices.shape[0], 2), dtype='int64')-1
    labels[:, 0] = range(vertices.shape[0])
    for i in range(seeds.shape[0]):
        labels[seeds[i]][1] = i
    # initiate AVLTree for binary search
    tree = FastAVLTree()
    # organisation of the tree will be
    # key: edge length; value: tuple of vertices (source, target)
    # add all neighbours of the voronoi seeds
    for v in seeds:
        add_neighbours(v, 0, graph, labels, tree)
    # Competetive fast marching starting from voronoi seeds
    while tree.count > 0:

        # pdb.set_trace()
        # pop the item with minimum edge length
        min_item = tree.pop_min()
        length = min_item[0]
        source = min_item[1][0]
        target = min_item[1][1]
        # if target no label yet (but source does!), assign label of source
        if labels[target][1] == -1:
            if labels[source][1] == -1:
                sys.exit('Source has no label, something went wrong!')
            else:
                # assign label of source to target
                labels[target][1] = labels[source][1]

        # test if labelling is complete
        if any(labels[:, 1] == -1):
            # if not, add neighbours of target to tree
            add_neighbours(target, length, graph, labels, tree)
        else:
            break

        # print 'tree '+str(tree.count)
        # print 'labels '+str(np.where(labels[:,1]==-1)[0].shape[0])

    return labels


def sample_simple(highres_data, labels):
    '''
    Computes the mean of data from highres mesh that is assigned to the same
    label (typical simple mesh vertices).
    '''
    # create new empty lowres data array
    lowres_data = np.empty((int(labels.max()+1), highres_data.shape[1]))
    # find all vertices on highres and mean
    for l in range(int(labels.max())):
        patch = np.where(labels == l)[0]
        patch_data = highres_data[patch]
        patch_mean = np.mean(patch_data, axis=0)
        lowres_data[l] = patch_mean

    return lowres_data
