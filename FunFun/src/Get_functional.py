def get_functionality(fungi, matrix_possitions, kofam_ontology, n_neigbors, epsilont=0.5):
    """
    The functuion functional content of an individual fungus.
    ----------
    fungi : Bio.Seq.Seq.id
        ITS ID.
    matrix_possitions : DataFrame
        Kmers counts dictionary.
    kofam_ontology : DataFrame
        Representation of group ortology for every fungi in base.
    n_neigbors : int
        K nearest neighbours.
    epsilont : float
        ε neighborhood.
    Returns
    -------
    dict_of_methabolic_function : dict
         Representation of group ortology for target fungi.
    """
    neighbor = matrix_possitions[fungi].sort_values()
    neighbor = neighbor[neighbor <= epsilont]

    if len(neighbor) > n_neigbors:
        
        neighbor = neighbor[: n_neigbors]

    if len(neighbor[neighbor == 0].index) > 0:
        
        dict_of_methabolic_function = kofam_ontology[neighbor[neighbor == 0].index].mean(axis=1)
        
    else:
        
        dict_of_methabolic_function = kofam_ontology[neighbor.index].mean(axis=1)

    return dict_of_methabolic_function