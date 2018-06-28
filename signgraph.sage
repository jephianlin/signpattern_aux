### a sign graph is a simple graph whose edges are labeled as "1" or "-1" to indicate their signs

### only sparse graph data structure allows edge label;
### use g=g.copy(sparse=True) to convert;


def rooted_spanning_tree(g,r=None):
    """
    Input:
        g: a simple graph;
        r: a vertex of g; default is g.vertices()[0];
    Output:
        an oriented spanning tree whose edges are directed to the root;
        (This function assume g is connected.)
    """
    if g.is_connected()==False:
        print "rooted_spanning_tree: The input graph should be connected."
        return None;
    V=g.vertices();
    n=len(V);
    if r==None:
        r=V[0];
    all_found=[r];
    front_wave=[r];
    to_be_added=[];
    st=DiGraph([V,[]]);
    ###Apply DFS
    while len(all_found)<n:
        for v in front_wave:
            for u in g.neighbors(v):
                if u not in to_be_added and u not in all_found:
                    st.add_edge(u,v);
                    all_found.append(u);
                    to_be_added.append(u);
        front_wave=to_be_added;
        to_be_added=[];
    return st;
    
def make_spanning_tree_positive(g,st=None,return_signature=False):
    """
    Input:
        g: a sign graph;
        st: a spanning tree of g;
        return_signature: if False, return the new sign graph h only;
            if True, return h, signature together;
    Output:
        a new sign graph which is sign similar to g,
        while all edges on st is positive;
    """
    if st==None:
        st=rooted_spanning_tree(g);
    roots=[v for v in st.vertices() if st.out_degree(v)==0];
    r=roots[0];
    tree_dict={u: v for u,v in st.edges(labels=False)};
    sign_signature={r:1};
    V=g.vertices();
    V.remove(r);
    ### get the sign_signature
    while len(V)>0:
        v=V.pop();
        path_v=[v];
        u=tree_dict[v];
        while u not in sign_signature.keys():
            path_v.append(u);
            V.remove(u);
            u=tree_dict[u];
        path_v.append(u);
        path_v.reverse();
        for i in range(1,len(path_v)):
            a=path_v[i-1];
            b=path_v[i];
            sign_signature[b]=sign_signature[a]*g.edge_label(a,b);
    #Diag=diagonal_matrix([sign_signature[v] for v in g.vertices()]);
    #A=g.weighted_adjacency_matrix();
    #B=Diag*A*Diag;
    h=g.copy();
    plus_side=[v for v in g.vertices() if sign_signature[v]==1];
    minus_side=[v for v in g.vertices() if sign_signature[v]==-1];
    for u in plus_side:
        for v in minus_side:
            if h.has_edge(u,v):
                h.set_edge_label(u,v,-h.edge_label(u,v)); 
    if return_signature:
        return h,sign_signature;
    else:
        return h;
    
def is_sign_bipartite(g,certificate=False):
    """
    Input:
        g: a sign graph;
        certificate: return a certificate or not;
    Output:
        True or False if g is a sign bipartite graph;
        certificate for True is a sign_signature (dictionary);
        certificate for False is a odd signed cycle (a list), which is a cycle with odd number of -1;
    """
    if g.is_connected():
        st=rooted_spanning_tree(g);
        simple_st=Graph(st);
        h,sig=make_spanning_tree_positive(g,st,True);
        for e in h.edges():
            if e[2]==-1:
                if certificate:
                    p=simple_st.all_paths(e[0],e[1])[0];
                    return p;
                else:
                    return False;
        if certificate:
            return sig;
        else:
            return True;
    else: #g has several components...
        for com in g.connected_components_subgraphs():
            if is_sign_bipartite(com)==False:
                if certificate:
                    return is_sign_bipartite(com,True);
                else:
                    return False;
        else:
            if certificate:
                sig={};
                for com in g.connected_components_subgraphs():
                    st=rooted_spanning_tree(com);
                    h,com_sig=make_spanning_tree_positive(com,st,True);
                    for key in com_sig.keys():
                        sig[key]=com_sig[key];
                return sig;
            else:
                return True;
