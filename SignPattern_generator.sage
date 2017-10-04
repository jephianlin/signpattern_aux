### digraph with loop

def dig_transpose(g):
    """
    Input:
        g: a digraph (possibly with loop)
    Output:
        return the transpose of g;
        that is, the digraph obtained from g by reversing all arcs;
    """
    h=DiGraph(0,loops=g.allows_loops());
    for v in g.vertices():
        h.add_vertex(v);
    for e in g.edges(labels=False):
        h.add_edge(e[1],e[0]);
    return h;

def is_canonical_transpose(g):
    """
    Input:
        g: a digraph without loop;
    Output:
        compute the canonical dig6_strings of g and dig_transpose(g), say stg and stg_trans;
        return stg<=stg_trans;
    """
    stg=g.canonical_label().dig6_string();
    stg_trans=dig_transpose(g).canonical_label().dig6_string();
    return stg<=stg_trans;
    
def all_loop_configurations(g,trans=True,return_graph=False):
    """
    Input:
        g: a digraph without loop;
    Output:
        generate all non-isomorphic digraphs whose underlying graph is g;
        if return_graph is True, return the graph format; otherwise, only the list of vertices with loop.
    """
    n=g.order();
    V=g.vertices();
    if trans: #here really assumes g is labeled on range(n);
        ans, cert = g.is_isomorphic(dig_transpose(g),certificate=True);
        if ans: #if the graph is isomorphic to itself, it create more flexibility to do permutation on graph;
            trans_per=[cert[i] for i in range(n)];
            gens=[[per(i) for i in range(n)] for per in g.automorphism_group().gens()];
            gens.append(trans_per);
            PerG=PermutationGroup(gens, domain=range(n));
        else:
            PerG=list(g.automorphism_group());
    else:
        PerG=list(g.automorphism_group());
    for k in range(n+1):
        yielded={};
        for com in Combinations(V,k):
            min_com=com;
            min_stg=sum([2^i for i in com]);
            for per in PerG:
                new_com=[per(i) for i in com];
                new_stg=sum([2^i for i in new_com]);
                #print com,per,new_com,new_stg;
                if new_stg<min_stg:
                    min_stg=new_stg;
                    min_com=new_com;
            if min_stg not in yielded.keys():
                if return_graph:
                    h=g.copy();
                    h.allow_loops(True);
                    for i in min_com:
                        h.add_edge(i,i);
                    yield h;
                else:
                    yield min_com;
                yielded[min_stg]=min_com;
    
def digraphs_with_loop(n,irreducible=True,trans=True,no_loop=False):
    """
    Input:
        n: the number of vertices;
        irreducible: whether to generate irreducible (strongly connected) digraphs only;
        trans: whether to consider the transpose of a digraph and itself as the same;
        no_loop: mainly for debug purpose, if True, return cases without loop only.
    Output:
        generate all loop digraphs on n vertices;
    """
    for g in digraphs(n):
        if (irreducible and g.is_strongly_connected()) or irreducible==False:
            if (trans and is_canonical_transpose(g)) or trans==False:
                g=g.canonical_label();
                stg=g.dig6_string();
                if no_loop:
                    yield g;
                else:
                    for loop_g in all_loop_configurations(g,trans,True):
                        #yield stg, loop_g;
                        yield loop_g;

def tuple_generator(k,n):
    """
    Input:
        k: a positive integer at least ;
        n: a positive integer;
    Output:
        a generator generating [0,...,0], [1,0,...,0], to ,[k-1,...,k-1];
    """
    a=[0]*n;
    yield a;
    counter=1;
    max_counter=k^n;
    while counter<max_counter:
        a[0]+=1;
        for i in range(n):
            if a[i]>=k:
                a[i]-=k;
                a[i+1]+=1;
        yield a;
        counter+=1;

def sign_pattern_generator(n,irreducible=True,trans=True,neg=True):
    """
    Input:
        n: the size of the sign pattern;
        irreducible: whether to generate irreducible patterns only;
        trans: whether to conside the transpose of a pattern and itself as the same;
        neg: whether to consider the negation of a pattern and itself as the same;
    Output:
        generate all n by n sign patterns;
    """
    for dig in digraphs_with_loop(n,irreducible,trans):
        E=dig.edges(labels=False);
        m=len(E);
        appeared=[];
        if neg:
            for tup in tuple_generator(2,m-1):
                entries=[1]+[2*i-1 for i in tup];
                mtx=zero_matrix(n);
                for i in range(m):
                    mtx[E[i]]=entries[i];
                ptn=SignPattern(mtx);
                similar_class=ptn.sign_similar_class(trans,neg);
                for app_ptn in appeared:
                    if app_ptn in similar_class:
                        break;
                else:
                    ptn=ptn.canonical_label();
                    yield ptn;
                    appeared.append(ptn);
                
        else:
            print "not yet implemented";
