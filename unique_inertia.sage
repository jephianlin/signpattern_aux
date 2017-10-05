def det_sign(mtx,print_counter=False):
    """
    Input:
        mtx: a matrix representing the sign pattern;
        print_counter: whether to print the counters for positive and negative terms (mainly for debugging);
    Output:
        return the sign (1,-1, or 0) of the determinant of mtx if it is a fixed sign, otherwise return -2;
    """
    plus_counter=0;
    minus_counter=0;
    n=mtx.dimensions()[0];
    for per in SymmetricGroup(n):
        this_term=per.sign()*prod([mtx[i,per(i+1)-1] for i in range(n)]);
        if this_term>0:
            plus_counter+=1;
        if this_term<0:
            minus_counter+=1;
        if plus_counter>=1 and minus_counter>=1:
            return -2; #I am trying to avoid 0 treated as False;
    if print_counter:
        print plus_counter,minus_counter;
    if plus_counter>0:
        return 1;
    if minus_counter>0:
        return -1;
    if plus_counter==0 and minus_counter==0:
        return 0;

def is_sign_nonsingular(mtx):
    """
    Input:
        mtx: a matrix representing the sign pattern;
    Output:
        whether the matrix is sign nonsingular or not;
        sign nonsingular means the determinant has a fixed sign and is 1 or -1;
    """
    return det_sign(mtx) in [1,-1];
        
def minor_sum_sign(mtx,k,print_counter=False):
    """
    Input:
        mtx: a matrix representing the sign pattern;
        k: the size of the minors considered;
        print_counter: whether to print the counters for positive and negative terms (mainly for debugging);
    Output:
        return the sign of S_k (the sum all of all k minors);
        -1,0,1 represent the signs, while -2 means not fixed sign;        
    """
    if k==0:
        return 1;
    plus_counter=0;
    minus_counter=0;
    n=mtx.dimensions()[0];
    for com in Combinations(range(n),k):
        this_term=det_sign(mtx[com,com]);
        if this_term==-2:
            return -2
        if this_term==1:
            plus_counter+=1;
        if this_term==-1:
            minus_counter+=1;
        if plus_counter>0 and minus_counter>0:
            return -2;
    if print_counter:
        print plus_counter,minus_counter;
    if plus_counter>0:
        return 1;
    if minus_counter>0:
        return -1;
    if plus_counter==0 and minus_counter==0:
        return 0;
        
def S_sequence(mtx):
    """
    Input:
        mtx: a matrix representing the sign pattern;
    Output:
        a sequence of the signs of [S0,S1,...,Sn];
    """
    n=mtx.dimensions()[0];
    return [minor_sum_sign(mtx,i) for i in range(n+1)];
    
def sign_changes(l):
    """
    Input:
        l: a list;
    Output:
        return the number of sign changes of l, ignoring zeros;
    """
    n=len(l);
    counter=0;
    for k in range(n):
        if l[k]!=0:
            a=l[k];
            break;
    for i in range(k+1,n):
        b=l[i];
        if b!=0:
            if a*b<0:
                counter+=1;
            a=b;
    return counter; 

def odd_terms_negation(l):
    """
    Input:
        l: a list;
    Output:
        the list obtained from l where the odd terms are negated;
    """
    n=len(l);
    new_l=copy(l);
    for i in range(n):
        if i%2==1:
            new_l[i]*=-1;
    return new_l;

def is_skew_tree_pattern(mtx):
    new_mtx=SignPattern(mtx).repr;
    if new_mtx==-new_mtx.transpose():
        n=new_mtx.dimensions()[0];
        for i in range(n):
            for j in range(n):
                if new_mtx[i,j]!=0:
                    new_mtx[i,j]=1;
        g=Graph(new_mtx);
        return g.is_tree();      
    else:
        return False;

def is_combinatorial_symmetric(mtx,same_sign=True):
    """
    Input:
        mtx: a matrix representing the sign pattern;
        same_sign: whether sign symmetric or just zero-nonzero symmetric;
    Output:
        return whether mtx is combinatorial symmetric or not;
    """
    new_mtx=SignPattern(mtx).repr;
    if same_sign:
        return (new_mtx==new_mtx.transpose());
    else:
        m,n=new_mtx.dimensions();
        for i in range(n):
            for j in range(m):
                if new_mtx[i,j]!=0:
                    new_mtx[i,j]=1;
        return (new_mtx==new_mtx.transpose());
        
def is_tree_pattern(mtx,same_sign=True):
    """
    Input:
        mtx: a matrix representing the sign pattern;
        same_sign: whether sign symmetric or just zero-nonzero symmetric;
    Output:
        either return the simple graph (which is a tree) of mtx or return False;
    """
    if is_combinatorial_symmetric(mtx,same_sign):
        new_mtx=copy(mtx);
        n=new_mtx.dimensions()[0];
        diag=[sign(new_mtx[i,i]) for i in range(n)];
        for i in range(n):
            new_mtx[i,i]=0;
        for i in range(n):
            for j in range(n):
                if new_mtx[i,j]!=0:
                    new_mtx[i,j]=1;
        g=Graph(new_mtx);
        for i in range(n):
            g.set_vertex(i,diag[i]);
        if g.is_tree():
            return g;
        else:
            return False;
    else:
        return False;
        
def UI_tree_pattern(g,print_process=False):
    """
    Input:
        g: a graph, whose vertices are associated with the signs;
           that is, g.get_vertex(i) = the sign;
           off-diagonal entries can be assume as 1 (or all positive) by similarity;
        return_inertia: if True, return the inertia if unique;
        print_process: if True, print the eliminating process;
    Output:
        return whether the corresponding sign pattern requires a unique inertia or not.
    """
    if g.size()==0:
        return True;
    else:
        V=g.vertices();
        for v in V:
            if g.degree(v)==1:
                u=g.neighbors(v)[0];
                ### Rule 1
                if g.get_vertex(v)==0:
                    if print_process:
                        print [v,u],"(1,1,0)";
                    h=g.copy();
                    h.delete_vertex(v);
                    h.delete_vertex(u);
                    return UI_tree_pattern(h,print_process);
                ### Rule 2
                if g.get_vertex(v)==1 and g.get_vertex(u) in [0,-1]:
                    if print_process:
                        print v,"(1,0,0)";
                    h=g.copy();
                    h.delete_vertex(v);
                    h.set_vertex(u,-1);
                    return UI_tree_pattern(h,print_process);
                if g.get_vertex(v)==-1 and g.get_vertex(u) in [0,1]:
                    if print_process:
                        print v,"(0,1,0)";
                    h=g.copy();
                    h.delete_vertex(v);
                    h.set_vertex(u,1);
                    return UI_tree_pattern(h,print_process);
        ### In the case that g has an edge yet no rules can be applied.
        return False;

def embed_subpattern(sub_ptn, ptn):
    """
    Input:
        sub_ptn, ptn: two sign patterns;
    Output:
        find a principle subpattern of ptn that is equivalent to sub_ptn;
        return True or False;
    """
    sim_cls=sub_ptn.sign_similar_class();
    k=sub_ptn.dim[0];
    n=ptn.dim[0];
    mtx=ptn.repr;
    for com in Combinations(range(n),k):
        mtx_com=mtx[com,com];
        ptn_com=SignPattern(mtx_com);
        if ptn_com.canonical_label() in sim_cls:
            return True;
    return False; 

def unique_inertia(mtx,S_seq=None):
    """
    Input:
        mtx: a matrix representing the sign pattern;
        S_seq: the S_sequence of mtx.  It will be computed if not givne;
    Output:
        Based on the known criteria, tell if the sign pattern requires a unique inertia or not;
        Scheme implemented:
            0) Trivial cases: if S_seq is a zero sequence, mtx is zero and return True;
            1) If last_nonzero==-2, return False;
            2) If symmetric tree sign pattern, return UI_tree_pattern(tree);
            3) If skew-symmetric tree sign pattern, return True;
            4) Separate S_seq to two sequences even_seq and odd_seq;
               Use Discartes' rule of sign to determine the number of pure imaginary eigenvalues;
            5) Embedded SAP. If an SAP of order k can be embedded in a sign pattern of order n<=2k-1, then return False.            
    """
    if S_seq==None:
        S_seq=S_sequence(mtx);
    n=len(S_seq)-1;
    ptn=SignPattern(mtx);
    
    ### 0) Trivial cases
    if S_seq==[0]*(n+1):
        return True;
    
    ### 1) Check last_nonzero;
    for i in range(n,-1,-1):
        if S_seq[i]!=0:
            last_nonzero=S_seq[i];
            break;
    if last_nonzero==-2: #the product of nonzero eigenvalues can be positive or negative, so not unique inertia
        return False;
    else: #sign nonsingular
        ### 2) Symmetric tree sign patterns
        g=is_tree_pattern(mtx);
        if isinstance(g,Graph):
            return UI_tree_pattern(g);        
        
        ### 3) Skew-symmetric tree sign patterns
        if is_skew_tree_pattern(mtx):
            return True;
            
        ### 4) Descartes' rule of sign (assuming the number of zero eigenvalues is fixed).
        even_seq=[];
        odd_seq=[];
        counter=0;
        for s in S_seq:
            if counter%2==0:
                even_seq.append(s);
            if counter%2==1:
                odd_seq.append(s);
            counter+=1;
        if -2 not in even_seq and (1 in even_seq or -1 in even_seq) and sign_changes(odd_terms_negation(even_seq))==0:
            return True;
        if -2 not in odd_seq and (1 in odd_seq or -1 in odd_seq) and sign_changes(odd_terms_negation(odd_seq))==0:
            return True;           
            
        ### 5) Embed SAP. (assuming n>=3.)
        T2=SignPattern(matrix(2,[-1,-1,1,1]));
        known_SAPs={i:[] for i in range(101)};
        known_SAPs[2]=[T2];
        for k in range(integer_ceil((n+1)/2),n): #k=n is unnecessary since SAP of order n has last_nonzero=-2
            for sub_ptn in known_SAPs[k]:
                if embed_subpattern(sub_ptn,ptn):
                    return False;
        return "No conclusion yet";
