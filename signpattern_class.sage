class SignPattern:
    def __init__(self,mtx):
        """
        Use SignPattern(mtx) to construct a sign pattern.
        When mtx is a digraph with arcs labelled by "+" or "-",
        translate it into a matrix representation.
        Use self.repr to call the matrix representation.
        """
        self.digraph_repr=None;
        if isinstance(mtx,sage.graphs.digraph.DiGraph):
            self.digraph_repr=mtx;
            n=mtx.order();
            m=n; # for defining self.dim later;
            ptn=zero_matrix(n);
            for i in range(n):
                for j in range(n):
                    try:
                        lbl=mtx.edge_label(i,j);
                    except LookupError:
                        lbl=0;
                    if lbl=="+":
                        ptn[i,j]=1;
                    if lbl=="-":
                        ptn[i,j]=-1;
        else:
            ptn=copy(mtx);
            m,n=ptn.dimensions();
            for i in range(m):
                for j in range(n):
                    if mtx[i,j]>0:
                        ptn[i,j]=1;
                    if mtx[i,j]<0:
                        ptn[i,j]=-1;
        self.repr=ptn;
        self.dim=m,n;
        
    def __repr__(self):
        return "%s by %s sign pattern"%self.dim;
        
    def __eq__(self,ptn_B):
        return self.repr==ptn_B.repr;

    def show(self):
        show(self.repr);
        
    def digraph(self):
        m,n=self.dim;
        assert m==n,"Only square sign pattern has digraph representation.";
        D=DiGraph(n,loops=True);
        ptn=self.repr;
        for i in range(m):
            for j in range(n):
                if ptn[i,j]>0:
                    D.add_edge(i,j,"+");
                if ptn[i,j]<0:
                    D.add_edge(i,j,"-");
        self.digraph_repr=D;
        return D;            
        
    def is_isomorphic(self,ptn_B):
        DA=self.digraph().canonical_label(edge_labels=True);
        DB=ptn_B.digraph().canonical_label(edge_labels=True);
        return DA==DB;

    def is_irreducible(self):
        if self.digraph_repr==None:
            self.digraph_repr=self.digraph();
        return self.digraph_repr.is_strongly_connected();

    def negation(self):
        return SignPattern(-self.repr);

    def transpose(self):
        return SignPattern(self.repr.transpose());

    def sign_similar(self,diag):
        """
        Input:
            diag: a list containing 1 and -1;
        Output:
            return the pattern obtained by diag*self*diag;
        """
        D=diagonal_matrix(diag);
        return SignPattern(D*self.repr*D);

    def canonical_label(self):
        return SignPattern(self.digraph().canonical_label(edge_labels=True));

    def sign_similar_class(self,trans=True,neg=True):
        m,n=self.dim;
        assert m==n,"Only square sign pattern has digraph representation.";
        all_canonical_ptns=[];
        for k in range(integer_floor(n/2)+1):
            for com in Combinations(range(n),k):
                diag=[1]*n;
                for i in com:
                    diag[i]=-1;
                new_ptn=self.sign_similar(diag).canonical_label();
                #print "diag",diag;
                #print new_ptn.repr;
                if new_ptn not in all_canonical_ptns:
                    all_canonical_ptns.append(new_ptn);
        if trans:
            #print "transpose starts";
            for new_ptn in all_canonical_ptns:
                trans_ptn=new_ptn.transpose().canonical_label();
                if trans_ptn not in all_canonical_ptns:
                    all_canonical_ptns.append(trans_ptn);
        if neg:
            #print "negation starts";
            for new_ptn in all_canonical_ptns:
                neg_ptn=new_ptn.negation().canonical_label();
                if neg_ptn not in all_canonical_ptns:
                    all_canonical_ptns.append(neg_ptn);
        return all_canonical_ptns;   
        
    def is_equivalent(self,ptn_B,trans=True,neg=True):         
        return (ptn_B.canonical_label() in self.sign_similar_class(trans,neg));
