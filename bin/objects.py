class Cell:
    """ Class describing a T-cell with A, B, G and D chains.

        Attributes
        ----------
        name : str
            Name of the cell. It usually includes the Plate, tissue and well.
        n_A : int
            Number of A chains.
        n_B : int
            Number of B chains.
        n_G : int
            Number of G chains.
        n_D : int
            Number of D chains.
        batch : str
            Name of the sequencing batch.
        AB_chain_names : list
            List with strings representing the locus and alleles of AB chains.
            E.g. ['A_1','B_1']
        GD_chain_names : list
            List with strings representing the locus and alleles of AB chains.
            E.g. ['G_1','D_1']
        A_chains : list
            Placeholder. List containing the A chains of the cell. Initialized with Nones
            and reeplaced with chain object when calling add_chain().
        B_chains : list
            Placeholder. List containing the B chains of the cell. Initialized with Nones
            and reeplaced with chain object when calling add_chain().
        G_chains : list
            Placeholder. List containing the G chains of the cell. Initialized with Nones
            and reeplaced with chain object when calling add_chain().
        D_chains : list
            Placeholder. List containing the D chains of the cell. Initialized with Nones
            and reeplaced with chain object when calling add_chain().

        Methods
        -------
        __init__(name=None,n_A=0,n_B=0,n_G=0,n_D=0)
            Creates the Cell object, initializing the attribute self.name with the given
            parameter. Calls functions to create chain placeholders.
        reset_AB(n_A=0,n_B=0)
            Initializes the attributes self.n_A and self.n_B and creates A, B chain placeholder
            lists. Then, calls internal function to initialize A, B chain names.
        reset_GD(n_G=0,n_D=0)
            Initializes the attributes self.n_G and self.n_D and creates A, B chain placeholder
            lists. Then, calls internal function to initialize G, D chain names.
        add_chain(chain)
            Wrapper. Tests for chain locus value and calls the locus-specific chain-adding
            function.
        _add_A_chain(chain)
            Internal function. Checks for a correct initialization of the A chain placeholder
            and adds the A chain in the first available position. If there are not availiable
            positions raises IndexError.
        _add_B_chain(chain)
            Internal function. Checks for a correct initialization of the B chain placeholder
            and adds the B chain in the first available position. If there are not availiable
            positions raises IndexError.
        _add_G_chain(chain)
            Internal function. Checks for a correct initialization of the G chain placeholder
            and adds the G chain in the first available position. If there are not availiable
            positions raises IndexError.
        _add_D_chain(chain)
            Internal function. Checks for a correct initialization of the D chain placeholder
            and adds the D chain in the first available position. If there are not availiable
            positions raises IndexError.
        _create_AB_chain_names()
            Internal function. Create name list for the AB chains based on the (previously
            initialized) chain list.
        _create_GD_chain_names()
            Internal function. Create name list for the GD chains based on the (previously
            initialized) chain list.
        add_batch(batch)
            Initialize the attribute self.batch with the given parameter.
    """
    def __init__(self, name=None,n_A=0,n_B=0,n_G=0,n_D=0):
        self.name = name
        self.reset_AB(n_A,n_B)
        self.reset_GD(n_G,n_D)

    def reset_AB(self,n_A=0,n_B=0):
        self.n_A = n_A
        self.n_B = n_B
        self.A_chains = n_A*[None]
        self.B_chains = n_B*[None]
        self._create_AB_chain_names()
    def reset_GD(self,n_G=0,n_D=0):
        self.n_G = n_G
        self.n_D = n_D
        self.G_chains = n_G*[None]
        self.D_chains = n_D*[None]
        self._create_GD_chain_names()

    def add_chain(self,chain):
        if chain.locus=='A':
            self._add_A_chain(chain)
        elif chain.locus=='B':
            self._add_B_chain(chain)
        elif chain.locus=='G':
            self._add_G_chain(chain)
        elif chain.locus=='D':
            self._add_D_chain(chain)
        else:
            raise ValueError\
            ("The chain has a invalid locus value. Has to be in ['A','B','C','D']")

    def _add_A_chain(self,chain):
        if len(self.A_chains)==0:
            raise AssertionError\
            ("No number of A chains specified. Run reset_A with non-zero parameters")
        for i in range(len(self.A_chains)):
            if self.A_chains[i] is None:
                self.A_chains[i] = chain
                break
            else:
                if i==len(self.A_chains)-1:
                    raise IndexError('This cell does not fit more A chains')
                else:
                    continue
    def _add_B_chain(self,chain):
        if len(self.B_chains)==0:
            raise AssertionError\
            ("No number of B chains specified. Run reset_B with non-zero parameters")
        for i in range(len(self.B_chains)):
            if self.B_chains[i] is None:
                self.B_chains[i] = chain
                break
            else:
                if i==len(self.B_chains)-1:
                    raise IndexError('This cell does not fit more B chains')
                else:
                    continue
    def _add_G_chain(self,chain):
        if len(self.G_chains)==0:
            raise AssertionError\
            ("No number of G chains specified. Run reset_G with non-zero parameters")
        for i in range(len(self.G_chains)):
            if self.G_chains[i] is None:
                self.G_chains[i] = chain
                break
            else:
                if i==len(self.G_chains)-1:
                    raise IndexError('This cell does not fit more G chains')
                else:
                    continue
    def _add_D_chain(self,chain):
        if len(self.D_chains)==0:
            raise AssertionError\
            ("No number of D chains specified. Run reset_D with non-zero parameters")
        for i in range(len(self.D_chains)):
            if self.D_chains[i] is None:
                self.D_chains[i] = chain
                break
            else:
                if i==len(self.D_chains)-1:
                    raise IndexError('This cell does not fit more D chains')
                else:
                    continue
    def _create_AB_chain_names(self):
        A_names = ['A_' + str(i) for i in list(range(1,self.n_A+1))]
        B_names = ['B_' + str(i) for i in list(range(1,self.n_B+1))]
        self.AB_chain_names = A_names + B_names

    def _create_GD_chain_names(self):
        G_names = ['G_' + str(i) for i in list(range(1,self.n_G+1))]
        D_names = ['D_' + str(i) for i in list(range(1,self.n_D+1))]
        self.GD_chain_names = G_names + D_names

    def add_batch(self,batch):
        self.batch = batch
#--------------------------------------------------------------------------------------------------#
class Chain:
    """ Abstract class for TCR chain holding all relevant characteristics.

        Attributes
        ----------
        locus : str
            Locus of chain ('A','B','G','D').
        allele : str
            A string composed of the locus followed by an underscore '_' and the number of the
            chain allele. E.g. 'A_1', 'G_2'.
        productive : str
            String representing a boolean. Wether the chain is productive or not.
        TPM : float
            Transcripts per million.
        stop_codon : str
            String representing a boolean. Wether the chain has a stop codon or not.
        in_frame : str
            String representing a boolean. ???.
        ID : str
            String composed of the V usage, a fraction of the CDR3nt chain and the J usage, all
            separated by underscores '_'.
        CDR3nt : str
            Nucleotide sequence of the TCR
        CDR3aa : str
            Amino acid sequence of the TCR.
        V : str
            V usage.
        J : str
            J usage.

        Methods
        -------
        __init__(locus,allele)
            Creates the Chain object, initializing the attributes self.locus and self.allele
            with the given parameters and the rest with None.
        fill_metadata(meta_dict)
            Reads a python dictionary and assigns its values to the class attributes.
    """
    def __init__(self,locus,allele):
        self.locus = locus
        self.allele = allele
        self.productive = None
        self.TPM = None
        self.stop_codon = None
        self.in_frame = None
        self.ID = None
        self.CDR3nt = None
        self.CDR3aa = None
        self.V = None
        self.J = None

    def fill_metadata(self,meta_dict):
        self.productive = meta_dict['Productive']
        self.V = meta_dict['V segment']
        self.J = meta_dict['J segment']
        self.CDR3aa = meta_dict['CDR3aa']
        self.CDR3nt = meta_dict['CDR3nt']
        self.TPM = float(meta_dict['TPM'])
        self.stop_codon = meta_dict['Stop codon']
        self.in_frame = meta_dict['In frame']
        self.ID = meta_dict['ID']
#-----------------------------------------------------------------------------------------------------
class AlphaChain(Chain):
    """ Wrapper of the abstract class with no modifications."""
    def __init(self):
        super().__init__()
class BetaChain(Chain):
    """ Wrapper of the abstract class that includes the D usage attribute and function."""
    def __init(self):
        super().__init__()
    def add_D_segment(self,D):
        self.D = D
class GammaChain(Chain):
    """ Wrapper of the abstract class with no modifications."""
    def __init(self):
        super().__init__()
class DeltaChain(Chain):
    """ Wrapper of the abstract class with no modifications."""
    def __init(self):
        super().__init__()
    def add_D_segment(self,D):
        self.D = D
#-----------------------------------------------------------------------------------------------------
def create_cell_from_AB(in_file):
    """ Reads the AB output from TraCeR assemble for a given cell and returns the Cell object.

      Reads the 'filtered_TCRs.txt' file coming from TraCeR assemble for a given cell,
      extracts the information for the available chains and writes it in a new Cell object
      as metadata.

      Parameters
      ----------
      in_file : string
        Path for the AB 'filtered_TCRs.txt' file.

      Raises
      ------
      ValueError
          If the locus of a chain is not 'A' or 'B'. This can only happen if the cell is not
          correctly initialized.

      Returns
      -------
      Cell
        Cell object with all the chains and their meta-data loaded.
    """
    # Read file
    with open(in_file,'r',encoding='utf8') as f:
        lines = f.readlines()
    # Create Cell object
    n_A = int(lines[3].strip()[-1])
    n_B = int(lines[4].strip()[-1])
    cell = Cell(lines[1].strip(),
                     n_A = n_A,
                     n_B = n_B)

    cont = 0 # Counter for the number of chains
    for i in range(len(lines)):
        if lines[i].startswith('##TRINITY'): # Detect chain data
            # Extract data from text chunk
            text_chunk = lines[i+1:i+11]
            i=i+10
            chain_dict = {}
            allele = cell.AB_chain_names[cont]
            locus = allele[0]
            for line in text_chunk:
                # Organize data fields in dictionary
                data_pair = line.strip().split(':\t')
                if len(data_pair)>1: # Exclude blank lines
                    chain_dict[data_pair[0]]=data_pair[1]
            # Create Chain object based on locus
            if locus=='A':
                chain = AlphaChain(locus,allele)
            elif locus=='B':
                chain = BetaChain(locus,allele)
                chain.add_D_segment(chain_dict['D segment'])
            else:
                raise ValueError("The locus has to be either 'A' or 'B'")

            chain.fill_metadata(chain_dict) # Fill in data from file
            cell.add_chain(chain) # Add chain to cell
            cont = cont+1 # Next chain
    return cell
#-----------------------------------------------------------------------------------------------------
def append_GD_data(in_file,cells):
    """ Reads the GD output from TraCeR assemble apends data to the right Cell objects.

      Reads the 'filtered_TCRs.txt' file coming from TraCeR assemble for a given cell,
      extracts the information for the available chains and writes it in an existing Cell
      object, retrieved from the cells dictionary.

      Parameters
      ----------
      in_file : string
        Path for the GD 'filtered_TCRs.txt' file.
      cells : dictionary
        Dictionary where the cells initialized with AB were stored.

      Raises
      ------
      ValueError
          If the locus of a chain is not 'G' or 'D'. This can only happen if the cell is not
          correctly initialized.
    """
    # Read file
    with open(in_file,'r',encoding='utf8') as f:
        lines = f.readlines()
    # Look for initialized Cell object
    n_G = int(lines[3].strip()[-1])
    n_D = int(lines[4].strip()[-1])
    name = lines[1].strip()
    cell = cells[name]
    cell.reset_GD(n_G,n_D)

    cont = 0 # Counter for the number of chains
    for i in range(len(lines)):
        if lines[i].startswith('##TRINITY'): # Detect chain data
            # Extract data from text chunk
            text_chunk = lines[i+1:i+11]
            i=i+10
            chain_dict = {}
            allele = cell.GD_chain_names[cont]
            locus = allele[0]
            for line in text_chunk:
                # Organize data fields in dictionary
                data_pair = line.strip().split(':\t')
                if len(data_pair)>1: # Exclude blank lines
                    chain_dict[data_pair[0]]=data_pair[1]
            # Create Chain object based on locus
            if locus=='G':
                chain = GammaChain(locus,allele)
            elif locus=='D':
                chain = DeltaChain(locus,allele)
                chain.add_D_segment(chain_dict['D segment'])
            else:
                raise ValueError("The locus has to be either 'G' or 'D'")

            chain.fill_metadata(chain_dict) # Fill in data from file
            cell.add_chain(chain) # Add chain to cell
            cont = cont+1 # Next chain
    cells[name] = cell # Update cell
