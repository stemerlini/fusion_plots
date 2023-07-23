import numpy as np
from urllib.request import urlopen

def read_NIST_data( url='', fname='' ):
#{{{
    """
    Read the atomic weight dataset from NIST from web or file.

    Parameters
    ----------
    url: string
        datasets from NIST can be obtained via queries, here we need the following
        http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some
    fname: string
        alternatively to directly reading the data from NIST, NIST also offers to 
        download an ascii file (which you need to do before running this script)
        
    Returns
    -------
    NIST_dataset: dict
    """

    # dictionary into which the data will be saved
    # note: keys must exactly correspond to identifiers on website / in file
    #       more clever way would be to create the dictionary dynamically
    #       reading (using) the identifiers from the website
    NIST_dataset = { 'Atomic Number':[],
                     'Atomic Symbol':[],
                     'Mass Number':[],
                     'Relative Atomic Mass':[]
                   }

    # if url is provided, read NIST dataset from webpage
    if len(url) > 0:
        # open connection to URL and check if everything is fine
        web_NIST = urlopen(url)
        if (web_NIST.getcode() != 200):
            print( 'ERROR: http-code while trying to read NIST data from web: {0}'.format(web_NIST.getcode()) )
            return -1
        # read website line by line
        for line in web_NIST:
            # decode line into readable string
            decoded_line = line.decode( "utf-8" )

            for label in NIST_dataset:
                if decoded_line.startswith( label ):

                    # get index of "=", starting to look after identifier text at beginning
                    # note: python gets the correct position (text at beginning is just ignored)
                    id_0 = decoded_line.index( "=", len(label) ) + 1

                    # extract value and if present, remove error-margin, given within brackets ()
                    if '(' in decoded_line[ id_0: ]:
                        id_1    = decoded_line.index( "(", len(label) )
                    else:
                        id_1    = len(decoded_line)
                    if label == 'Atomic Symbol':
                        val = decoded_line[id_0:id_1]
                    else:
                        val = float(decoded_line[id_0:id_1])

                    # add value to dictionary
                    NIST_dataset[ label ].append(val)
         

    # if filename is provided during function call, read NIST dataset from file
    if len(fname) > 0:
        # loop through file line by line
        for line in open( fname ):
            if "=" in line:
                tag, value = line.split( '=' )
                tag   = tag.strip()
                value = value.strip()
                # remove parantheses in values which indicate error
                if '(' in value:
                    value = value.split( '(' )[0]
                if tag in NIST_dataset.keys():
                    NIST_dataset[ tag ].append( value )

            # check if list lengths is the same for each key
            if line == '\n':
                lengths = []
                # get the lengths of the arrays for each key
                for key in NIST_dataset:
                    lengths.append( len( NIST_dataset[key] ) )
                # check if all the lengths are equal
                if not all( elem == lengths[0] for elem in lengths ):
                    print( 'WARNING: there is an error in extracting the data for the file' )
                    print( '         and sorting it into a dictionary' )

    return NIST_dataset
#}}}


def get_mass_number( NIST_dataset ):
#{{{
    """
    Extracts the mass number from the NIST dataset.

    Parameters
    ----------
    NIST_dataset: dict

    Returns
    -------
    mass_number: numpy array
    """

    atomic_mass = np.asarray( NIST_dataset['Relative Atomic Mass'], dtype=np.float32 )
    mass_number = np.around( atomic_mass )

    return mass_number
#}}}


def get_binding_energy( NIST_dataset, norm=True):
#{{{
    """

    Parameters
    ----------
    NIST_dataset: dict
    norm: boolean

    Returns
    -------
    binding_energy: numpy array
    """

    # atomic mass unit in units of MeV/c^2
    amu         = 931.494095
    # masses of particle in units of u
    neutron_u   = 1.00866491588
    proton_u    = 1.00727646688
    electron_u  = 0.00054858

    atomic_mass     = np.asarray( NIST_dataset['Relative Atomic Mass'], dtype=np.float32 )
    atomic_number   = np.asarray( NIST_dataset['Atomic Number'], dtype=np.float32 )
    mass_number     = get_mass_number( NIST_dataset )
    n_neutrons      = mass_number - atomic_number
    mass_defect     = ( ( n_neutrons*neutron_u + atomic_number*proton_u + atomic_number*electron_u ) - atomic_mass ) * amu

    if norm:
        bind_energy = mass_defect/mass_number
    else:
        bind_energy = mass_defect

    return bind_energy
#}}}